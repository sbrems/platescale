from __future__ import print_function
from __future__ import division
import os
import astrom
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from subprocess import call
#import skimage.feature
from scipy import ndimage
#from scipy import signal
from scipy.stats import sigmaclip
import astrom.plots
import astrom.analysis
from astrom import misc
from astrom import sextract
from astrom import psf_fitting
from astrom import calc_platescale
from collections import Counter
import numpy as np
from astropy.io import fits
from astrom.parameters import *

import ipdb

pd.options.mode.chained_assignment = None  # default='warn'

def do():
    '''This code gives you the platescale of imgaes. It needs a file parameters.py which contains a number of parameters,
    of which the most important are the catalogue with the true star positions (and pm if available), aswell as the 
    initial guess of the platescale and the True nort.
    The idea is to take the median of the (rotated) images and sextract the sources in this file. This is then used to
    crosscorrelate with the input catalogue, where here at each position a gauss-like source is put. In order to do so,
    we need an initial platescale and true north to make the images of the same shape (the routine cannot scale or rotate).
    The found position is then used to sextract all stars in the single images and connect them to a star catalogue. You
    can check wether this worked in the output images.
    Once this is done, all possible connections between the sources in each image are done, and we compare the measured
    distances and angles aswell as the ones from the catalogue. This gives us the statistics for the output.
    The important pandas-tables are:
    source: containing the coordinates of the input table modified with date and initial platescale/true north
    sex_coords[im_number]: contains the Identified star and the px-coordinates and ra/dec
    astrometry: gives all the combinations between the different targets and their distance and angle measured and from the catalogue. len= (n_images * n_stars*(n_stars+1)/2 )
    astrometry_grouped: grouped by the different combinations
    results: gives the final platescale    
    '''
    dir_default = os.getcwd()
    astrom.make_dirs([dir_out,dir_temp],delete_old=True)
    filename_cube, data_cube, header_cube =  astrom.read_fits(dir_data)

    data_med = np.median(data_cube,axis=0)
    

    #get the objects via sextractor. The output file is called objects.fits
    os.chdir(dir_out)
    if agpm:
        fn_median ='median_'+target+'_agpm.fits'
    else:
        fn_median ='median_'+target+'.fits'
    fits.writeto(fn_median,data_med,header_cube[0],output_verify = 'ignore')
    call(['sex',fn_median,'-c',dir_data+'default.sex','-checkimage_name',target+'_median_objects.fits'])
    data_objects_med = fits.getdata(target+'_median_objects.fits')
    os.chdir(dir_default)
    ###################################################
    ####create the source list#########################
    ###################################################

    mjd_obs = header_cube[0]['MJD-OBS']
    print(mjd_obs)
    #quick sanity check
    if header_cube[-1]['MJD-OBS'] - mjd_obs > 1:
        raise ValueError('Observations not from the same day:',mjd_obs, header_cube[-1]['MJD-OBS'])
    source  = astrom.misc.get_catalogue_data(fn_source,mjd_obs=mjd_obs)

    ###################################################
    #####create an artificial map with the sources in##
    #####and match the souces with the ones found by###
    #####sectractor####################################
    ###################################################
    source  = astrom.plots.make_artificial_map(source,data_objects_med)

    ############################################
    ####sextract the sources####################
    ############################################
    n_images = data_cube.shape[0]
    sex_coords  = sextract.coordinates(data_cube,header_cube,target)
    
    #make a reference psf to use to center the stars
    psf_stars = psf_fitting.make_psf(data_cube,sex_coords, n_rms = 2, conv_gauss=conv_gauss, keepfr = keepfr)
    #find all targets via a ccmap and the psf just created)
    #sex_coords_ccmap = psf_fitting.find_via_ccmap(data_cube,psf_stars)
    #now refine the stellar positions with the new psf. old ones in x/y_image_sex column
    sex_coords = psf_fitting.refine_fit(data_cube,psf_stars,sex_coords,conv_gauss=conv_gauss)
    
    ############################################
    #now do the astrometry######################
    ############################################
    #connect the sextracted coordinates to the ones found in the catalogue. Also make some
    #sanity checks, e.g. more sources were found than in catalogue. remove those and save in
    #excluded
    sex_coords, excluded, multiples = astrom.analysis.connect_sources(sex_coords,source,
                                                                      n_images=n_images)

    #########################################################
    #do the statistics. e.g. measure all the distances#######
    #########################################################

    os.chdir(dir_out)
    #in astrometry the single distances/angles between all possible combinations are saved.
    #in results these will be grouped
    astrometry = astrom.calc_platescale.single_matches(data_cube,sex_coords,
                                             n_images=n_images)

    ################################################################
    #now combine all measurements that have the same start and end##
    ################################################################
    
    astrometry_grouped, astrometry_grouped_cliped, results = astrom.calc_platescale.multiple_matches(astrometry,sigma_outliers)


    #now ignore the star completely where the fraction of bad scales was higher than ignore_frac            pd.DataFrame([len(results_entries)*[np.nan]],columns=results_entries)                                      
    #stars_found = list(set(astrometry_grouped[star_id+'1']+astrometry_grouped[star_id+'2']))                                                    
    connections_tot = Counter(dict(astrometry_grouped.groupby(star_id+'1').size())) +\
                      Counter(dict(astrometry_grouped.groupby(star_id+'2').size()))
    connections_good = Counter(dict(astrometry_grouped_cliped.groupby(star_id+'1').size())) +\
                       Counter(dict(astrometry_grouped_cliped.groupby(star_id+'2').size()))
    astrometry_good = astrometry_grouped
    stars_removed = []
    for star in connections_tot.keys():
        if connections_good[star]/connections_tot[star] < ignore_frac:
            print('Removing all measurements of star ',star, ' as it only had ', connections_good[star],\
                  ' of ',connections_tot[star],' good measurements. This is worse than ignore_frac=',ignore_frac )
            stars_removed.append(star)
            astrometry_good = astrometry_good[astrometry_good[star_id+'1'] != star]
            astrometry_good = astrometry_good[astrometry_good[star_id+'2'] != star]


    
    print('The final pixel_scale and true north are:\n',results)
    os.chdir(dir_out)
    #save the results
    astrometry.to_csv('astrometry.csv',index=False)
    results.to_csv('results.csv',index=False)

        
    #compare to andres results
    andre=np.array([[614.463,548.146,556.571,527.867,412.971],[565.113,426.196,578.740,468.201,588.342]])
    
    if False:
        #only for trap_med
        sex_coords[0]['x_andre'] = andre[0]
        sex_coords[0]['y_andre'] = andre[1]
        temp = []
        for first in range(len(sex_coords[0])):
            for second in range(first+1,len(sex_coords[0])):
                dist_meas,ang_meas= misc.dist_ang((sex_coords[0]['y_andre'][first],
                                                   sex_coords[0]['x_andre'][first]),
                                                  (sex_coords[0]['y_andre'][second],
                                                   sex_coords[0]['x_andre'][second]),
                                                  spherical=False)
                temp.append(dist_meas)
        astrometry['dis_andre'] = temp

        dis_an = []
        ang_an  = []
        for first in range(andre.shape[1]):
            for second in range(first+1,andre.shape[1]):
                dis_an.append(misc.dist_ang((andre[1][first],andre[0][first]), (andre[1][second],andre[0][second]))[0])
                ang_an.append(misc.dist_ang((andre[1][first],andre[0][first]), (andre[1][second],andre[0][second]))[1])
        astrometry['ps_andre'] =astrometry['dis_cat_mas'] /  astrometry['dis_andre']     
        #print(ii,comb,np.mean(astrometry.ps_andre[blubb]),np.std(astrometry.ps_andre[blubb]),np.median(astrometry.ps_andre[blubb]))
        #plot the median image
        plt.imshow(np.sqrt(data_med),interpolation='none')

        plt.savefig(dir_temp+'found_sources_1.svg')
    print('Finished platescale.py successfully!')
    ipdb.set_trace()