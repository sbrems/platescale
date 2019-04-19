import os
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# from subprocess import call
# import skimage.feature
# from scipy import ndimage
# from scipy import signal
# from scipy.stats import sigmaclip
import numpy as np
# from astropy.io import fits
import astrom
from astrom import misc, analysis, sextract,\
    psf_fitting, calc_platescale, plots  # , parameters
from astrom.parameters import sigma_outliers, conv_gauss, keepfr, dir_cat,\
    keepfr_med, true_north_guess, pxscl_init


pd.options.mode.chained_assignment = None  # default='warn'


def do(verbose=True, plot=True, delete_old=True,
       use_mags=True, dir_data2=None):
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
    if dir_data2 is not None:
        dir_data = dir_data2
        dir_out = os.path.join(dir_data, 'results')
        dir_temp = os.path.join(dir_out, 'temp')
    else:
        dir_data = astrom.parameters.dir_data
        dir_out = astrom.parameters.dir_out
        dir_temp = astrom.parameters.dir_temp
    dir_default = os.getcwd()
    astrom.make_dirs([dir_out, dir_temp],
                     delete_old=delete_old, verbose=verbose)
    filename_cube, data_cube, header_cube = astrom.read_fits(
        dir_data, verbose=verbose)
    n_images = data_cube.shape[0]
    if n_images == 1:
        print('Only one image found. Skipping analysis of single frames')
    else:
        print('{} images found. Processing median and single frame analysis'.format(
            n_images))
    # determine some parameters from the fits files
    target, mjd_obs, agpm, coord_head, rot_header = misc.get_headerparams(
        header_cube[0], verbose=verbose)
    # determine more (standard) parameters as filenames (if None in parameters file)
    mjd_source, fn_sextr_med, fn_sextr_sgl, fn_source = misc.get_catfiles(target, agpm,
                                                                          dir_cat, dir_temp,
                                                                          verbose=verbose)

    data_med = np.median(data_cube, axis=0)

    ###################################################
    ####get the coordinates of the objects found ######
    ####in the median and create the objects.fits######
    ####file for the artificial map####################
    median_sources, data_objects_med = sextract.coordinates(data_med, header_cube, target,
                                                            fn_sextr_sgl, fn_sextr_med,
                                                            dir_data, dir_out, dir_temp,
                                                            med=True, verbose=verbose)

    ###################################################
    ####create the source list#########################
    ###################################################

    # quick sanity check
    if header_cube[-1]['MJD-OBS'] - mjd_obs > 1:
        raise ValueError('Observations not from the same day:',
                         mjd_obs, header_cube[-1]['MJD-OBS'])

    source_cat = misc.get_catalogue_data(fn_source, mjd_source, coord_head,
                                         dir_out,
                                         imsize=np.max(data_med.shape)*pxscl_init/1000,
                                         mjd_obs=mjd_obs,
                                         verbose=verbose, plot=plot)

    ###################################################
    #####create an artificial map with the sources in##
    #####and match the sources with the ones found by###
    #####sextractor####################################
    ###################################################

    source_cat = plots.make_artificial_map(source_cat, data_objects_med,
                                           rot_header,
                                           verbose=verbose, plot=plot,
                                           use_mags=use_mags,
                                           dir_out=dir_out)

    
    # Only use the sources found in the median image further. Find them
    source_med, ignored_med = analysis.connect_median_sources(source_cat,
                                                              median_sources,
                                                              verbose=verbose)
    if len(source_med) <= 1:
        raise TooFewStarsFoundError
    # plot the median image with all the sources found and connected
    if plot:
        plots.median_stars(data_med, source_med, ignored_med, source_cat,
                           dir_out=dir_out)
    ############################################
    ####sextract the sources####################
    ############################################

    sex_coords = sextract.coordinates(data_cube, header_cube, target, fn_sextr_sgl, fn_sextr_med,
                                      dir_data, dir_out, dir_temp,
                                      med=False, verbose=verbose)
    psf_med = psf_fitting.make_psf(
        np.array([data_med]), [source_med], dir_temp,
        fn_out='med_cube', keepfr=keepfr_med)
    source_med = psf_fitting.refine_fit(np.array([data_med]), psf_med, [source_med], ign_ims=[],
                                        conv_gauss=conv_gauss, verbose=verbose, plot=plot)[0]
    # excluded is excluded full images pd frame. Ignored is a list with len of n_images containing
    # single, ignored sources.i_ign im is a list whith the numbers of the ignored frames as they
    # are f.e. in data_cube. Do this only if there is more than one image in the cube.
    if n_images == 1:
        excluded = []
        ignored = []
        nr_ign_im = []
        psf_stars = psf_med
        sex_coords = [source_med, ]
    else:
        sex_coords, excluded, ignored, nr_ign_im = analysis.compare_to_med(sex_coords,
                                                                           source_med,
                                                                           ignored_med,
                                                                           verbose=verbose)

        # make a reference psf to use to center the stars
        psf_stars = psf_fitting.make_psf(data_cube, sex_coords, dir_temp,
                                         ign_ims=nr_ign_im, n_rms=2,
                                         conv_gauss=conv_gauss,
                                         keepfr=keepfr, fn_out='psf_cube')
        # find all targets via a ccmap and the psf just created)
        # sex_coords_ccmap = psf_fitting.find_via_ccmap(data_cube,psf_stars,target, dir_data, dir_temp)
        # refine the stellar positions with the new psf. old ones in x/y_image_sex column
        sex_coords = psf_fitting.refine_fit(data_cube, psf_stars, sex_coords, ign_ims=nr_ign_im,
                                            conv_gauss=conv_gauss, verbose=verbose)

    ############################################
    #now do the astrometry######################
    ############################################
    # connect the sextracted coordinates to the ones found in the catalogue. Also make some
    # sanity checks, e.g. more sources were found than in catalogue. remove those and save in
    # excluded
    # sex_coords, excluded, multiples = analysis.connect_sources(sex_coords,source_med,
    #                                                                   n_images=n_images)

    #########################################################
    #do the statistics. e.g. measure all the distances#######
    #########################################################
 
    os.chdir(dir_out)
    # in astrometry the single distances/angles between all possible combinations are saved.
    # main_id1 is the one coming first in the alphabet
    # in results these will be grouped.
    astrometry = calc_platescale.single_matches(data_cube, sex_coords, nr_ign_im,
                                                rot_header,
                                                verbose=verbose, plot=plot)
    # if only stars with weight 0 found (e.g. they were too close to the border), ignore them
    if np.max(astrometry['weight']) <= 0.:
        raise TooFewStarsFoundError
    ################################################################
    #now combine all measurements that have the same start and end##
    ################################################################

    astrometry_grouped, astrometry_grouped_cliped, results = calc_platescale.multiple_matches(
        astrometry, sigma_outliers, )
    # add some parameters to be saved
    results['target'] = target
    results['agpm'] = agpm
    results['mjd_obs'] = mjd_obs
    results['nr_used_med_stars'] = np.sum(np.isfinite(
        source_med['shift_error']))
    results['nr_connections_tot'] = len(np.where(np.array(astrometry['weight'])
                                                 > 0)[0])
    results['rot_header'] = rot_header
    results['true_north_guess'] = true_north_guess

    if plot and n_images > 1:
        # make some plots about the statistics
        plots.distance_distributions(astrometry, astrometry_grouped,
                                     dir_out=dir_out)
    # now ignore the star completely where the fraction of bad scales was higher than ignore_frac
    astrometry_good = calc_platescale.ignore_stars_bad_connections(astrometry_grouped,
                                                                   astrometry_grouped_cliped,
                                                                   verbose=verbose)
    # pd.DataFrame([len(results_entries)*[np.nan]],columns=results_entries)
    # stars_found = list(set(astrometry_grouped[star_id+'1']+astrometry_grouped[star_id+'2']))

    print('The final pixel_scale and true north are:\n', results)
    os.chdir(dir_default)

    # save the results
    astrometry.to_csv(os.path.join(dir_out, 'astrometry.csv'), index=False)
    results.to_csv(os.path.join(dir_out, 'results.csv'), index=False)
    print('Saved results to {}'.format(os.path.join(dir_out, 'results.csv')))

    if False:
        # compare to andres results
        andre = np.array([[614.463, 548.146, 556.571, 527.867, 412.971], [
            565.113, 426.196, 578.740, 468.201, 588.342]])
        # only for trap_med
        sex_coords[0]['x_andre'] = andre[0]
        sex_coords[0]['y_andre'] = andre[1]
        temp = []
        for first in range(len(sex_coords[0])):
            for second in range(first + 1, len(sex_coords[0])):
                dist_meas, ang_meas = misc.dist_ang((sex_coords[0]['y_andre'][first],
                                                     sex_coords[0]['x_andre'][first]),
                                                    (sex_coords[0]['y_andre'][second],
                                                     sex_coords[0]['x_andre'][second]),
                                                    spherical=False)
                temp.append(dist_meas)
        astrometry['dis_andre'] = temp

        dis_an = []
        ang_an = []
        for first in range(andre.shape[1]):
            for second in range(first + 1, andre.shape[1]):
                dis_an.append(misc.dist_ang((andre[1][first], andre[0][first]),
                                            (andre[1][second], andre[0][second]))[0])
                ang_an.append(misc.dist_ang((andre[1][first], andre[0][first]),
                                            (andre[1][second], andre[0][second]))[1])
        astrometry['ps_andre'] = astrometry['dis_cat_mas'] / \
            astrometry['dis_andre']
        plt.imshow(np.sqrt(data_med), interpolation='none')

        plt.savefig(os.path.join(dir_temp, 'found_sources_1.svg'))
    print('Finished platescale.py successfully!')

    return results


class TooFewStarsFoundError(Exception):
    '''Exception if too few stars were found in the image'''
    pass
