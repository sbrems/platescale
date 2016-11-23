from __future__ import print_function,division
import os
import numpy as np
import pandas as pd
import math
from astropy.io import fits
from subprocess import call
import subprocess
#from scipy.signal import convolve2d
from parameters import *


def coordinates(data_cube,header_cube,target,fn_sextr_sgl,fn_sextr_med,med=False,verbose=True):
    '''###################################################
    ####get sextractor to get all the coordinates#######
    ####################################################
    Run sextractor on each single image to get the coordinates of the sources
    If there was a source detected in the median image, mark save the coordinates
    This is done very inefficiently, but we have time :)
    If med=True, doing the same with different output names 
    (_med instead of _{image_number}) and returning also the med image'''

    os.chdir(dir_out)
    FNULL = open(os.devnull,'w') #suppress output
    if not med:
        n_images = data_cube.shape[0]
        if verbose: print('Sextracting from ',n_images,' images.')
        for kk in range(n_images):
            fits.writeto('temp.fits',data_cube[kk,:,:],output_verify = 'ignore')
            call(['sex','temp.fits','-c',fn_sextr_sgl,
                  '-checkimage_name',dir_temp+target+'_objects_'+'{:03}'.format(kk)+'.fits',
                  '-catalog_name',dir_temp+target+'_'+'{:03}'.format(kk)+'.sex'],stdout=FNULL,
                 stderr=subprocess.STDOUT)
            call(['rm','temp.fits'])
        os.chdir(dir_temp)
        all_objects = np.full(data_cube.shape, np.nan)
        sex_coords = [None] * n_images
        for kk in range(n_images):
            all_objects[kk,:,:] = fits.getdata(target+'_objects_'+'{:03}'.format(kk)+'.fits')
            sex_coords[kk] =pd.read_csv(target+'_'+'{:03}'.format(kk)+'.sex', header=None,
                                        delimiter=r'\s+',comment='#',
                                        names=('mag_auto','x_image','y_image','rms_A_image','rms_B_image'),
                                        dtype={'x_image':np.float64,'y_image':np.float64},
                                        skip_blank_lines = True)
            #subtract 1 as sextractor starts at (1,1) and not (0,0)
            sex_coords[kk]['x_image'] -= 1
            sex_coords[kk]['y_image'] -= 1
        fits.writeto(dir_out+target+'_objects_cube.fits',all_objects,header_cube[0],output_verify = 'ignore')

    if med:
        if verbose: print('Processing the median image.')
        fits.writeto('temp.fits',data_cube,output_verify = 'ignore')
        call(['sex','temp.fits','-c',fn_sextr_med,
              '-checkimage_name',dir_temp+target+'_objects_med.fits',
              '-catalog_name',dir_temp+target+'_med.sex'],stdout=FNULL,
             stderr=subprocess.STDOUT)
        call(['rm','temp.fits'])
        #now that we have all data in files, read it in
        os.chdir(dir_temp)
        all_objects = np.full(data_cube.shape, np.nan)
        sex_coords = [None]
        all_objects = fits.getdata(dir_temp+target+'_objects_med.fits')
        sex_coords =pd.read_csv(dir_temp+target+'_med.sex', header=None,
                                delimiter=r'\s+',comment='#',
                                names=('mag_auto','x_image','y_image','rms_A_image','rms_B_image'),
                                dtype={'x_image':np.float64,'y_image':np.float64},
                                skip_blank_lines = True)
        #subtract 1 as sextractor starts at (1,1) and not (0,0)
        sex_coords['x_image'] -= 1
        sex_coords['y_image'] -= 1
        #overwrite the just created fits file with the header
        fits.writeto(dir_out+target+'_objects_med.fits',all_objects,header_cube[0],output_verify = 'ignore')
    os.chdir(dir_data)
    if med:
        return sex_coords,all_objects
    else:
        return sex_coords
    
