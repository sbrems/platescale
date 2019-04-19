import os
import numpy as np
import pandas as pd
# import math
from astropy.io import fits
from subprocess import call
import subprocess
# from scipy.signal import convolve2d
from .parameters import cut_size, center_sextr


def coordinates(data_cube, header_cube, target, fn_sextr_sgl, fn_sextr_med,
                dir_data, dir_out, dir_temp,
                med=False, verbose=True):
    '''###################################################
    ####get sextractor to get all the coordinates#######
    ####################################################
    Run sextractor on each single image to get the coordinates of the sources
    If there was a source detected in the median image, mark save the
    coordinates. This is done very inefficiently, but we have time :)
    If med=True, doing the same with different output names
    (_med instead of _{image_number}) and returning also the med image'''

    os.chdir(dir_temp)
    FNULL = open(os.devnull, 'w')  # suppress output
    if not med:
        n_images, y_size, x_size = data_cube.shape
        if verbose:
            print('Sextracting from ', n_images, ' images.')
        for kk in range(n_images):
            fits.writeto(
                'temp.fits', data_cube[kk, :, :], output_verify='ignore')
            call(['sex', 'temp.fits', '-c', fn_sextr_sgl,
                  '-checkimage_name', os.path.join(dir_temp,
                                                   target + '_objects_{:03}.fits'.format(kk)),
                  '-catalog_name', os.path.join(dir_temp, target + '_{:03}.sex'.format(kk))],
                 stdout=FNULL, stderr=subprocess.STDOUT)
            call(['rm', 'temp.fits'])
        all_objects = np.full(data_cube.shape, np.nan)
        sex_coords = [None] * n_images
        for kk in range(n_images):
            all_objects[kk, :, :] = fits.getdata(
                target + '_objects_{:03}.fits'.format(kk))
            if center_sextr == 'peak':
                sex_coords[kk] = pd.read_csv(target + '_{:03}.sex'.format(kk), header=None,
                                             delimiter=r'\s+', comment='#',
                                             names=('mag_auto', 'x_image', 'y_image',
                                                    'x_auto_not_used', 'y_auto_not_used',
                                                    'rms_A_image', 'rms_B_image'),
                                             dtype={'x_image': np.float64,
                                                    'y_image': np.float64},
                                             skip_blank_lines=True)
            elif center_sextr == 'auto':
                sex_coords[kk] = pd.read_csv(target + '_' + '{:03}'.format(kk) + '.sex', header=None,
                                             delimiter=r'\s+', comment='#',
                                             names=('mag_auto', 'x_peak_not_used', 'y_peak_not_used',
                                                    'x_image', 'y_image', 'rms_A_image', 'rms_B_image'),
                                             dtype={'x_image': np.float64,
                                                    'y_image': np.float64},
                                             skip_blank_lines=True)
            else:
                raise ValueError('Invalid value "', center_sextr,
                                 '"for center_sextr. Choose peak or auto.')
            # subtract 1 as sextractor starts at (1,1) and not (0,0)
            sex_coords[kk]['x_image'] -= 1
            sex_coords[kk]['y_image'] -= 1
            # remove the ones too close to the border
            sex_coords[kk] = sex_coords[kk][(sex_coords[kk].x_image > cut_size + 4) &
                                            (sex_coords[kk].y_image > cut_size + 4) &
                                            (sex_coords[kk].x_image < (x_size - cut_size - 4)) &
                                            (sex_coords[kk].y_image < (y_size - cut_size - 4))]
            sex_coords[kk] = sex_coords[kk].sort_values(
                by='mag_auto').reset_index(drop=True)

        fits.writeto(os.path.join(dir_out, target + '_objects_cube.fits'),
                     all_objects, header_cube[0], output_verify='ignore')

    if med:
        if verbose:
            print('Processing the median image.')
        fits.writeto('temp.fits', data_cube, output_verify='ignore')
        call(['sex', 'temp.fits', '-c', fn_sextr_med,
              '-checkimage_name', os.path.join(dir_temp,
                                               target + '_objects_med.fits'),
              '-catalog_name', os.path.join(dir_temp, target + '_med.sex')],
             stdout=FNULL, stderr=subprocess.STDOUT)

        call(['rm', 'temp.fits'])
        # now that we have all data in files, read it in
        all_objects = np.full(data_cube.shape, np.nan)
        sex_coords = [None]
        all_objects = fits.getdata(os.path.join(
            dir_temp, target + '_objects_med.fits'))
        sex_coords = pd.read_csv(os.path.join(dir_temp, target + '_med.sex'),
                                 header=None,
                                 delimiter=r'\s+', comment='#',
                                 names=('number', 'mag_auto', 'x_peak_image', 'y_peak_image',
                                        'x_image', 'y_image', 'fwhm', 'class_star'),
                                 dtype={'x_image': np.float64,
                                        'y_image': np.float64},
                                 skip_blank_lines=True)
        # subtract 1 as sextractor starts at (1,1) and not (0,0). convert deg to arcsec
        sex_coords['fwhm'] *= 3600
        sex_coords['x_image'] -= 1
        sex_coords['y_image'] -= 1
        # remove sources with the wrong FWHM
        sex_coords = sex_coords[np.logical_and(sex_coords['fwhm'] >= 0.08,
                                               sex_coords['fwhm'] <= .25)]
        # overwrite the just created fits file with
        sex_coords = sex_coords.sort_values(
            by='mag_auto').reset_index(drop=True)
        # make a peakmap without the ignored sources
        peakmap = np.full(all_objects.shape, 0)
        for x, y, z in zip(sex_coords['x_image'], sex_coords['y_image'],
                           sex_coords['mag_auto']):
            x, y = np.int(np.round(x)), np.int(np.round(y))
            peakmap[y, x] = 10**(-z)
        # overwrite the checkimage
        fits.writeto(os.path.join(dir_out, target + '_objects_med.fits'),
                     peakmap, header_cube[0], output_verify='ignore')
    os.chdir(dir_data)
    if med:
        return sex_coords, peakmap
    else:
        return sex_coords
