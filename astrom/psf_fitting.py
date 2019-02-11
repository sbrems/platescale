
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import fftpack
# from scipy import interpolate
from scipy import signal
# from scipy.optimize import minimize
import scipy.optimize as opt
from scipy.optimize import OptimizeWarning
from subprocess import call
import skimage.feature
import warnings
import os
import subprocess
from .misc import twoD_Gaussian
from .parameters import cut_size, min_dist_to_border


def make_psf(data_cube, sex_coords, dir_temp, ign_ims=[], n_rms=2,
             conv_gauss=True, keepfr=0.7, fn_out='psf_cube', verbose=True):
    '''Makes a psf out of the stars in the images found by fitting and median combining them.n_rms
    gives the image size in rms of the ellipse parameter of the brightest star found by sextractor
    to which the images are cropped.'''
    x_dim = data_cube.shape[2]
    y_dim = data_cube.shape[1]
    process_images = [x for x in range(data_cube.shape[0]) if x not in ign_ims]
    n_images = len(sex_coords)
    n_stars = 0
    # get rms (fwhm/2) for 4micron at 8.2m telescope
    # cut_size = int(np.ceil(1.44*4e-6 /8.2 /math.pi*360*3600*1000/pxscl_init/2./2.)*2)
    for i_image in range(n_images):
        n_stars += len(sex_coords[i_image])
        sex_coords[i_image].reset_index(inplace=True, drop=True)

    i_max = sex_coords[0]['mag_auto'].idxmin()
# uncomment for sextractor values
#    cut_size = int(round(n_rms* math.sqrt(math.pow(sex_coords[0]['rms_A_image'][i_max],2) +\
#                                      math.pow(sex_coords[0]['rms_B_image'][i_max],2))))
    # make cut_size even
    # cut_size = int(np.ceil(cut_size / 2.) * 2)
    # Now cut the image around this region and update x/y_dim in
    # case it was at the border
    x_star = int(sex_coords[0]['x_image'][i_max])
    y_star = int(sex_coords[0]['y_image'][i_max])
    brightest_star = data_cube[process_images[0],
                               y_star - cut_size: y_star + cut_size + 1,
                               x_star - cut_size: x_star + cut_size + 1]
    x_dim = brightest_star.shape[1]
    y_dim = brightest_star.shape[0]
    data_stars = np.full((n_stars, y_dim, x_dim), np.nan, dtype=np.float64)
    shifts = np.full((n_stars, 2), np.nan, dtype=np.float64)
    data_stars[0, :, :] = brightest_star
    shifts[0, :] = (0, 0)

    del brightest_star
    # now cut out the other stars
    i_star_tot = 1  # already have the ref-star
    max_shift = 4
    n_removed_stars = 0  # count stars could not be fit properly
    if verbose:
        print('Getting PSFs from ', n_images, 'images')

    for i_image, i_use_image in zip(list(range(n_images)), process_images):
        # if verbose: print('Getting PSF for image ',i_image+1,' of ',n_images)
        for i_star in range(len(sex_coords[i_image])):
            use_star = True  # set to true if fitting worked. otherwise dont use the star
            if i_image != 0 or i_star != i_max:
                x_cen = int(round(sex_coords[i_image]['x_image'][i_star]))
                y_cen = int(round(sex_coords[i_image]['y_image'][i_star]))
                cut_star = data_cube[i_use_image,
                                     y_cen - cut_size:y_cen + cut_size + 1,
                                     x_cen - cut_size:x_cen + cut_size + 1].copy()
                # if the cut_size is too close to border, dont use for ref-psf generation
                if (y_cen < min_dist_to_border) or\
                   (x_cen < min_dist_to_border) or\
                   (data_cube.shape[1] - min_dist_to_border < y_cen) or\
                   (data_cube.shape[2] - min_dist_to_border < x_cen):
                    use_star = False

                if conv_gauss and use_star:
                    shifts[i_star_tot, :] = find_shift(data_stars[0, :, :], cut_star, n_rms=n_rms,
                                                       plotpath='results/temp/aligned_to_psf_{:03}.pdf'.format(i_star_tot))[0]
                elif use_star:
                    found_shift = find_shift_upsampling(data_stars[0, :, :], cut_star,
                                                        max_shift=max_shift, try_smaller_fov=False)[0]
                    if np.sqrt(found_shift[0]**2 + found_shift[1]**2) > max_shift:
                        # ignore the star
                        use_star = False
                    else:
                        shifts[i_star_tot, :] = found_shift
                        use_star = True
                # shift the image star region +-3px and cut the star dont shift whole image cause
                # fft_shift can't treat nans
                if use_star:
                    shifted_data = fft_shift(data_cube[i_use_image,
                                                       y_cen - cut_size - 3:y_cen + cut_size + 3 + 1,
                                                       x_cen - cut_size - 3:x_cen + cut_size + 3 + 1],
                                             shifts[i_star_tot][1], shifts[i_star_tot][0])[3:-3, 3:-3]
                    # check if only infs are left
                    if np.count_nonzero(np.isnan(shifted_data)) == shifted_data.size:
                        raise ValueError(
                            'Only nan array detected when shifting image at make psf!')
                    data_stars[i_star_tot, :, :] = (
                        data_stars[0].max() / shifted_data.max()) * shifted_data

                    i_star_tot += 1
                else:
                    n_removed_stars += 1
                    shifts = np.delete(shifts, i_star_tot, 0)
                    data_stars = np.delete(data_stars, i_star_tot, 0)
                    n_stars -= 1
    # scale up all to the same brightness
    pix_max = np.float64(np.nanmax(data_stars))
    for i_star in range(n_stars):
        data_stars[i_star] *= pix_max / np.max(data_stars[i_star])
    first_med = np.median(data_stars, axis=0)
    ###########################################################################################
    ############Only keep keepfr of the images with the least residuals to the median######
    ###########################################################################################
    residuals = np.full((n_stars), np.nan)
    for i_star in range(n_stars):
        residuals[i_star] = np.sum((data_stars[i_star] - first_med)**2)

    sorted_resid = np.argsort(residuals)
    selected_resid = sorted_resid[0:int(round(keepfr * n_stars))]
    selected_stars = data_stars[selected_resid, :, :]

    second_med = np.median(selected_stars, axis=0)

    fits.writeto(os.path.join(dir_temp, fn_out + '_w_median.fits'),
                 np.concatenate((data_stars, [first_med]), axis=0))
    fits.writeto(os.path.join(dir_temp, fn_out + '_selected_w_median.fits'),
                 np.concatenate((selected_stars, [second_med]), axis=0))
    return second_med


def refine_fit(data_cube, data_psf, sex_coords, ign_ims=[], conv_gauss=True, verbose=True, plot=True):
    '''Resets the x_image and y_image from sextractor via psf. Use convolution with gaussian fitting
    if conv_gauss==True. Else use upsampling and find the max of crosscorrelation.
    Sex coords has to be a list of pandas data frames.'''
    # only use the images which where not discarded before
    data_cube = data_cube[[x for x in range(
        data_cube.shape[0]) if x not in ign_ims], :, :]
    x_dim = data_cube.shape[2]
    y_dim = data_cube.shape[1]
    n_images = data_cube.shape[0]
    n_stars = 0
    for i_image in range(n_images):
        n_stars += len(sex_coords[i_image])
    # make sure data_psf has an uneven shape
    if data_psf.shape[0] % 2 == 0:
        data_psf = data_psf[1::, ]
    if data_psf.shape[1] % 2 == 0:
        data_psf = data_psf[::, 1::]
    size_y = data_psf.shape[0]
    size_x = data_psf.shape[1]

    i_stars_tot = 0
    for i_image in range(n_images):
        shifts = np.full(
            (len(sex_coords[i_image]), 2), np.nan, dtype=np.float64)
        errors = np.full((len(sex_coords[i_image])), np.nan, dtype=np.float64)
        # copy the old coords
        sex_coords[i_image]['x_image_sex'] = sex_coords[i_image]['x_image']
        sex_coords[i_image]['y_image_sex'] = sex_coords[i_image]['y_image']
        for i_star in range(len(sex_coords[i_image])):
            # first cut the image around the found position
            x_old = int(round(sex_coords[i_image]['x_image'][i_star]))
            y_old = int(round(sex_coords[i_image]['y_image'][i_star]))

            # remove stars which are too close to the border
            if ((x_dim - x_old) <= size_x) or\
               (x_old <= size_x) or\
               ((y_dim - y_old) <= size_y) or\
               (y_old <= size_y):
                print('Not fitting star {} in image {} as it is too close to the \
border. Be cautious! Its weight is put to 0. e.g. not used in robust values.\
'.format(i_star, i_image))
                sex_coords[i_image]['x_image'][i_star] = sex_coords[i_image]['x_image'][i_star]
                sex_coords[i_image]['y_image'][i_star] = sex_coords[i_image]['y_image'][i_star]
                errors[i_star] = np.inf
                shifts[i_star,:] = [np.nan, np.nan]

            else:  # star far enough from edges
                cut_star = data_cube[i_image,
                                     y_old - int(np.floor(size_y / 2.)):y_old +
                                     int(np.ceil(size_y / 2.)),
                                     x_old - int(np.floor(size_x / 2.)):x_old +
                                     int(np.ceil(size_x / 2.))]
                if conv_gauss:
                    if plot:
                        plotpath = 'results/temp/psf_aligned_{:03}.pdf'.format(i_star)
                    else:
                        plotpath = None
                    shifts[i_star, :], errors[i_star] = find_shift(data_psf, cut_star,
                                                                   plotpath=plotpath)
                else:
                    shifts[i_star, :], errors[i_star] = find_shift_upsampling(
                        data_psf, cut_star)
    #                shifts[i_star,1] *= -1 #correct for different y convention
                # update the pd frame
                sex_coords[i_image]['x_image'][i_star] = x_old - \
                    shifts[i_star, 0]
                sex_coords[i_image]['y_image'][i_star] = y_old - \
                    shifts[i_star, 1]
                sex_coords[i_image]['shift_error'] = errors
            # if verbose: print('Shifts found in image nr',i_image,':\n',shifts)

    return sex_coords


def find_shift(ref_frame, shift_frame, plotpath=None, n_rms=2, convolve=True):
    '''Find the shift in an image by gaussian fitting of the convolution map. Res in x,y
    of the direction in which the shift image is offset (e.g. move shift_image by -res to align).
    perr is the error.
    If you want the plots give a plotpath.'''
    # first make images non-negative
    shift_frame -= np.min(shift_frame)
    ref_frame -= np.min(ref_frame)
    conv = signal.correlate2d(ref_frame, shift_frame, mode='same')
    x = conv.shape[1]
    y = conv.shape[0]
    xx, yy = np.meshgrid(list(range(x)), list(range(y)))
    # use max_values as initial guesses and
    ymax, xmax = np.unravel_index(np.nanargmax(conv), conv.shape)
    initial_guess = (np.nanmax(conv), ymax, xmax,
                     x / n_rms, y / n_rms, 0, 0)
    with warnings.catch_warnings():
        warnings.simplefilter('error', OptimizeWarning)
        try:
            popt, pcov = opt.curve_fit(
                twoD_Gaussian, (xx, yy), conv.flatten(), p0=initial_guess)
            res = (popt[1] - x / 2., popt[2] - y / 2.)
            perr = np.sqrt(np.diag(pcov))
            if popt[0] <= 0:
                raise OptimizeWarning(
                    'Fit was negative what it shouldnt be:', popt[0])
        except:  # if it fails, try a smaller FoV
            print('Fit failed. Trying a smaller FoV...')
            quarter_size = int(round(conv.shape[0] / 3))
            conv = conv[quarter_size:-quarter_size, quarter_size:-quarter_size]
            x = conv.shape[1]
            y = conv.shape[0]
            xx, yy = np.meshgrid(list(range(x)), list(range(y)))
            # use max_values as initial guesses and
            ymax, xmax = np.unravel_index(np.nanargmax(conv), conv.shape)
            initial_guess = (np.nanmax(conv), ymax, xmax,
                             x / n_rms, y / n_rms, 0, 0)
            try:
                popt, pcov = opt.curve_fit(
                    twoD_Gaussian, (xx, yy), conv.flatten(), p0=initial_guess)
                # bounds = ([0,ymax-2,xmax-2,0,0,0,0],
                #          [np.inf,ymax+2,xmax+2,np.inf,np.inf,np.inf,np.inf]))
                res = (popt[1] - x / 2., popt[2] - y / 2.)
                perr = np.sqrt(np.diag(pcov))
                if (not 0 < popt[1] < x) or (not 0 < popt[2] < y):
                    raise ValueError(
                        'Also smaller field didnt work.')

            except:
                print('Trying upsampling')
                y, x = ref_frame.shape
                popt[1:3], perr = skimage.feature.register_translation(ref_frame, shift_frame,
                                                                       upsample_factor=20, space='real')[0:2]
                res = (popt[1] - x / 2., popt[2] - y / 2.)
                plotpath = None

    # make some plots
    if plotpath is not None:
        data_fitted = twoD_Gaussian((xx, yy), *popt)

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # ax1.hold(True) #depreciated
        ax1.imshow(conv, cmap=plt.cm.jet, origin='bottom',
                   extent=(xx.min(), xx.max(), yy.min(), yy.max()), interpolation='none')
        ax1.contour(xx, yy, data_fitted.reshape(y, x), 20, colors='w')
        ax1.plot(res[0] + x / 2., res[1] + x / 2., 'ro')
        ax1.set_title('Correlation map')

        ax2.imshow(ref_frame)
        ax2.set_title('ref_frame (=psf?)')

        ax3.imshow(shift_frame)
        ax3.set_title('shift_frame (=cut star?)')

        plt.savefig(plotpath)
        plt.close('all')

    return res, perr


# shift_detection with register_translation
def find_shift_upsampling(ref_frame, shift_frame, max_shift=4,
                          try_smaller_fov=True):
    # make images non-negative. There seems to be a bug else
    ref_frame -= np.min(ref_frame)
    shift_frame -= np.min(shift_frame)
    res, perr = skimage.feature.register_translation(ref_frame, shift_frame,
                                                     upsample_factor=30,
                                                     space='real')[0:2]
    res = res[::-1]
    # if shift too big, downscale shift_frame
    if (math.sqrt(math.pow(res[0], 2) + math.pow(res[1], 2)) >= max_shift) & (try_smaller_fov):
        print('Shift of %s considered too big. Trying smaller FoV now.' % (res))
        quarter_size = int(round(shift_frame.shape[0] / 4))
        shift_frame_q = shift_frame[quarter_size: -quarter_size,
                                    quarter_size: -quarter_size]
        ref_frame_q = ref_frame[quarter_size: -quarter_size,
                                quarter_size: -quarter_size]
        res, perr = skimage.feature.register_translation(ref_frame_q, shift_frame_q,
                                                         upsample_factor=20, space='real')[0:2]
        res = res[::-1]
        if math.sqrt(math.pow(res[0], 2) + math.pow(res[1], 2)) >= max_shift:
            raise ValueError(
                'Shift of %s still too big...maybe try with gaussian?' % (res))
        else:
            print('That worked. Shift now:', res)

    return res, perr


# shift detection with cubic interpolation (didnt work well):

# def find_shift(ref_frame,shift_frame):
#     '''Find the shift in an image by cubic interpolation of the convolution map'''
#     conv = signal.correlate2d(ref_frame,shift_frame)
#     xx = range(conv.shape[1])
#     yy = range(conv.shape[0])
#     f = interpolate.interp2d(xx,yy,-conv,kind='cubic',fill_value=np.nan)
#     f2 = lambda x: f(*x)
#     bnds=((0,ref_frame.shape[1]),(0,ref_frame.shape[0]))
#     arr = np.full((300,300),np.nan)
#     for i,ii in enumerate(np.linspace(0,45,300)):
#         for j,jj in enumerate(np.linspace(0,45,300)):
#             arr[i,j] = f2((ii,jj))
#     res =minimize(f2,(ref_frame.shape[1]/2.,ref_frame.shape[0]/2.),
#                    method='SLSQP',bounds=bnds)
#     ipdb.set_trace()
#     if not res['success']:#that means something failed
#         raise ValueError('Aborting. Somehow the fitting for the correlation of psf\
# didnt work. Error message:',res['message'])

#     return res['x']


def fft_shift(in_frame, dx, dy):
    f_frame = fftpack.fft2(in_frame)
    N = fftpack.fftfreq(in_frame.shape[0])
    v = np.exp(-2j * np.pi * dx * N)
    u = np.exp(-2j * np.pi * dy * N)
    f_frame = f_frame * u
    f_frame = (f_frame.T * v).T

    return np.real(fftpack.ifft2(f_frame))


def find_via_ccmap(data_cube, psf, target, dir_data, dir_temp):
    print('Finding sources via ccmap. This should exclude fake multiples.')
    if len(data_cube.shape) == 3:
        n_images = data_cube.shape[0]
    elif len(data_cube.shape) == 2:
        n_images = 1
        data_cube = [data_cube, ]
    else:
        raise ValueError(
            'in find_via_ccmap:data_cube has unknown shape: ', data_cube.shape)

    for kk in range(n_images):
        conv = signal.correlate2d(data_cube[kk, :, :], psf, mode='same')
        fits.writeto('temp.fits', conv, output_verify='ignore')
        FNULL = open(os.devnull, 'w')
        call(['sex', 'temp.fits', '-c', dir_data + 'ccmap.sex',
              '-checkimage_name', dir_temp + target +
              '_ccmap_objects_' + '{:03}'.format(kk) + '.fits',
              '-catalog_name', dir_temp + target + '_ccmap_' + '{:03}'.format(kk) + '.sex'], stdout=FNULL,
             stderr=subprocess.STDOUT)
        call(['rm', 'temp.fits'])
    sex_coords_ccmap = [None] * n_images

    # now read in the data and subtract 1 as sextractor stars counting at (1,1)
    for kk in range(n_images):
        sex_coords_ccmap[kk] = pd.read_csv(target + '_' + '{:03}'.format(kk) + '.sex', header=None,
                                           delimiter=r'\s+', comment='#',
                                           names=(
            'mag_auto', 'x_image', 'y_image', 'rms_A_image', 'rms_B_image'),
            dtype={'x_image': np.float64,
                   'y_image': np.float64},
            skip_blank_lines=True)

        sex_coords_ccmap[kk]['x_image'] -= 1
        sex_coords_ccmap[kk]['y_image'] -= 1

    return sex_coords_ccmap
