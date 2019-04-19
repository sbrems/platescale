import os
from astropy.coordinates import SkyCoord
import numpy as np

# dont forget to change sexparams file in sextract.py

#target = 'trapezium'
#agpm = False
# dir_par = "/home/sbrems/Documents/NACO/astrometry/self/" #parent folder to data
# fname_psf = target+'_sat.fits' #set to None if you don't want positioning by crosscorr but by sextractor

# if agpm == False:
#    dir_data = dir_par+target+'/'
# else:
#    dir_data = dir_par+target+'_agpm/'
dir_data = os.getcwd()

dir_out = os.path.join(dir_data, 'results')
dir_temp = os.path.join(dir_out, 'temp')
# directory of calibration files
dir_cat = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cat')

coord_trapez = SkyCoord('05h35m15.84s', '-05d23m22.60s', frame='icrs')
coord_47tuc_sphere = SkyCoord('00h23m58.12s', '-72d05m30.19s', frame='icrs')
coord_47tuc_wolfgang = SkyCoord('00h24m06.30s', '-72d04m59.60s', frame='icrs')
query_around_head_coord = 10.  # in arcsec how far stars in catalog may be away from header position+fov. inf for whole catalog
# parameters of source coordinate file. Comment out for standard nomenclature dependent on target
# give julian date of the coordinates
# mjd_source = 55850  #Close12: 55850
# fn_source = dir_data+'trapezium_close12.csv'  #filename where the coordinates are stored
# fn_sextr_med = dir_data+'median.sex' #filename with the sextractor parameters for the median image
# fn_sextr_sgl = None #dir_data+'single_image.sex' #....for the single images. Might be different
# column names to use
ra = 'RA'
dec = 'DEC'
pm_ra = 'PMRA'
pm_dec = 'PMDEC'
# MAG_USE (MAG_K for close12, MAG_F475W for McLaughlin06, BB fit for Bellini)
mag = 'MAG_USE'
star_id = 'MAIN_ID'


# guess initial pixel scale and true north to correctly identify the sources
# if agpm:
#    pxscl_init = 27.2
#    rot_init   = -1.6317
#  if target == 'trapezium_med':
#      rot_init = -9.56
# else:
pxscl_init = 27.19  # mas/px
true_north_guess = .6
# set True to use header Value, e.g. arcsin(CD2_1/CD1_1), 0.0 else assuming
# it is already rotated to north (according to header value)
use_rot_header = True
mag_lim_identification = 14  # dont use stars fainter for identifying the
# global FoV. But use fainter ones later
# guess where the true north is (EoW) and add it to the header-value/0 if use_rot_header == False. This is a guess for the error of the instrument and only helps identifying the constellatiions.
conv_gauss = False  # choose if the fitting shall be done by convolution and fitting a
# gaussian to the convolution map(=True) or by upsampling (=False) upsampling
# seems better
# keeping the number of star images when making the reference psf. in (0,1]
keepfr = .4
keepfr_med = .9
fwhm = 20  # set fwhm in px for source detection. Bigger is more robust but may lead to more errors
search_rad = 8  # set radius where a source is still accounted as match in px
sigma_outliers = 5.  # give tolerance for sigma clipping
# ignore stars with less good connections as the fraction given (e.g. = 0 keeps all stars,even if only bad connections). Good connetction
ignore_frac = 0.7
max_mvmnt = 5.  # maximal assumed movement (in px) of the stars in the cc images compared to the median image. If the value is too big neighbouring stars will intefere. if too small, tracking of stars doesnt work
# if only one star is found in single image but two in median in search_rad, use the brighter one if the contrast is higher than this given. Else ignore the source.
min_contrast = 0.8
# get rms (fwhm/2) for 4micron at 8.2m telescope
# size to cut around stars for psf fitting. has to be even
cut_size = int(np.ceil(1.44 * 4.e-6 / 8.2 / np.pi * 360 *
                       3600 * 1000 / pxscl_init / 2. / 2.)) * 2
min_dist_to_border = 20 # in px. Min distance to border for reference psf generation
# 'peak' to use brightest pixel,'auto' to use sextractors auto detection
center_sextr = 'peak'
