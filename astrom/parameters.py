import os
from astropy.coordinates import SkyCoord
import numpy as np

###dont forget to change sexparams file in sextract.py

#target = 'trapezium'
#agpm = False
#dir_par = "/home/sbrems/Documents/NACO/astrometry/self/" #parent folder to data
#fname_psf = target+'_sat.fits' #set to None if you don't want positioning by crosscorr but by sextractor

#if agpm == False:
#    dir_data = dir_par+target+'/'
#else:
#    dir_data = dir_par+target+'_agpm/'
dir_data = os.getcwd()+'/'

dir_out = dir_data+'results/'
dir_temp= dir_out +'temp/'

dir_cat = os.path.dirname(os.path.realpath(__file__))+'/cat/' #directory of calibration files

coord_trapez =       SkyCoord('05h35m15.84s','-05d23m22.60s',frame='icrs')
coord_47tuc_sphere=  SkyCoord('00h23m58.12s','-72d05m30.19s',frame='icrs')
coord_47tuc_wolfgang=SkyCoord('00h24m06.30s','-72d04m59.60s',frame='icrs')
query_around_head_coord = 30. #in arcsec how far stars in catalog may be away from header position. inf for whole catalog
#parameters of source coordinate file. Comment out for standard nomenclature dependent on target
#give julian date of the coordinates
#mjd_source = 55850  #Close12: 55850 
#fn_source = dir_data+'trapezium_close12.csv'  #filename where the coordinates are stored
#fn_sextr_med = dir_data+'median.sex' #filename with the sextractor parameters for the median image
#fn_sextr_sgl = None #dir_data+'single_image.sex' #....for the single images. Might be different
#column names to use
ra ='RA'
dec='DEC'
pm_ra='PMRA'
pm_dec='PMDEC'
mag= 'MAG_USE' #MAG_USE (MAG_K for close12, MAG_F475W for McLaughlin06)
star_id='MAIN_ID'


#guess initial pixel scale and true north to correctly identify the sources
#if agpm:
#    pxscl_init = 27.2
#    rot_init   = -1.6317
  #  if target == 'trapezium_med':
  #      rot_init = -9.56
#else:
pxscl_init = 27.19 #mas/px
#rot_init   = -0.6
use_rot_header = True #set True to use header Value, e.g. arcsin(CD2_1/CD1_1)
rot_offset_guess = -1. #add this to the header Value. This is a guess for the error of the instrument
conv_gauss = False #choose if the fitting shall be done by convolution and fitting a 
                   #gaussian to the convolution map(=True) or by upsampling (=False) upsampling
                   #seems better
keepfr = .4 #keeping the number of star images when making the reference psf. in (0,1]
fwhm = 20 #set fwhm in px for source detection. Bigger is more robust
search_rad =20 #set radius where a source is still accounted as match in px
sigma_outliers = 5. #give tolerance for sigma clipping
ignore_frac = 0.7 #ignore stars with less good connections as the fraction given (e.g. = 0 keeps all stars,even if only bad connections). Good connetction 
max_mvmnt = 5. #maximal assumed movement of the stars in the cc images compared to the median image. If the value is too big neighbouring stars will intefere. if too small, tracking of stars doesnt work
min_contrast = 1.5 #if only one star is found in single image but two in median in search_rad, use the brighter one if the contrast is higher than this given. Else ignore the source.
#get rms (fwhm/2) for 4micron at 8.2m telescope
cut_size = int(np.ceil(1.44*4e-6 /8.2 /np.pi*360*3600*1000/pxscl_init/2./2.)*2) #size to cut around stars for psf fitting. has to be even
center_sextr = 'peak' #'peak' to use brightest pixel,'auto' to use sextractors auto detection
