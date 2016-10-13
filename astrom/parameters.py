import os

###dont forget to change sexparams file in sextract.py

target = 'trapezium_med'
agpm = True
dir_par = os.getcwd() #parent folder to data
fname_psf = target+'_sat.fits' #set to None if you don't want positioning by crosscorr but by sextractor

if agpm == False:
    dir_data = dir_par+'/'+target+'/'
else:
    dir_data = dir_par+'/'+target+'_agpm/'

dir_out = dir_data+'results/'
dir_temp= dir_out +'temp/'


#parameters of source coordinate file
#give julian date of the coordinates
mjd_source = 55850 #Close12: 55850 
fn_source = dir_data+'trapezium_close12.csv'  #filename where the coordinates are stored
#column names to use
ra ='RA'
dec='DEC'
pm_ra='PMRA'
pm_dec='PMDEC'
mag='MAG_K'
star_id='MAIN_ID'


#guess initial pixel scale and true north to correctly identify the sources
if agpm:
    pxscl_init = 27.2
    rot_init   = -1.6317
    if target == 'trapezium_med':
        rot_init = -9.56
else:
    pxscl_init = 27.2
    rot_init   = -0.56

conv_gauss = False #choose if the fitting shall be done by convolution and fitting a 
                   #gaussian to the convolution map(=True) or by upsampling (=False) upsampling
                   #seems better
keepfr = .6 #keeping the number of star images when making the reference psf. in (0,1]
fwhm = 20 #set fwhm in px for source detection. Bigger is more robust
search_rad = 20 #set radius where a source is still accounted as match in px
sigma_outliers = 5. #give tolerance for sigma clipping
ignore_frac = 0.7 #ignore stars with less good connections as the fraction given (e.g. = 0 keeps all stars,even if only bad connections). Good connetction 
