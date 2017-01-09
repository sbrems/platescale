from __future__ import print_function,division
import numpy as np
import matplotlib.pyplot as plt
import math
from parameters import *
import dd_dec_converter as dc
from shutil import copy2
from astropy.coordinates import SkyCoord
import pandas as pd
import misc

def makeGaussian(size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def dist_ang(a,b,spherical=False):
    '''Quick distance and angle calculator for 2x(y,x)-coordinates.
    Returns dist and angle. Angle is calculated around the first
    point given with respect to the positive x-axis, then 
    counter-clockwise in radians. If spherical
    x[0] is RA and x[1] is DEC. Assuming them given in deg, return angle
    in deg, dist in mas.'''
    if spherical:
        aa = SkyCoord(a[0],a[1],frame='icrs',unit='deg')
        bb = SkyCoord(b[0],b[1],frame='icrs',unit='deg')
        dist = aa.separation(bb).mas
        angle = aa.position_angle(bb).deg #in east of north
    else:
        delta_ra = b[0]-a[0]
        delta_dec= b[1]-a[1]
        dist  = math.sqrt(delta_ra**2 + delta_dec**2)
        angle = math.atan2(delta_ra,delta_dec)%(2*np.pi)
    return dist,angle


def rotate_coords(x,y,angle,rad = False,center=[0.,0.]):
    '''Rotate the coordinates around (0,0) by the angle given. 
    Angle in degrees, x,y are np.arrays. rad = True uses
    radians'''

    if not rad:#convert deg to rad
        angle *= 2*math.pi/360.
    x = np.array(x)
    y = np.array(y)
    x_dum = x - center[0]
    y_dum = y - center[1]
    x = x_dum * math.cos(angle) - y_dum * math.sin(angle)
    y = x_dum * math.sin(angle) + y_dum * math.cos(angle)
    x += center[0]
    y += center[1]
    return x,y
    
def get_catalogue_data(fn_source,mjd_source,coord_head,mjd_obs=None,verbose=True,plot=True):
    #Reads the positions in RA and DEC from the catalogue given. if mjd is given it calculates
    #the current position based on the pm found in that catalogue
    
    #read in the coords and convert them
    if verbose: print('Converting coordinates. Assuming hmsdms and mas/yr given in',fn_source,pm_ra,pm_dec)
    source = pd.read_csv(fn_source,delimiter =',')
    for numeric in [mag,pm_ra,pm_dec]:#convert to floats
        source[numeric] = pd.to_numeric(source[numeric],errors='coerce')
    rem_coord = []

    if len(str(source[ra][00]).strip().split(' ')) == 3: #assuming hms/dms
        if verbose: print('Converting coordinates. Assuming hms/dms and mas/yr given in',fn_source,pm_ra,pm_dec)
        coord_fmt = 'hmsdms'
    else: #assuming deg
        if verbose: print('Converting coordinates. Assuming deg and mas/yr given in',fn_source,pm_ra,pm_dec) 
        coord_fmt = 'deg'
    for ii in xrange(len(source[ra])):
        if coord_fmt == 'hmsdms':
            coord_dum = SkyCoord(dc.hmsToDeg(source[ra][ii].strip(), sep=' '), \
                                 dc.dmsToDeg(source[dec][ii].strip(),sep=' '), \
                                 frame='icrs',unit='deg')
        elif coord_fmt == 'deg':
            coord_dum = SkyCoord(source[ra][ii],source[dec][ii],frame='icrs',unit='deg')
        else: raise ValueError('Unknown coordinate format. Use hmsdms with space separation or deg')

        if coord_head.separation(coord_dum).arcsec <= query_around_head_coord:
            source[ra][ii] = coord_dum.ra.deg
            source[dec][ii]= coord_dum.dec.deg
        else: 
            rem_coord.append(ii)
                    
    source = source.drop(rem_coord)
    source.reset_index(inplace=True)
    

    #check for validity of mjd_obs and set to source date (zero pm) if none is given
    if mjd_obs != None:
        if mjd_obs <= 57300 or mjd_obs >= 65000:
            raise ValueError('mjd of ',mjd_obs,' smaller than 57300! This should not be the case for ISPY')
    else:
        mjd_obs = mjd_source
    
    #calculate the current positions and make some plots
    ra_now =[]
    dec_now=[]
    for ii in range(len(source[ra])):
        if np.isfinite(source[pm_ra][ii]):
            ra_now.append(source[ra][ii]+ (mjd_obs - mjd_source)/365.2425 * source[pm_ra][ii] / (1000.*60.*60.))
            dec_now.append(source[dec][ii]+ (mjd_obs - mjd_source)/365.2425 * source[pm_dec][ii] / (1000.*60.*60.))
            if (abs(ra_now[ii]-source[ra][ii] >= 1)) | (abs(dec_now[ii]-source[dec][ii] >= 1)):
                raise ValueError('PM too high for target', source[star_id])
        else:
            ra_now.append(source[ra][ii])
            dec_now.append(source[dec][ii])
            
    ra_now = np.array(ra_now)
    dec_now= np.array(dec_now)
    source['ra_now'] = ra_now
    source['dec_now']= dec_now
#    source['ra_now_rot'],source['dec_now_rot'] = misc.rotate_coords(ra_now,dec_now,-rot_init,\
#                                                    center=[np.mean(ra_now),np.mean(dec_now)])
    #plot the result
    if plot:
        plt.xlabel('RA [deg]')
        plt.ylabel('DEC [deg]')
        plt.scatter(source[ra],source[dec],label = 'Orig. Pos.',color='blue')
        plt.scatter(source['ra_now'],source['dec_now'], color='red',label='With PM',alpha = 0.3)
       # plt.scatter(source['ra_now_rot'],source['dec_now_rot'], color='k',label='With PM and rotated',alpha = 0.3)
        plt.legend()
        plt.savefig(dir_out+'catalogue_stars.svg')
        plt.close('all')
    
    return source
    
def reject_outliers(data, m=2):
    '''returns the mx sigma clipped 1-d array'''
    return data[abs(data - np.median(data)) < m * np.std(data)]

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

def get_headerparams(header,verbose=True):
    targetname = header['OBJECT'].lower()
    target = None
    coord_head = SkyCoord(header['CRVAL1'],header['CRVAL2'],'icrs',unit='deg')
    while target == None:
        if ('trapezium' in targetname) or (coord_head.separation(coord_trapez).arcsec < 30) :
            target = 'trapezium'
        elif ('47tuc' in targetname) \
             or (coord_head.separation(coord_47tuc_sphere).arcsec < 30)\
             or (coord_head.separation(coord_47tuc_wolfgang).arcsec < 30):
            target = '47tuc'
        else:
            print('ATTENTION!!! Target <',targetname,'> unknown. Please type trapezium or 47tuc .')
            targetname = str(raw_input())
    if 'apo_165' in header['Apodizer'].lower():
        agpm=True
    else:
        agpm=False
    mjd_obs=None
    mjd_obs = header['MJD-OBS']
    
    if use_rot_header:
        rot_header = np.rad2deg(np.arcsin(header['CD2_1']/header['CD1_1']))
    else:
        rot_header = 0.0

    if verbose: print('Found the following configuration:\n AGPM: %s \n target: %s \n mjd_obs: %s \n rot_init: %s' \
                      %(agpm,target,mjd_obs,rot_header))
    

    return target,mjd_obs,agpm,coord_head,rot_header

def get_catfiles(target,agpm,dir_cat,dir_temp,verbose=True):
    try: fn_sextr_med
    except NameError: fn_sextr_med = dir_cat+target+'_med.sex'
    try: fn_sextr_sgl
    except NameError: fn_sextr_sgl = dir_cat+target+'_single.sex'
    try: fn_source
    except NameError: fn_source = dir_cat+target+'.csv'
    try: mjd_source
    except NameError:
        if target == 'trapezium':
            mjd_source = 55850.
        if target == '47tuc':
            mjd_source = 52369.5
    if verbose: print('Using and copying the following parameterfiles: \n sextractor median:%s \n sextractor single: %s \n source catalog: %s \n mjd of source cat: %s' %(fn_sextr_med,fn_sextr_sgl,fn_source,mjd_source))
    for fn in [fn_sextr_med,fn_sextr_sgl,fn_source]:
        copy2(fn,dir_temp+os.path.basename(fn))
        fn = os.path.basename(fn)
    copy2(dir_cat+'sex.param',dir_temp+'sex.param')
    return mjd_source,fn_sextr_med,fn_sextr_sgl,fn_source
