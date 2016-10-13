from __future__ import print_function
from __future__ import division
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy import signal
from astropy.io import fits
import pandas as pd
import misc
from parameters import *


def make_artificial_map(source,data_objects):
    '''creates an artificial map with the sources in it'''
    print('Matching the sources...')
    offset = 3*fwhm
    degtopix = 60*60*1000/pxscl_init
    ra_min = min(source['ra_now_rot'])
    dec_min= min(source['dec_now_rot'])
    x_dim = int(round((max(source['ra_now_rot'])-ra_min) *  degtopix + 2*offset))
    y_dim = int(round((max(source['dec_now_rot'])-dec_min)* degtopix + 2*offset))
    source_map = np.zeros((y_dim,x_dim))
    gauss =misc.makeGaussian(3*fwhm,fwhm = fwhm,center = None)
    
    #put the gaussian to the coordinates. Scale with mag in logspace. Real fluxes would overkill.
    x_im =[]
    y_im =[]
    for ii in range(len(source['ra_now_rot'])):
        x = offset + (source['ra_now_rot'][ii] - ra_min) * degtopix
        y = offset + (source['dec_now_rot'][ii]- dec_min)* degtopix
        x_im.append(x)
        y_im.append(y)
        if np.isfinite(source[mag][ii]):
            source_map[int(y-np.floor(0.5*gauss.shape[0])):int(y+np.ceil(0.5*gauss.shape[0])),
                       int(x-np.floor(0.5*gauss.shape[1])):int(x+np.ceil(0.5*gauss.shape[1]))] +=\
                    gauss*(max(source[mag]+1-source[mag][ii]))
        else:
            source_map[int(y-np.floor(0.5*gauss.shape[0])):int(y+np.ceil(0.5*gauss.shape[0])),
                       int(x-np.floor(0.5*gauss.shape[1])):int(x+np.ceil(0.5*gauss.shape[1]))] +=\
                    gauss* 1
    source['x_im_rot']=x_im
    source['y_im_rot']=y_im
    #flip in right direction to align with image and rotate it
    source_map = np.fliplr(source_map)
    source['x_im_rot'] =source_map.shape[1] -source['x_im_rot']
    corr = signal.fftconvolve(source_map,data_objects[::-1,::-1],mode='same')
    y_corr,x_corr = np.unravel_index(np.argmax(corr),corr.shape)

    #plot the results
    print('Plotting the results')
    
    fig, (ax_source,ax_data,ax_corr) = plt.subplots(1,3)
    ax_source.imshow(source_map)    
    ax_source.plot(x_corr,y_corr,'ro')
    ax_source.add_patch(
        mpatches.Rectangle((x_corr-data_objects.shape[1]/2.,y_corr-data_objects.shape[0]/2.),\
                           data_objects.shape[1],data_objects.shape[0],alpha=0.3,fc='r'))
    ax_source.set_title('Catalog objects')
    ax_data.imshow(np.sqrt(data_objects))
    ax_data.set_xlim((0,data_objects.shape[1]))
    ax_data.set_ylim((0,data_objects.shape[0]))
    ax_data.set_title('Sextracted median')
    ax_corr.imshow(corr)
    ax_corr.set_title('Correlated images')
    
    #add all the sources to second plot
    x_data = []
    y_data = []
    for ii in range(len(source['x_im_rot'])):
        x_data.append(source['x_im_rot'][ii] - x_corr + data_objects.shape[1]/2.)                     
        y_data.append(source['y_im_rot'][ii] - y_corr + data_objects.shape[0]/2.)
        ax_data.plot(x_data[ii],y_data[ii],'ro',alpha= 0.3, ms= 4)
        ax_data.annotate(source[star_id][ii],xy=(x_data[ii],y_data[ii]),size=6)
    fig.savefig(dir_out+'correlated_image.svg')
    source['x_data_rot'] = x_data
    source['y_data_rot'] = y_data
    plt.close('all')

    return source


