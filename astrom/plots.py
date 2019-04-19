
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy import signal
# from astropy.io import fits
# import pandas as pd
from . import misc
from .parameters import pxscl_init, fwhm, mag, star_id,\
    true_north_guess, mag_lim_identification


def make_artificial_map(source, data_objects, rot_header, verbose=True,
                        plot=True, use_mags=True, dir_out=None):
    '''creates an artificial map with the sources in it'''
    print('Matching the sources...')
    rot_init = rot_header - true_north_guess
    offset = 10 * fwhm
    degtopix = 60 * 60 * 1000 / pxscl_init
    ra_min = min(source['ra_now'])
    dec_min = min(source['dec_now'])
    x_dim = max(int(round((max(source['ra_now']) - ra_min) *
                          np.cos(dec_min * np.pi / 180.) *
                          degtopix + 2 * offset)),
                int(round((max(source['dec_now']) - dec_min) * degtopix +
                          2 * offset)))
    y_dim = x_dim
    source_map = np.zeros((y_dim, x_dim))
    gauss = misc.makeGaussian(3 * fwhm, fwhm=fwhm, center=None)
    # gauss = np.array([[1,]])

    # put the gaussian to the coordinates. Scale with mag in logspace. Real fluxes would only match the brightest source.
    x_im = []
    y_im = []
    x_im_rot = []
    y_im_rot = []
    for ii in range(len(source['ra_now'])):
        x = offset + (source['ra_now'][ii] - ra_min) * \
            np.cos(source['dec_now'][ii] * np.pi / 180.) * degtopix
        y = offset + (source['dec_now'][ii] - dec_min) * degtopix
        x_rot, y_rot = misc.rotate_coords(
            x, y, -rot_init, center=[int(x_dim / 2.), int(y_dim / 2.)])

        x_im.append(x)
        y_im.append(y)
        x_im_rot.append(x_rot)
        y_im_rot.append(y_rot)
        if use_mags and (source[mag][ii] <= mag_lim_identification):
            source_map[int(y_rot - np.floor(0.5 * gauss.shape[0])):int(y_rot + np.ceil(0.5 * gauss.shape[0])),
                       int(x_rot - np.floor(0.5 * gauss.shape[1])):int(x_rot + np.ceil(0.5 * gauss.shape[1]))] +=\
                np.array(
                    gauss * (10.**((max(source[mag]) - source[mag][ii]) / 2.5)))
            # gauss*(max(source[mag]+1-source[mag][ii]))
        else:
            source_map[int(y_rot - np.floor(0.5 * gauss.shape[0])):
                       int(y_rot + np.ceil(0.5 * gauss.shape[0])),
                       int(x_rot - np.floor(0.5 * gauss.shape[1])):
                       int(x_rot + np.ceil(0.5 * gauss.shape[1]))] +=\
                gauss * 1
    source_map = source_map**0.1
    source['x_im_rot'] = x_im_rot
    source['y_im_rot'] = y_im_rot
    # flip in right direction to align with image and rotate it
    source_map = np.fliplr(source_map)
    source['x_im_rot'] = source_map.shape[1] - source['x_im_rot']
    data_objects[data_objects < 0] = 0.

    plt.imshow(source_map, origin='lower')
    plt.savefig(os.path.join(dir_out, 'artificial_source_map.pdf'))
    plt.close()
    
    corr = signal.fftconvolve(
        source_map, (data_objects[::-1, ::-1]**0.1), mode='same')
    cut_edges = 100 # dont take points too close to corner, in px
    y_corr, x_corr = np.array(np.unravel_index(np.argmax(corr[cut_edges:-cut_edges,
                                                     cut_edges:-cut_edges]),
                                      corr[cut_edges:-cut_edges,
                                           cut_edges:-cut_edges].shape)) + cut_edges
    # plot the results
    if verbose:
        print('Plotting the results')
    if plot:
        if verbose:
            print('Plotting the results')
        fig, (ax_source, ax_data, ax_corr) = plt.subplots(1, 3)
        ax_source.imshow(
            (source_map[::-1, :] + np.min(source_map))**0.5)
        ax_source.plot(x_corr, source_map.shape[0] - y_corr, 'ro', alpha=0.3)
        ax_source.add_patch(
            mpatches.Rectangle((x_corr - data_objects.shape[1] / 2.,
                                source_map.shape[0] - (y_corr) -
                                data_objects.shape[0] / 2.),
                               data_objects.shape[1], data_objects.shape[0],
                               alpha=0.3, fc='r'))
        ax_source.set_title('Catalog objects')
        ax_data.imshow(np.log10(data_objects))
        ax_data.set_xlim((0, data_objects.shape[1]))
        ax_data.set_ylim((0, data_objects.shape[0]))
        ax_data.set_title('Sextracted median')
        ax_corr.imshow(corr[::-1, :])
        ax_corr.set_title('Correlated images')

    # add all the sources to second plot
    x_data = []
    y_data = []
    for ii in range(len(source['x_im_rot'])):
        x_data.append(source['x_im_rot'][ii] - x_corr +
                      data_objects.shape[1] / 2.)
        y_data.append(source['y_im_rot'][ii] - y_corr +
                      data_objects.shape[0] / 2.)
        if plot:
            ax_data.plot(x_data[ii], y_data[ii], 'ro', alpha=0.3, ms=4)
            ax_data.annotate(source[star_id][ii], xy=(
                x_data[ii], y_data[ii]), size=6)
    if plot:
        fig.savefig(os.path.join(dir_out, 'correlated_image.pdf'))
        plt.close('all')
    source['x_data_rot'] = x_data
    source['y_data_rot'] = y_data

    return source


def median_stars(data_med, source_median, ignored_med,
                 source_cat, dir_out=None):
    '''making some nice plot on the median image with the sources found
    and the sources from the catalogue. Also the ignored ones are plotted.
    Input: median_image (array),connected_source_list, ingored_med,
    (of sextraceted median_sources),source_cat. Nothing returned.'''
    fig, im = plt.subplots()
    im.imshow(np.sqrt(data_med), interpolation='none')
    im.set_xlim((0, data_med.shape[1]))
    im.set_ylim((0, data_med.shape[0]))
    # plot the found sources with their names and positions. connect them to the catalogue ones.
    im.scatter(source_cat['x_data_rot'], source_cat['y_data_rot'],
               edgecolors='red', facecolors='none', marker='o', s=7, label='catalogue')
    im.scatter(source_median['x_image'], source_median['y_image'],
               edgecolors='green', facecolors='none', marker='o', s=7, label='sextracted')
    for idx in source_median.index:
        im.plot((source_median['x_data_rot'][idx], source_median['x_image'][idx]),
                (source_median['y_data_rot'][idx],
                 source_median['y_image'][idx]),
                'black', linewidth=1)
    im.scatter(ignored_med['x_image'], ignored_med['y_image'],
               color='orange', marker='x', s=15, label='ignored sextr')
    im.legend(loc=3)
    plt.savefig(os.path.join(dir_out, 'sextracted_from_median.pdf'))
    plt.close('all')


def distance_distributions(astrometry, astrometry_grouped, dir_out=None):
    '''Make some plots showing the different distributions for the platescale/rotation'''
    print('Plotting the distributions for the different conections')
    group = astrometry.groupby([star_id + '1', star_id + '2'])
    combinations = list(group.groups.keys())
    # define the colors
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(combinations))))
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    for ii, comb in enumerate(combinations):
        col = next(color)
        this_group = group.get_group(comb)
        n_tg = len(this_group)
        for i_tg in range(n_tg):
            ax1.scatter(ii + np.float(i_tg) / n_tg, this_group['platescale'].iloc[i_tg], c=col,
                        label=str(comb))
            ax2.scatter(ii + np.float(i_tg) / n_tg, this_group['ang_diff'].iloc[i_tg], c=col,
                        label=str(comb))
    ax1.set_title('Platescales for the different connections')
    ax2.set_title('True north for the different connections')
    ax1.set_xlabel('Combination nr')
    ax2.set_xlabel('Combination nr')
    ax1.set_ylabel('Platescale [mas/px]')
    ax2.set_ylabel('True north [deg]')
    ax1.grid(True)
    ax2.grid(True)
    ax1.set_xticklabels(combinations, rotation=20)
    plt.savefig(os.path.join(dir_out, 'scatter_plot.pdf'))
    plt.close('all')
    # also do some histogram plots
    for ii, comb in enumerate(combinations):
        fig = plt.figure()
        hist = fig.add_subplot(111)
        this_group = group.get_group(comb)
        n_tg = len(this_group)
        hist.hist(this_group['platescale'], 20, facecolor=col)
        hist.set_ylabel('Number')
        hist.set_xlabel('Platescale')
        hist.grid(True)
        hist.set_title('Platescale for ' + str(comb).replace("'", '').replace("(", '').replace(")", '') +
                       ' (#' + str(int(n_tg)) + ')')
        plt.savefig(os.path.join(dir_out,
                                 'histogram_' + str(comb).replace("'", '').replace("(", '').replace(")", '').replace(",", "_") + '.pdf'))
        plt.close('all')
    # histogramplot for all combinations
    fig = plt.figure()
    hist_all = fig.add_subplot(111)
    hist_all = fig.add_subplot(111)
    hist_all.hist(astrometry['platescale'], 20, facecolor='red')
    hist_all.set_ylabel('Number')
    hist_all.set_xlabel('Platescale')
    hist_all.grid(True)
    hist_all.set_title('Platescale for all combinations ' +
                       '(#' + str(len(astrometry)) + ')')
    plt.savefig(os.path.join(dir_out, 'histogram_all.pdf'))
    plt.close('all')
