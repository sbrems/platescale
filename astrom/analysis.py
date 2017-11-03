# import os
import numpy as np
import pandas as pd
from . import misc
import itertools
from .parameters import *


def compare_to_med(sex_coords, source_med, ignored_med, verbose=True):
    '''Compare the stellar positions to the positions in the median image.
    If no correspondance in median image found => discard image (prob bad quality)
    If multiple are found: error. This shouldnt happen. make max mvmnt smaller.
    Return the connected properties so statistics can be made.'''
    n_images = len(sex_coords)
    if verbose:
        print('Connecting stars to stars in median image for ', n_images, ' images.')
    ign_im = []
    ignored = []

    for i_im in range(n_images):
        idz_src = np.array(sex_coords[i_im].index.values)
        idz_med = np.array(source_med.index.values)
        idz_ign = np.array(ignored_med.index.values)
        sex2med = {}
        ign_in_this_im = []

        # use the unshifted images as there are no fine shifted for source_med (yet)
        for i_src in idz_src:
            idx_med = []
            idx_ign = []
            x_src = sex_coords[i_im]['x_image'][i_src]
            y_src = sex_coords[i_im]['y_image'][i_src]
            # for i_med,xy_med in zip(idz_med,source_med[['x_image','y_image']].as_matrix()):
            for i_med in idz_med:
                #x_med,y_med = xy_med
                x_med, y_med = source_med['x_image'][i_med], source_med['y_image'][i_med]
                dist = np.sqrt((x_src - x_med)**2 + (y_src - y_med)**2)
                if dist <= max_mvmnt:
                    idx_med.append(i_med)
            # check if in ignored sources
            if len(idx_med) == 0:
                for i_ign, xy_ign in zip(ignored_med.index.values, ignored_med[['x_image', 'y_image']].as_matrix()):
                    x_ign, y_ign = xy_ign
                    dist = np.sqrt((x_src - x_ign)**2 + (y_src - y_ign)**2)
                    if dist <= max_mvmnt:
                        idx_ign.append(i_ign)
                if len(idx_ign) == 1:
                    # an ignored source found. everything ok.
                    ign_in_this_im.append(idx_ign[0])
                    idz_ign = np.delete(idz_ign, np.where(idz_ign == idx_ign))
                elif len(idx_ign) > 1:
                    # multiples found? look if also multiples found in image
                    raise ValueError('There were multiple ignored sources found in image_nr/src_nr: ',
                                     i_im, i_src)
                else:  # len(idx_ign) ==0
                    # idx_ign == 0:still no source was found => ignore the image. Source was prob split
                    if verbose:
                        print('No corresponding source found in median image. Maybe source %d (%d,%d)was split? Ignoring image nr %d.' % (
                            i_src, y_src, y_src, i_im))
                    ign_im.append(i_im)
                    break

            # now differentiate the different cases in the found sources:
            elif len(idx_med) > 1:  # multiple sources found => See if there are also multiple in image
                multi_im = []
                for ii_src in idz_src:
                    x_src2, y_src2 = sex_coords[i_im][[
                        'x_image', 'y_image']].loc[ii_src]
                    dist = np.sqrt((x_src - x_src2)**2 + (y_src - y_src2)**2)
                    if dist <= search_rad:
                        multi_im.append(ii_src)
                if len(multi_im) != len(idx_med):
                    if len(multi_im) == 1 and len(idx_med) == 2:
                        # assume the brighter source in the cat is the one found
                        mag0 = source_med['mag_auto'].loc[idx_med[0]]
                        mag1 = source_med['mag_auto'].loc[idx_med[1]]
                        idx_med_use = -999
                        if mag0 < (mag1 - min_contrast):  # 0 brighter
                            idx_med_use = idx_med[0]
                            idx_med_notuse = idx_med[1]
                        elif mag1 < (mag0 - min_contrast):  # 1 brighter
                            idx_med_use = idx_med[1]
                            idx_med_notuse = idx_med[0]
                        if idx_med_use != -999:
                            sex2med[multi_im[0]] = idx_med_use
                            # remove both found indices
                            if verbose:
                                print('Found sources %s in image nr %d but %s in median image. As %s is brighter (mag %d vs %d in median image), assuming this as true source.' % (
                                    multi_im, i_im, idx_med, idx_med_use, mag0, mag1))

                            idz_src = [
                                x for x in idz_src if x not in list(sex2med.keys())]
                            idz_med = [x for x in idz_med if x not in idx_med]

                            # break
                    else:
                        if verbose:
                            print('Only sources %s found in image nr %d, but %s in median image. Ignoring image.' % (
                                multi_im, i_im, idx_med))
                        ign_im.append(i_im)
                        break
                    # raise ValueError('There were multiple sources found in image_nr/src_nr: ',
                    #                 i_im,i_src)
                else:  # len(multi_im) == len(idx_med) >1v:#good case. minimize distance
                    list_im = np.hstack((np.transpose([multi_im, ]),
                                         sex_coords[i_im].loc[multi_im].as_matrix(columns=['x_image', 'y_image'])))
                    list_med = np.hstack((np.transpose([idx_med, ]),
                                          source_med.loc[idx_med].as_matrix(columns=['x_image', 'y_image'])))
                    dict_match = match_smallest_dist(list_im, list_med)
                    sex2med.update(dict_match)
                    # remove the found indizes from the lists
                    idz_src = [
                        x for x in idz_src if x not in list(sex2med.keys())]
                    idz_med = [x for x in idz_med if x not in list(
                        sex2med.values())]
            elif len(idx_med) == 1:  # the good case
                sex2med[int(i_src)] = int(idx_med[0])
                idz_med = np.delete(idz_med, np.where(idz_med == idx_med))
            else:  # alrdy have all cases (I hope)
                raise ValueError('You shouldnt be here...bug in code???')
        # now save the results for this image

        if i_im not in ign_im:
            ignored.append([i_im, ign_in_this_im])

            sex_dum1 = source_med.loc[list(sex2med.values())]
            sex_dum2 = sex_coords[i_im].loc[list(sex2med.keys())]
            sex_dum2 = sex_dum2.set_index(
                [[sex2med[x] for x in sex_dum2.index]])
            sex_coords[i_im] = sex_dum1.join(
                sex_dum2, how='inner', lsuffix='_in_med')
    excluded = [sex_coords[x] for x in ign_im]

    sex_coords = [sex_coords[x]
                  for x in range(len(sex_coords)) if x not in ign_im]

    return sex_coords, excluded, ignored, ign_im


def connect_sources(sex_coords, source_med, n_images=None, verbose=True):
    if verbose:
        print('Connecting and verifying matches of ', n_images, ' images')
    excluded = pd.DataFrame()
    for aa in range(n_images):
        # set values to nan and then to the value
        for val in ['x_data_rot', 'y_data_rot', star_id,
                    'offset', 'ra_now', 'dec_now']:
            sex_coords[aa][val] = [np.nan] * len(sex_coords[aa]['x_image'])
        for bb in range(len(sex_coords[aa]['x_image'])):
            for cc in range(len(source['x_data_rot'])):
                dist, dummy_angle = misc.dist_ang((sex_coords[aa]['y_image'][bb], sex_coords[aa]['x_image'][bb]),
                                                  (source['y_data_rot'][cc], source['x_data_rot'][cc]))
                if dist <= search_rad:
                    # new entry
                    if not np.isfinite(sex_coords[aa]['x_data_rot'][bb]):
                        sex_coords[aa]['x_data_rot'][bb] = source['x_data_rot'][cc]
                        sex_coords[aa]['y_data_rot'][bb] = source['y_data_rot'][cc]
                        sex_coords[aa]['offset'][bb] = dist
                        sex_coords[aa][star_id][bb] = source[star_id][cc]
                        sex_coords[aa]['ra_now'][bb] = source['ra_now'][cc]
                        sex_coords[aa]['dec_now'][bb] = source['dec_now'][cc]
                    elif dist < sex_coords[aa]['offset'][bb]:  # if closer replace
                        sex_coords[aa]['x_data_rot'][bb] = source['x_data_rot'][cc]
                        sex_coords[aa]['y_data_rot'][bb] = source['y_data_rot'][cc]
                        sex_coords[aa]['offset'][bb] = dist
                        sex_coords[aa][star_id][bb] = source[star_id][cc]
                        sex_coords[aa]['ra_now'][bb] = source['ra_now'][cc]
                        sex_coords[aa]['dec_now'][bb] = source['dec_now'][cc]
        # move stars where nothing in the catalogue was found in table excluded
        temp = sex_coords[aa][~np.isfinite(sex_coords[aa]['x_data_rot'])]
        temp['image'] = aa
        excluded = excluded.append(temp, ignore_index=True)
        sex_coords[aa] = sex_coords[aa][np.isfinite(
            sex_coords[aa]['x_data_rot'])]
        sex_coords[aa].reset_index(inplace=True, drop=True)
        # also check that each star is in image is matched to only one cat star. Always use closer ones
        # return the multiples
        multiples = sex_coords[aa][sex_coords[aa].duplicated(
            subset=('x_data_rot', 'y_data_rot'), keep=False)]
        multiples.reset_index(inplace=True, drop=True)
        if len(multiples) >= 1:
            # check if there are as many stars in catalogue around as duplicates
            i_stars_around = []  # for stars around
            dists = []
            # give the different distances to all valid catalogue points around
            for cc in range(len(source['x_data_rot'])):
                dist, dummy_angle = misc.dist_ang((multiples['y_data_rot'][0], multiples['x_data_rot'][0]),
                                                  (source['y_data_rot'][cc], source['x_data_rot'][cc]))
                if dist <= search_rad:
                    i_stars_around.append(cc)
                    dists.append(dist)
            if (len(multiples) > len(i_stars_around)):
                if verbose:
                    print('ATTENTION!!\nThere are not enough stars in catalogue. Found ',
                          len(multiples), ' : ', '\n', multiples[[
                              'x_image', 'y_image']],
                          ' in image ', aa, ' but only these in catalogue: ', [
                              source[star_id][kk] for kk in i_stars_around], ' \n Maybe',
                          ' sextractor split one source or search rad too small or',
                          ' there are different multiples!?\n Ignoring those.')
                excluded = excluded.append(sex_coords[aa].iloc[np.where(sex_coords[aa]['x_data_rot'] == multiples['x_data_rot'][0])],
                                           ignore_index=True)
                excluded['image'].fillna(aa, inplace=True)
                sex_coords[aa] = sex_coords[aa].iloc[np.where(
                    sex_coords[aa]['x_data_rot'] != multiples['x_data_rot'][0])]
                sex_coords[aa].reset_index(inplace=True, drop=True)
                excluded.reset_index(inplace=True, drop=True)

            elif len(multiples) > 2:
                raise ValueError('More than 2 multiples (not implemented): ', '\n', multiples[['x_image', 'y_image']],
                                 ' in image ', aa, ' but only these in catalogue: ', [source[star_id][kk] for kk in i_stars_around])

            else:  # choose the closest ones in that sense that the maximum distance is minimal and each star has one counterp.
                dist_max = []
                raise ValueError('NOT IMPLEMENTED!!!MULTIPLE DETECTION!')
                # create all possible combinations for 2 same
                for icat1 in range(len(i_stars_around)):
                    for icat2 in range(len(i_stars_around) - 1):
                        if icat1 != icat2:
                            dist_max.append([np.max(misc.dist_ang((multiples['y_image'][0], multiples['x_imag'][0]),
                                                                  (source['y_data_rot'][i_stars_around[icat1]], source['x_data_rot'][i_stars_around[icat1]]))[0],
                                                    misc.dist_ang((multiples['y_image'][1], multiples['x_image'][1]),
                                                                  (source['y_data_rot'][i_stars_around[icat2]], source['x_data_rot'][i_stars_around[icat2]]))[0]), icat1, icat2])

        # sort alphabetically so we always get the same star as primary
        sex_coords[aa] = sex_coords[aa].sort_values(['MAIN_ID']).reset_index()

    return sex_coords, excluded, multiples


def connect_median_sources(source_cat, median_sources, verbose=True):
    '''Find the corresponding source in the catalog to each source found in the median image 
    and only return those. Do this by checking in an area and then
    making all possible source-source combinations
    and find the one with the minimum distance between all that were found in
    that area. Also return the ones where
    no counterpart in catalogue was found.'''

    if verbose:
        print('Finding the sources in the catalogue')

    idz_cat = np.array(source_cat.index.values)
    idz_med = np.array(median_sources.index.values)
    med_ignore = []  # for sources found in image but not in catologue
    med2cat = {}
    for i_med in idz_med:
        # check this for the case of multiples
        if i_med not in list(med2cat.keys()):
            id_cat = []
            for i_cat in idz_cat:
                dist = np.sqrt((source_cat.x_data_rot[i_cat] - median_sources.x_image[i_med])**2 +
                               (source_cat.y_data_rot[i_cat] - median_sources.y_image[i_med])**2)
                if dist < search_rad:
                    id_cat.append(i_cat)

            if len(id_cat) == 0:  # no source in cat found:
                if verbose:
                    print('WARNING!There has been no source in the catalogue nearby ' +
                          'the one found at coordinates x,y:', median_sources.x_image[i_med],
                          median_sources.y_image[i_med])
                    print('Ignoring this source')
                med_ignore.append(i_med)
                idz_med = np.delete(idz_med, np.where(idz_med == i_med))

            elif len(id_cat) == 1:  # one source found
                med2cat[i_med] = id_cat[0]
                # remove the found indizes from the lists
                idz_cat = [x for x in idz_cat if x not in list(
                    med2cat.values())]
                idz_med = [x for x in idz_med if x not in list(med2cat.keys())]
            elif len(id_cat) >= 2:  # multiple sources found in area
                # if there is only one source nearby in the data, ignore it
                multi_med = []
                for ii_med in idz_med:
                    dist = np.sqrt((median_sources.x_image[ii_med] - median_sources.x_image[i_med])**2 +
                                   (median_sources.y_image[ii_med] - median_sources.y_image[i_med])**2)
                    if dist < search_rad:
                        multi_med.append(ii_med)
                if len(multi_med) < len(id_cat):
                    # kick out the faintest ones, if the contrast is high enough
                    magorder = source_cat.iloc[id_cat]['MAG_USE'].sort_values(
                    ).index
                    magorder = magorder[:len(multi_med)]
                    for magidx in range(len(magorder) - 1):
                        if -(source_cat.iloc[magorder[magidx]]['MAG_USE'] -
                             source_cat.iloc[magorder[magidx + 1]]['MAG_USE']) < min_contrast:
                            raise ValueError('There are multiple equally bright sources in the catalogue at the corresponding ' +
                                             'position in median position. Please check from hand which is the right one and ' +
                                             'remove the other from the catalogue.x/y med/cat:', median_sources.x_image[
                                                 multi_med],
                                             median_sources.y_image[multi_med],
                                             source_cat.x_data_rot[id_cat],
                                             source_cat.y_data_rot[id_cat])

                    # The sources brightness is sufficiently different.
                    # Use brightest ones
                    if verbose:
                        print('Removed faint sources in same area as contrast was \
high enough')
                    id_cat = magorder

                # same number of sources...minimize their distances
                if len(multi_med) == len(id_cat):
                    list_med = np.hstack((np.transpose([multi_med, ]),
                                          median_sources.loc[multi_med].as_matrix(columns=['x_image', 'y_image'])))
                    list_cat = np.hstack((np.transpose([id_cat, ]),
                                          source_cat.loc[id_cat].as_matrix(columns=['x_data_rot', 'y_data_rot'])))
                    dict_match = match_smallest_dist(list_med, list_cat)
                    med2cat.update(dict_match)
                    # remove the found indizes from the lists
                    idz_cat = [x for x in idz_cat if x not in list(
                        med2cat.values())]
                    idz_med = [
                        x for x in idz_med if x not in list(med2cat.keys())]
                else:  # len(multi_med) > len(id_cat)
                    print('There were sources in the catalogue found. \
However, not as many as there were sources found in the median image \
 ({} vs {}). Ignoring them.'.format(
                        len(id_cat), len(multi_med)),
                        median_sources.x_image[multi_med],
                        median_sources.y_image[multi_med],
                        source_cat.x_data_rot[id_cat],
                        source_cat.y_data_rot[id_cat])
                    for imulti_med in multi_med:
                        med_ignore.append(imulti_med)
                        idz_med = np.delete(idz_med, np.where(idz_med == imulti_med))
                    for iid_cat in id_cat:
                        idz_cat = np.delete(idz_cat, np.where(idz_cat == iid_cat))

    # make the catalogue with only the sources found in the median image
    source_med_dum1 = source_cat.loc[list(med2cat.values())]
    source_med_dum2 = median_sources.loc[list(med2cat.keys())]
    source_med_dum2 = source_med_dum2.set_index(
        [[med2cat[x] for x in source_med_dum2.index]])
    source_med = source_med_dum1.join(source_med_dum2, how='inner')
    ignored_med = median_sources.loc[med_ignore]
    return source_med, ignored_med


def match_smallest_dist(list1, list2):
    '''Give to lists with 3xn dim: index, x_coord,y_coord.
    This code returns a dictionary converting from list 1 index
    to list 2 index in that way, that the total connection distances
    are minimized. Each point is only used once.len(list2)!>=len(list1)'''
    # #make a list with all possible combinations. a list with lists of all combinations.
    list1 = np.array(list1)
    list2 = np.array(list2)
    l12l2 = {}
    combinations = [list(zip(x, list1))
                    for x in itertools.permutations(list2, len(list1))]
    comb_len = []
    for comb in combinations:
        tot_length = 0.
        for connection in comb:
            length = np.sqrt((connection[0][1] - connection[1][1])**2 +
                             (connection[0][2] - connection[1][2])**2)
            tot_length += length
        comb_len.append(tot_length)
    # now choose the minimum and get the index (second min, should be only one value)
    i_min_comb = np.min(np.where(comb_len == np.min(comb_len)))
    # add this to the dict
    for connection in combinations[i_min_comb]:
        l12l2[int(connection[1][0])] = int(connection[0][0])

    return l12l2
