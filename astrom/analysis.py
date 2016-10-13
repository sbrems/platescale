from __future__ import print_function,division
import os
import numpy as np
import pandas as pd
import misc
from parameters import *

def connect_sources(sex_coords,source,n_images=None):
    print('Connecting and verifying matches of ', n_images,' images')
    excluded = pd.DataFrame()
    for aa in range(n_images):
        #set values to nan and then to the value
        for val in ['x_cat_rot','y_cat_rot',star_id,'offset','ra_now','dec_now']:
                    sex_coords[aa][val] = [np.nan] * len(sex_coords[aa]['x_image'])
        for bb in range(len(sex_coords[aa]['x_image'])):
            for cc in range(len(source['x_data_rot'])):
                dist,dummy_angle = misc.dist_ang((sex_coords[aa]['y_image'][bb],sex_coords[aa]['x_image'][bb]),\
                                     (source['y_data_rot'][cc],source['x_data_rot'][cc]))
                if dist <= search_rad:
                    if not np.isfinite(sex_coords[aa]['x_cat_rot'][bb]): #new entry
                        sex_coords[aa]['x_cat_rot'][bb] = source['x_data_rot'][cc]
                        sex_coords[aa]['y_cat_rot'][bb] = source['y_data_rot'][cc]
                        sex_coords[aa]['offset'][bb]= dist
                        sex_coords[aa][star_id][bb] = source[star_id][cc]
                        sex_coords[aa]['ra_now'][bb]= source['ra_now'][cc]
                        sex_coords[aa]['dec_now'][bb]=source['dec_now'][cc]
                    elif dist < sex_coords[aa]['offset'][bb]: #if closer replace
                        sex_coords[aa]['x_cat_rot'][bb] = source['x_data_rot'][cc]
                        sex_coords[aa]['y_cat_rot'][bb] = source['y_data_rot'][cc]
                        sex_coords[aa]['offset'][bb]  = dist
                        sex_coords[aa][star_id][bb] = source[star_id][cc]
                        sex_coords[aa]['ra_now'][bb]= source['ra_now'][cc]
                        sex_coords[aa]['dec_now'][bb]=source['dec_now'][cc]
        #move stars where nothing in the catalogue was found in table excluded
        temp = sex_coords[aa][~np.isfinite(sex_coords[aa]['x_cat_rot'])]
        temp['image'] = aa
        excluded = excluded.append(temp,ignore_index=True)
        sex_coords[aa] = sex_coords[aa][np.isfinite(sex_coords[aa]['x_cat_rot'])]
        sex_coords[aa].reset_index(inplace=True,drop= True)
        #also check that each star is in image is matched to only one cat star. Always use closer ones
        #return the multiples
        multiples = sex_coords[aa][sex_coords[aa].duplicated(subset=('x_cat_rot','y_cat_rot'),keep=False)]
        multiples.reset_index(inplace=True,drop=True)
        if len(multiples) >= 1:
        #check if there are as many stars in catalogue around as duplicates 
            i_stars_around = [] #for stars around 
            dists = []
            #give the different distances to all valid catalogue points around
            for cc in range(len(source['x_data_rot'])):
                dist,dummy_angle = misc.dist_ang((multiples['y_cat_rot'][0],multiples['x_cat_rot'][0]),\
                                     (source['y_data_rot'][cc] ,source['x_data_rot'][cc] ))
                if dist <= search_rad:
                    i_stars_around.append(cc)
                    dists.append(dist)
            if (len(multiples) > len(i_stars_around)):
                print('ATTENTION!!\nThere are not enough stars in catalogue. Found ',\
                                 len(multiples),' : ','\n', multiples[['x_image','y_image']],\
                                 ' in image ',aa,' but only these in catalogue: ',[source[star_id][kk] for kk in i_stars_around],' \n Maybe',\
                                 ' sextractor split one source or search rad too small or',\
                                 ' there are different multiples!?\n Ignoring those.')
                excluded = excluded.append(sex_coords[aa].iloc[np.where(sex_coords[aa]['x_cat_rot'] == multiples['x_cat_rot'][0])],\
                                           ignore_index = True)
                excluded['image'].fillna(aa,inplace=True)
                sex_coords[aa] =sex_coords[aa].iloc[np.where(sex_coords[aa]['x_cat_rot'] != multiples['x_cat_rot'][0])]
                sex_coords[aa].reset_index(inplace=True,drop=True)
                excluded.reset_index(inplace=True,drop=True)

            elif len(multiples) >2:
                raise ValueError('More than 2 multiples (not implemented): ','\n', multiples[['x_image','y_image']],\
                                 ' in image ',aa,' but only these in catalogue: ',[source[star_id][kk] for kk in i_stars_around])
                
            else:#choose the closest ones in that sense that the maximum distance is minimal and each star has one counterp.
                dist_max = []
                raise ValueError('NOT IMPLEMENTED!!!MULTIPLE DETECTION!')
                #create all possible combinations for 2 same
                for icat1 in range(len(i_stars_around)):
                    for icat2 in range(len(i_stars_around)-1):
                        if icat1 != icat2:
                            dist_max.append([np.max(misc.dist_ang((multiples['y_image'][0],multiples['x_imag'][0]),
                                (source['y_data_rot'][i_stars_around[icat1]],source['x_data_rot'][i_stars_around[icat1]]))[0],
                                misc.dist_ang((multiples['y_image'][1],multiples['x_image'][1]),
                                              (source['y_data_rot'][i_stars_around[icat2]],source['x_data_rot'][i_stars_around[icat2]]))[0]),icat1,icat2])
                
                        ipdb.set_trace()
        
        #sort alphabetically so we always get the same star as primary
        sex_coords[aa] = sex_coords[aa].sort_values(['MAIN_ID']).reset_index()
    
    return sex_coords,excluded,multiples
