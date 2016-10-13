from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import sigmaclip
from parameters import *
import misc

def single_matches(data_cube,sex_coords,n_images=None):
    '''plotting the single matches and doing some sorting in the pd.DataFrames (the astrometry one)'''
   
    print('Plotting and solving platescale for', n_images,' images')
    entries = [star_id+'1',star_id+'2','dis_cat_mas','dis_meas','shift_error','ang_cat','ang_meas','platescale','ang_diff','image_nr'] 
    astrometry = pd.DataFrame(columns=entries)
    for i_im in range(n_images):
        fig,im_matched = plt.subplots()
        im_matched.imshow(np.sqrt(data_cube[i_im,::]))
        im_matched.set_xlim((0,data_cube.shape[2]))
        im_matched.set_ylim((0,data_cube.shape[1]))
        #now calc the distances
        for first in range(len(sex_coords[i_im])):
            for second in range(first+1,len(sex_coords[i_im])):
                entry=[]
                entry.append(sex_coords[i_im][star_id][first])
                entry.append(sex_coords[i_im][star_id][second])
                dist_cat_deg,ang_cat_deg = misc.dist_ang((sex_coords[i_im]['ra_now'][first],\
                                                          sex_coords[i_im]['dec_now'][first]),\
                                                         (sex_coords[i_im]['ra_now'][second],\
                                                          sex_coords[i_im]['dec_now'][second]),
                                                         spherical=True)
                dist_meas,ang_meas= misc.dist_ang((-sex_coords[i_im]['x_image'][first],
                                                   sex_coords[i_im]['y_image'][first]),
                                                  (-sex_coords[i_im]['x_image'][second],
                                                   sex_coords[i_im]['y_image'][second]))
                #- needed to get east of north and not west of north
                ang_meas = np.rad2deg(ang_meas)
                shift_error = np.sqrt(sex_coords[i_im].shift_error[first]**2+
                                      sex_coords[i_im].shift_error[second]**2)
                #- for different definition of angles
                #ang_meas = (-1*ang_meas)#%360. dont do modulo to keep direction info4stead 359
                ang_cat_deg = ang_cat_deg#%360.
                entry.append(dist_cat_deg)
                entry.append(dist_meas)
                entry.append(shift_error)
                entry.append(ang_cat_deg)
                entry.append(ang_meas)
                entry.append(dist_cat_deg/dist_meas)
                #note: ang_diff contains correction for initial guess.
                ang_diff = (ang_meas - ang_cat_deg)%360
                if ang_diff > 180: ang_diff -=360 #make it betwenn -+180
                entry.append(ang_diff) 
                entry.append(i_im)
                astrometry = astrometry.append(pd.DataFrame([entry],columns=entries),ignore_index=True)
                astrometry['weight'] = astrometry.dis_meas / astrometry.shift_error
                #plot the platescale measured for each comb
                im_matched.plot([sex_coords[i_im]['x_image'][first],sex_coords[i_im]['x_image'][second]],\
                                [sex_coords[i_im]['y_image'][first],sex_coords[i_im]['y_image'][second]],\
                                'k--')
                im_matched.annotate(str('{:.2f}'.format(astrometry.iloc[-1]['platescale']))+'\n'+\
                                    str('{:.1f}'.format(astrometry.iloc[-1]['ang_diff'])),
                                    xy=( (sex_coords[i_im]['x_image'][first]+sex_coords[i_im]['x_image'][second])/2,\
                                         (sex_coords[i_im]['y_image'][first]+sex_coords[i_im]['y_image'][second])/2),\
                                    size=3,color='k')
            #continue plotting
            scat_cat = im_matched.scatter(sex_coords[i_im]['x_cat_rot'][first],sex_coords[i_im]['y_cat_rot'][first],\
                                          edgecolors ='red', facecolors='none',marker='o',s=6)
            scat_sex = im_matched.scatter(sex_coords[i_im]['x_image'][first],sex_coords[i_im]['y_image'][first],\
                                          edgecolors='black',facecolors='none',marker='o',s=6)
            im_matched.plot([sex_coords[i_im]['x_cat_rot'][first],sex_coords[i_im]['x_image'][first]],\
                            [sex_coords[i_im]['y_cat_rot'][first],sex_coords[i_im]['y_image'][first]],\
                            'r--')
            im_matched.annotate(sex_coords[i_im][star_id][first],xy=(sex_coords[i_im]['x_image'][first],\
                                                                     sex_coords[i_im]['y_image'][first]),\
                                size=6,color='red')
            
        im_matched.legend((scat_cat,scat_sex),('Catalogue','Found position'))
        plt.savefig('found_sources_'+'{:03}'.format(i_im)+'.svg')
        plt.close('all')

    return astrometry



def multiple_matches(astrometry,sigma_outliers):
    '''Give it the single matches and the sigma for clipping. This algorithm combines 
    the different measurements between the same targets and also does sigma clipping.
    Results w/ and w/o sigma clipping are returned.
    Also the final pxscl and rotation is stored and returned in results'''

    results_entries = ['pxscl_all','pxscl_all_err','rot_all','rot_all_err','pxscl_robust',
                       'pxscl_robust_err', 'rot_robust','rot_robust_err','pxscl_rob_weighted',
                       'rot_rob_weighted']
    results = pd.DataFrame([len(results_entries)*[np.nan]],columns=results_entries)
    
    #calculate the results
    results['pxscl_all']     = np.mean(astrometry['platescale'])
    results['pxscl_all_err'] = np.std(astrometry['platescale'])
    results['rot_all']       = np.mean(astrometry['ang_diff'])
    results['rot_all_err']   = np.std(astrometry['ang_diff'])

    #Do the robust calculation. I.e kick out all outliers
    rob_idz = astrometry[(astrometry.platescale <= np.median(astrometry.platescale)+
                                 sigma_outliers*np.std(astrometry.platescale)) & 
                                (astrometry.platescale <= np.median(astrometry.platescale)+
                                 sigma_outliers*np.std(astrometry.platescale))].index
    rob_idz_rot = astrometry[(astrometry.ang_diff <= np.median(astrometry.ang_diff)+
                                 sigma_outliers*np.std(astrometry.ang_diff)) & 
                                (astrometry.ang_diff <= np.median(astrometry.ang_diff)+
                                 sigma_outliers*np.std(astrometry.ang_diff))].index
    results['pxscl_robust']     = np.mean(astrometry.platescale[rob_idz])
    results['pxscl_robust_err'] = np.std(astrometry.platescale[rob_idz])
    results['rot_robust']       = np.mean(astrometry.ang_diff[rob_idz_rot])
    results['rot_robust_err']   = np.std(astrometry.ang_diff[rob_idz_rot])
    
    
    #do the same weighted
    results['pxscl_rob_weighted'],results['pxscl_rob_weighted_err'] =\
            misc.weighted_avg_and_std(astrometry.platescale[rob_idz],
                                      weights=astrometry.weight[rob_idz])
    results['rot_rob_weighted'],results['rot_rob_weighted_err'] =\
            misc.weighted_avg_and_std(astrometry.ang_diff[rob_idz_rot],
                                      weights=astrometry.weight[rob_idz_rot])

 

    
    #clip the means of the platescale and then ignore the ones above and below               

    astrometry_grouped = astrometry.groupby([star_id+'1',star_id+'2']).mean().add_suffix('_mean').reset_index()
    astrometry_grouped_cliped = astrometry_grouped[abs(astrometry_grouped['platescale_mean'] -np.median(astrometry_grouped['platescale_mean']))\
                             <= (sigma_outliers * np.std(astrometry_grouped['platescale_mean']))]
    

    return astrometry_grouped, astrometry_grouped_cliped,results