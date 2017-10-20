# a quick rapper to run all fields in the subdirs of the dir dir_sub

import os
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# returns to dir def in the end
def run(dir_def=os.getcwd(),
        dir_sub=os.getcwd()):

    # list all subdirs excluding files
    dirs_to_run = [d for d in sorted(os.listdir(dir_sub))
                   if os.path.isdir(os.path.join(dir_sub, d))]

    for idir, direct in enumerate(dirs_to_run):
        dir_cur = os.path.join(dir_sub, direct)
        os.chdir(dir_cur)
        from .platescale import do as do_platescale
        resdum = do_platescale(dir_data2=dir_cur)
        resdum['directory'] = direct
        if idir == 0:
            results = resdum
        else:
            results = results.append(resdum)

    os.chdir(dir_def)
    results = results.sort_values('mjd_obs')
    results = results.reset_index()
    results.to_csv(os.path.join(dir_def, 'results_subfolders.csv'),
                   index=False)

    # make some plots
    res_agpm = results.iloc[np.where(results.agpm == True)[0]]
    res_sat  = results.iloc[np.where(results.agpm == False)[0]]
    plt.errorbar(res_agpm.mjd_obs, res_agpm.rot_rob_weighted,
                 yerr=res_agpm.rot_rob_weighted_err,
                 label='AGPM', marker='x')
    plt.errorbar(res_sat.mjd_obs, res_sat.rot_rob_weighted,
                 yerr=res_sat.rot_rob_weighted_err,
                 label='saturated', marker='x')
    plt.title('True North (East of West)')
    plt.ylabel('True North [deg]')
    plt.xlabel('Modified Julian date')
    plt.legend()
    plt.savefig('true_north_rob_weighted.pdf', overwrite=True)
    plt.close('all')

    plt.errorbar(res_agpm.mjd_obs, res_agpm.pxscl_rob_weighted,
                 yerr=res_agpm.pxscl_rob_weighted_err,
                 label='AGPM', marker='x')
    plt.errorbar(res_sat.mjd_obs, res_sat.pxscl_rob_weighted,
                 yerr=res_sat.pxscl_rob_weighted_err,
                 label='saturated', marker='x')
    plt.title('Platescale')
    plt.ylabel('Platescale [mas/px]')
    plt.xlabel('Modified Julian date')
    plt.legend()
    plt.savefig(os.path.join(dir_def, 'pxscl_rob_weighted.pdf'),
                overwrite=True)
    plt.close('all')
    
    
    print('Done running all {} subfolders on {}. Saved to {}\
'.format(len(dirs_to_run),
         dir_sub,
         dir_def))
    
    import ipdb; ipdb.set_trace()

