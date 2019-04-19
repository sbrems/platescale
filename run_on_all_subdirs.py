# a quick rapper to run all fields in the subdirs of the dir dir_sub

import os
import pandas as pd
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
from platescale import do as do_platescale
from platescale import TooFewStarsFoundError
from my_stats import weighted_avg_and_std
from tqdm import tqdm
import warnings

import ipdb


# returns to dir def in the end
def run(dir_def=os.getcwd(),
        dir_sub=os.getcwd(),
        only_load_savedata=False):
    # list all subdirs excluding files
    dirs_to_run = [d for d in sorted(os.listdir(dir_sub))
                   if os.path.isdir(os.path.join(dir_sub, d))]

    # 2017-05-19T10:02:17.135 (mjd = 57892.41825387732) is a good example
    # the bad obs will not be plotted later
    bad_obs = [
        58057.99636672,  # 2017-10-31T23:54:46.085
        58057.99636672,  # 2017-10-31T23:55:15.808
        58058.00108815,   # 2017-11-01T00:01:34.016
        58057.99602724,  # 2017-10-31T23:54:16.754
        58058.00075385,  # 2017-11-01T00:01:05.133
        58377.3987758,  # 2018-09-16T09:34:14.229
        58057.00062587,  # 2017-10-31T00:00:54.075
        57995.426999,  # 2017-08-30T10:14:52.681
        57995.427516,  # 2017-08-30T10:15:37.413
        57995.427516,  # 2017-08-30T10:16:20.747
        58377.3987758, # 2018-09-16T09:34:14.229
        57995.426497,  # 2017-08-30T10:14:09.349
    ]
    dummy = [
        # bad results
        58276.41307139,
        58410.18972717,
        58410.18972717,
        58317.415908,   # 2018-07-18 09:58:54.419
        58276.413759,  # 2018-07-18 09:58:54.451
        57948.43758987,  # 2017-07-14T10:30:07.765

        # bad images
        57701.00367119213,  # 2016-11-09T00:05:17.191 # very bad
        #57830.02246298611,  # 2017-03-18T00:32:20.802 # elongated
        #57830.02287944444,  # 2017-03-18T00:32:56.784 # elongated
        #57830.024864212966,  # 2017-03-18T00:35:48.268 # elongated
        #57830.02528041667,  # 2017-03-18T00:36:24.228 # elongated
        #57830.03030775463,  # 2017-03-18T00:43:38.590 # elongated
        57948.4366744213,  # 2017-07-14T10:28:48.670
        # 57995.425997893515 - 57995.42801790509   # 2017-08-30T10:13:26.218-2017-08-30T10:16:20.747  # have problems with stripes and corners
        58057.995355,  #   2017-10-31T23:53:18.698
        58057.99671074,  # 2017-10-31T23:55:15.808
        58058.01508667824, # 2017-11-01T00:21:43.489  # heavilly elongated
        58058.01555386574, # 2017-11-01T00:22:23.854  #heavilly elongated
        58276.41341312,  # 2018-06-07T09:55:18.894
        58317.4155746,  # 2018-07-18T09:58:25.645
        58376.38211902778, # 2018-09-15T09:10:15.084
        58376.38245853009, # 2018-09-15T09:10:44.417
        58376.3827952662,  # 2018-09-15T09:11:13.511
        58410.19006140046,  # 2018-10-19T04:33:41.305
        58410.189727164354,  # 2018-10-19T04:33:12.427
        58410.18938724537,  # 2018-10-19T04:32:43.058
        58449.05086319445,  # 2018-11-27T01:13:14.580 # stripes in top right+agpm
        58449.05119777778,  # 2018-11-27T01:13:43.488 # stripes in top right+agpm
        58276.41375933,  # 2018-06-07T09:55:48.778
        58317.41590763,  # 2018-07-18T09:58:54.419
        58449.05153923,  # 2018-11-27T01:14:12.989
        58449.05188973,  # 2018-11-27T01:14:43.273
    ]

    if not only_load_savedata:
        failed = []
        warnings.filterwarnings("ignore")
        for idir, direct in enumerate(tqdm(dirs_to_run)):
            dir_cur = os.path.join(dir_sub, direct)
            os.chdir(dir_cur)
            try:
                resdum = do_platescale(dir_data2=dir_cur)
            except TooFewStarsFoundError:
                failed.append(direct)
            resdum['directory'] = direct
            if idir == 0:
                results = resdum
            else:
                results = results.append(resdum)
        results = results.reset_index()

        os.chdir(dir_def)
        results = results.sort_values('mjd_obs')
        results = results.reset_index()
        results.to_csv(os.path.join(dir_sub, 'results_subfolders.csv'),
                       index=False)
        with open(os.path.join(dir_sub, 'failed_folders.txt'), 'w+') as f:
            for item in failed:
                f.write("{}\n".format(item))
            f.write('####FILES IGNORED FOR PLOTTING###')
            for item in bad_obs:
                tt = Time(item, format='mjd')
                f.write('{}; {}'.format(tt.mjd, tt.isot))
    else:  # load data
        results = pd.read_csv(os.path.join(dir_sub, 'results_subfolders.csv'))

    # make some plots


    sns.set_style('darkgrid')
    sns.set_context('paper')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    plt.close('all')
    list_valid = [~np.any(np.isclose(resi, bad_obs, rtol=1e-20, atol=1e-5)) for resi in results['mjd_obs']]
    print('Excluded obs: {} (of {}, should be {} discarded, but maybe some are in failed)'.format(
        np.sum(np.invert(list_valid)),len(results), len(bad_obs)))

    res_agpm = results[(results['agpm'] == True) & (list_valid)]
    res_sat  = results[(results['agpm'] == False) & (list_valid)]

    agpm_mean_std = weighted_avg_and_std(
        res_agpm.pxscl_rob_weighted,
        (1./res_agpm.pxscl_rob_weighted_err)**2)
    sat_mean_std = weighted_avg_and_std(
        res_sat.pxscl_rob_weighted,
        (1./res_sat.pxscl_rob_weighted_err)**2)

    agpm_rot_mean_std = weighted_avg_and_std(
        res_agpm.rot_rob_weighted,
        (1./res_agpm.rot_rob_weighted_err)**2)
    sat_rot_mean_std = weighted_avg_and_std(
        res_sat.rot_rob_weighted,
        (1./res_sat.rot_rob_weighted_err)**2)
    
    plt.errorbar(res_agpm.mjd_obs, res_agpm.rot_rob_weighted,
                 yerr=res_agpm.rot_rob_weighted_err,
                 label=r'AGPM (${:.2f}\pm{:.2f})^\circ$'.format(
                     agpm_rot_mean_std[0], agpm_rot_mean_std[1]),
                 capsize=0, fmt='x', color='b')
    plt.errorbar(res_sat.mjd_obs, res_sat.rot_rob_weighted,
                 yerr=res_sat.rot_rob_weighted_err,
                 label='satPSF (${:.2f}\pm{:.2f})^\circ$'.format(
                     sat_rot_mean_std[0], sat_rot_mean_std[1]),
                 fmt='x', capsize=0, color='g')

    plt.axhline(agpm_rot_mean_std[0], color='b')
    plt.axhspan(agpm_rot_mean_std[0] - agpm_rot_mean_std[1],
                agpm_rot_mean_std[0] + agpm_rot_mean_std[1],
                color='b', alpha=.3)
    plt.axhline(sat_rot_mean_std[0], color='g', )
    plt.axhspan(sat_rot_mean_std[0] - sat_rot_mean_std[1],
                sat_rot_mean_std[0] + sat_rot_mean_std[1],
                color='g', alpha=.3)

    
    # plt.title('True North (East of West)')
    plt.ylabel('True North [deg]')
    plt.xlabel('Modified Julian date')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(dir_sub, 'true_north_rob_weighted.pdf'),
                pad_inches=0.,
                overwrite=True)
    plt.close('all')

    plt.errorbar(res_agpm.mjd_obs, res_agpm.pxscl_rob_weighted,
                 yerr=res_agpm.pxscl_rob_weighted_err,
                 label=r'AGPM (${:.3f}\pm{:.3f}$) mas/px'.format(
                     agpm_mean_std[0], agpm_mean_std[1]),
                 fmt='x', color='b', capsize=0)
    plt.errorbar(res_sat.mjd_obs, res_sat.pxscl_rob_weighted,
                 yerr=res_sat.pxscl_rob_weighted_err,
                 label=r'satPSF (${:.3f}\pm{:.3f})$ mas/px'.format(
                     sat_mean_std[0], sat_mean_std[1]),
                 fmt='x', color='g', capsize=0)
    
    plt.axhline(agpm_mean_std[0], color='b')
    plt.axhspan(agpm_mean_std[0] - agpm_mean_std[1],
                agpm_mean_std[0] + agpm_mean_std[1],
                color='b', alpha=.3)
    plt.axhline(sat_mean_std[0], color='g', )
    plt.axhspan(sat_mean_std[0] - sat_mean_std[1],
                sat_mean_std[0] + sat_mean_std[1],
                color='g', alpha=.3)
    
    # plt.title('Platescale')
    plt.ylabel('Platescale [mas/px]')
    plt.xlabel('Modified Julian date')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(dir_sub, 'pxscl_rob_weighted.pdf'),
                pad_inches=0.,
                overwrite=True)
    plt.close('all')

    plt.scatter(res_agpm.rot_rob_weighted, res_agpm.rot_header,
                label='AGPM', color='g')
    plt.scatter(res_sat.rot_rob_weighted, res_sat.rot_header,
                label='satPSF', color='b')

    plt.axhline(sat_mean_std[0])
    plt.axhspan(sat_mean_std[0] - sat_mean_std[1],
                sat_mean_std[0] + sat_mean_std[1],)

    # plt.title('True north vs. header rotation')
    plt.xlabel('True north offset found [deg]')
    plt.ylabel('Rotation given in header [deg]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(dir_sub, 'header_rot_vs_true_north.pdf'),
                pad_inches=0.,
                overwrite=True)
    plt.close('all')
    
    print('Done running all {} subfolders on {} where results are saved'.format(
        len(dirs_to_run),
         dir_sub,))

    import ipdb; ipdb.set_trace()
