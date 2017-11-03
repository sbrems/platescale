from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy. units as u
import numpy as np
from starclass import Star
import os
import ipdb


###############################################


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx
####################################################


pfnametcolmag = 'BB_HSTcolorTempconverter.csv'
tfn = 'NGC0104c.pm'
vosalim = 900000  # Split by this number. Vosa can handle only
# ~10000 at  a time
pmtable = Table.read(tfn, format='ascii.csv',
                     delimiter='\s', comment='#')
# apply suggested corrections
zerocoord = SkyCoord(['00:24:05.71 -72:04:52.70'],
                     unit=[u.hourangle, u.deg], frame='icrs')
pmtable['DEC_cor'] = ((pmtable['y_M'] - pmtable['Delta_y'] - 4950.07) *
                      40 * u.mas + zerocoord.dec).to('deg')
pmtable['RA_cor'] = (-((pmtable['x_M'] - pmtable['Delta_x'] - 4950.79) *
                       40 * u.mas * u.deg /
                       np.cos(np.deg2rad(pmtable['DEC_cor']))) +
                     zerocoord.ra).to('deg')
pmtable['RA_err'] = [0, ] * len(pmtable) * u.deg
pmtable['DEC_err'] = [0, ] * len(pmtable) * u.deg
pmtable['MAG_EXTRAPOL'] = 15 * (pmtable['m_F814W'] - pmtable['m_F606W']) +\
    pmtable['m_F814W']
pmtable['MAG_EXTRAPOL_err'] = np.sqrt(pmtable['rms_F606W(mag)']**2 +
                                      pmtable['rms_F814W(mag)']**2)
##########################################
print('Fitting BB to all the stars')
# make a master colortable. This is possible as only two values are given
# and so the flux only depends on the difference of the two
pfnametcolmag = 'BB_HSTcolorTempconverter.csv'
if os.path.exists(pfnametcolmag):
    print('Found and using color catalog {}'.format(pfnametcolmag))
    tcolmag = Table.read(pfnametcolmag, delimiter=',',
                         format='ascii.csv')
else:
    print('Making the color master')
    color = pmtable['m_F814W'] - pmtable['m_F606W']
    colstars = []
    coltemps = []
    colLp_mags = []
    colbins = np.linspace(np.min(color), np.max(color), 50)
    for icol, col in enumerate(colbins):
        print('Color {}. ({}/{})'.format(col, icol, len(colbins)))
        colstars.append(Star('color_'+str(col)))
        colstars[-1].mag = [[0., 0.01, 'HST_ACS_WFC.F606W_77'],
                            [col, 0.01, 'HST_ACS_WFC.F814W_77']]
        temparea = colstars[-1].get_temp_via_bb()
        colstars[-1].temperature = temparea[0] * u.K
        if temparea[0] == 1000.:
            import ipdb;ipdb.set_trace()
        colstars[-1].skyarea = temparea[1] * u.sr
        coltemps.append(temparea[0])
        colLp_mags.append(colstars[-1].temp2mag_from_filter('Paranal_NACO.Lp'))
    tcolmag = Table([colbins, coltemps*u.K, col - colLp_mags],
                    names=['Color', 'Temperature', '814Angs-Lp'])

    tcolmag.write(pfnametcolmag, delimiter=',',
                  format='ascii')

stars = []
temps = []
Lp_mags = []
nstars = len(pmtable)
for istar, m606, m814, m606err, m814err in zip(
        range(len(pmtable)), pmtable['m_F606W'], pmtable['m_F814W'],
        pmtable['rms_F606W(mag)'],
        pmtable['rms_F814W(mag)']):
    
    print('Star {}/{} (mags {}, {})'.format(istar, nstars, m606, m814))
    idx = find_nearest(tcolmag['Color'], m606-m814)
    temps.append(tcolmag['Temperature'][idx])
    Lp_mags.append(m814 - tcolmag['814Angs-Lp'][idx])
    print('Temp,Lp mag = {}, {}'.format(temps[-1], Lp_mags[-1]))

# this is when you want to get the flux for each individually
# for istar, m606, m814, m606err, m814err in zip(
#     range(len(pmtable)), pmtable['m_F606W'], pmtable['m_F814W'],
#         pmtable['rms_F606W(mag)'],
#         pmtable['rms_F814W(mag)']):
    
#     stars.append(Star(pmtable['ID'][istar]))
#     print('Star {}/{} (mags {}, {})'.format(istar, nstars, m606, m814))
#     stars[-1].mag = [[m606, m606err, 'HST_ACS_WFC.F606W_77'],
#                      [m814, m814err, 'HST_ACS_WFC.F814W_77']]
#     temparea = stars[-1].get_temp_via_bb()
#     stars[-1].temperature = temparea[0] * u.K
#     stars[-1].skyarea = temparea[1] * u.sr
#     temps.append(temparea[0])
#     Lp_mags.append(stars[-1].temp2mag_from_filter('Paranal_NACO.Lp'))
#     print('Temp,Lp mag = {}, {}'.format(temps[-1], Lp_mags[-1]))
    
pmtable['BB_Lpmag'] = Lp_mags
pmtable['BBtemps'] = temps
pmtable.sort('BB_Lpmag')
# decide what to save
outnames = ['MAIN_ID', 'RA', 'DEC', 'RA_err', 'DEC_err',
            'PMRA', 'PMRA_err', 'PMDEC', 'PMDEC_err',
            'N_used', 'MAG_EXTRAPOL',
            'm_F606W', 'm_F814W', 'm_F606W_err', 'm_F814W_err',
            'BB_Lpmag', 'BBtemps']
catnames = ['ID', 'RA_cor', 'DEC_cor', 'RA_err', 'DEC_err',
            'pmx', '1sig_pmx', 'pmy', '1sig_pmy',
            'N_used', 'MAG_EXTRAPOL',
            'm_F606W', 'm_F814W', 'rms_F606W(mag)', 'rms_F814W(mag)',
            'BB_Lpmag', 'BBtemps']


pmtable[catnames].write('47tuc_from_{}.csv'.format(tfn.split('.')[0]),
                        format='ascii', names=outnames,
                        delimiter=',', overwrite=True)


# Also save the results for vosa
vosalim = min(len(pmtable), vosalim)
dfilname = {'m_F606W': 'HST/ACS_WFC.F606W_77',
            'm_F814W': 'HST/ACS_WFC.F814W_77',
            'MAG_EXTRAPOL': 'Paranal/NACO.Lp'}
dusefit = {'m_F606W': '---',
           'm_F814W': '---',
           'MAG_EXTRAPOL': 'nofit'}
derror = {'m_F606W': 'rms_F606W(mag)',
          'm_F814W': 'rms_F814W(mag)',
          'MAG_EXTRAPOL': 'MAG_EXTRAPOL_err'}


print('Making the vosa output for the {} brightest stars (F814W <= {})'.format(
    vosalim, pmtable['m_F814W'][vosalim - 1]))
# splitting in multiple files for vosa
ifile = 0
while ifile * vosalim < len(pmtable):
    fvosa = open('bellini_vosa_{}.csv'.format(ifile), 'w')
    fvosa.write('#catalog from bellini for vosa (part ifile) \n')
    for itar, tarcol in enumerate(pmtable[ifile * vosalim:
                                          (ifile + 1) * vosalim]):
        if itar % 1000 == 0:
            print('Writing line {}/{}'.format(itar, len(pmtable)), end='\r')
        for filt in ['m_F606W', 'm_F814W', 'MAG_EXTRAPOL']:
            fvosa.write('{} {} {} {} {} {} {} {} {} {}\n'.format(
                tarcol['ID'], tarcol['RA_cor'], tarcol['DEC_cor'],
                '---', '0.0', dfilname[filt], tarcol[filt],
                tarcol[derror[filt]], dusefit[filt],
                '---'))
    print('Done writing part number {}'.format(ifile))

    fvosa.close()
    ifile += 1
print('Done writing all files')
ipdb.set_trace()
