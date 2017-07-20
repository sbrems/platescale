from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy. units as u
import numpy as np
import ipdb

tfn = 'NGC0104c.pm'
vosalim = 9000  # only take the xx brightest stars. Vosa can handle only
# ~10000 at  a time
pmtable = Table.read(tfn, format='ascii.csv',
                     delimiter='\s', comment='#')
# apply suggested corrections
zerocoord = SkyCoord(['00:24:05.71 -72:04:52.70'],
                     unit=[u.hourangle, u.deg], frame='icrs')
pmtable['RA_cor'] = ((pmtable['x_M'] - pmtable['Delta_x']) * u.mas +
                     zerocoord.ra).to('deg')
pmtable['DEC_cor'] = ((pmtable['y_M'] - pmtable['Delta_y']) * u.mas +
                      zerocoord.dec).to('deg')
pmtable['RA_err'] = [0, ] * len(pmtable) * u.deg
pmtable['DEC_err'] = [0, ] * len(pmtable) * u.deg
pmtable['MAG_EXTRAPOL'] = 15 * (pmtable['m_F814W'] - pmtable['m_F606W']) +\
    pmtable['m_F814W']
pmtable['MAG_EXTRAPOL_err'] = np.sqrt(pmtable['rms_F606W(mag)']**2 +
                                      pmtable['rms_F814W(mag)']**2)

print('Fitting BB to all the stars')

stars = []
temps = []
Lp_mags = []
nstars = len(pmtable)
for istar, m606, m814, m606err, m814err in enumerate(
        pmtable['m_F606W', 'm_F814W', 'rms_F606W', 'rms_F814W']):
    stars.append(Star(pmtable['ID']))
    print('Star {}/{}'.format(istar,nstars))
    stars[-1].mag = [[m606, m606err, 'HST_ACS_WFC.F606W_77'],
                     [m814, m814err, 'HST_ACS_WFC.F814W_77']]
    temparea = stars[-1].get_temp_via_bb()
    stars[-1].temperature = temparea[0]
    stars[-1].skyarea = temparea[1]
    temps.append(temparea[0])
    Lp_mags.append(stars[-1].temp2mag_from_filter('Paranal_NACO.Lp'))
    



# decide what to save
outnames = ['MAIN_ID', 'RA', 'DEC', 'RA_err', 'DEC_err',
            'PMRA', 'PMRA_err', 'PMDEC', 'PMDEC_err',
            'N_used', 'MAG_EXTRAPOL',
            'm_F606W', 'm_F814W', 'm_F606W_err', 'm_F814W_err']
catnames = ['ID', 'RA_cor', 'DEC_cor', 'RA_err', 'DEC_err',
            'pmx', '1sig_pmx', 'pmy', '1sig_pmy',
            'N_used', 'MAG_EXTRAPOL',
            'm_F606W', 'm_F814W', 'rms_F606W(mag)', 'rms_F814W(mag)']



pmtable[catnames].write('47tuc_from_{}.csv'.format(tfn.split('.')[0]),
                        format='ascii', names=outnames,
                        delimiter=',', overwrite=True)



# Also save the results for vosa
pmtable.sort('m_F814W')
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
                '---', '---', dfilname[filt], tarcol[filt],
                tarcol[derror[filt]], dusefit[filt],
                'Av:0.0/3.0'))
    print('Done writing part number {}'.format(ifile))

    fvosa.close()
    ifile += 1
print('Done writing all files')
###############################################ü
