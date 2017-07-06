from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy. units as u
import ipdb

tfn = 'NGC0104c.pm'

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
pmtable['MAG_USE'] = 15 * (pmtable['m_F814W'] - pmtable['m_F606W']) +\
    pmtable['m_F814W']

# decide what to save
outnames = ['MAIN_ID', 'RA', 'DEC', 'RA_err', 'DEC_err',
            'PMRA', 'PMRA_err', 'PMDEC', 'PMDEC_err',
            'N_used', 'MAG_USE',
            'm_F606W', 'm_F814W']
catnames = ['ID', 'RA_cor', 'DEC_cor', 'RA_err', 'DEC_err',
            'pmx', '1sig_pmx', 'pmy', '1sig_pmy',
            'N_used', 'MAG_USE',
            'm_F606W', 'm_F814W']
ipdb.set_trace()
pmtable[catnames].write('47tuc_from_{}.csv'.format(tfn.split('.')[0]),
                        format='ascii', names=outnames,
                        delimiter=',')
