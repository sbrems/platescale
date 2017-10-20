# gives an estimate about the differential chromatic diffraction based on
# Sturms Diploma thiesis chapter 4, which is based on Stone (1996)
import astropy.units as u
from astropy.units import cds
from astropy.time import Time
import numpy as np

cds.enable()


# Input parameters/ header keywords
# temperature in degrees Celsius
hktemp = {'ESO': ['HIERARCH ESO TEL AMBI TEMP', u.deg_C]}
# air pressure in mmHg
hkpresss = {'ESO': ['HIERARCH ESO TEL AMBI PRES START', u.hPa]}
# watervapour pressure
hkpressw = {'ESO': ['calc_man', 1]}
# relative humidity in %
hkrel_humidity = {'ESO': 'HIERARCH ESO TEL AMBI RHUM'}
hktime = {'ESO': ['MJD-OBS', 'mjd']}
dsites = {'Paranal': [-24.6272, -70.4039] * u.deg}


def total_refract(header, wavel, site='Paranal', headtype='ESO'):
    '''Get refractional index based on Owens1976 using water and air pressure.
    Nomenclature follows Sturms thesis'''
    # read from the fitsheader
    temp = header[hktemp[headtype][0]] * hktemp[headtype][1]
    presss = header[hkpresss[headtype][0]] * hkpresss[headtype][1]
    if hkpressw[headtype][0] == 'calc_man':
        rel_humidity = header[hkrel_humidity[headtype]]
        pressw = temp2waterpressure(temp, rel_humidity)
    else:
        pressw = header[hkpressw[headtype][0]] * hkpressw[headtype][1]
    time = Time(header[hktime[headtype][0]], format=hktime[headtype][1])

    # get help numbers
    beta = 0.001245 * ((273.15 * u.deg_C + temp) / (273.15 * u.deg_C))
    T = 273.15 * u.deg_C + temp
    Ps = 1.333224 * (presss - pressw)
    Pw = 1.333224 * (pressw)
    sigma = 10**4 / wavel  # sign wrong in Sturm
    Ds = (1 + Ps.to('mmHg') / cds.mmHg * (57.90e-8 - 9.325e-4 * u.deg_C / T +
                                          0.25844 * u.deg_C**2 / T**2)) * Ps.to('mmHg') / T\
        / cds.mmHg * u.deg_C
    Di = -2.37321e-3 + 2.23366 * u. Celsius / T -\
        710.79207 * u.deg_C**2 / T**2 +\
        7.75141e4 * u.deg_C**3 / T**3
    Dw = (1 + Pw.to('mmHg') / cds.mmHg * (1 + 3.7e-4 * Pw.to('mmHg') / cds.mmHg) * Di) *\
        Pw.to('mmHg') / T / cds.mmHg * u.deg_C

    # finally get the refractional index n = gamma +1
    n0 = ((2371.34 + 683939.7 / (130 - sigma**2 * u.Angstrom**2) +
           4574.3 / (38.9 - sigma**2 * u.Angstrom**2)) * Ds +
          (6487.31 + 58.05 * sigma**2 * u.Angstrom**2 -
              0.71150 * sigma**4 * u.Angstrom**4 +
              0.08851 * sigma**6 * u.Angstrom**6) * Dw)\
        * 10**-8 + 1
    gamma = n0 - 1
    # die zenithdistanze
    ra = header['RA'] * u.deg
    dec = header['DEC'] * u.deg
    z0 = zenithdistance(time, ra, dec, site)

    # and the total refraction_index
    refract = gamma * (1 - beta) * np.tan(z0) -\
        gamma * (beta - gamma / 2) * np.tan(z0)**3
    import ipdb
    ipdb.set_trace()
    return refract


def temp2waterpressure(temp, rel_humidity):
    '''Based on the Magnus formula between -45 and 60 deg Celsius
     above water surfaces.'''
    sat_pressure = 6.112 * u.hPa * np.exp(17.62 * temp /
                                          (243.12 * u.deg_C + temp))
    return rel_humidity * sat_pressure


def zenithdistance(time, ra, dec, site):
    lat, lon = dsites[site]
    lst = time.sidereal_time('mean', lon)
    print('LST = {}'.format(lst))
    z0 = np.arccos(np.sin(np.deg2rad(dec)) * np.sin(np.deg2rad(lat)) +
                   np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(lat)) *
                   np.cos(np.deg2rad(lst - ra)))
    return z0
