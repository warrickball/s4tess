import numpy as np

TAU = 2*np.pi
CLIGHT = 299792.458  # km/s

def two_sines_plus_const(x, a0, a1, p1, a2, p2):
    return a0 + a1*np.sin(TAU*(x/360.0 + p1)) + a2*np.sin(TAU*(x/180.0 + p2))


def mean_v(lon):
    popt = [1.37622546, 18.8420428, 0.577075134, 7.93274147, -0.0148034692]
    return two_sines_plus_const(lon, *popt)


def sigma_v(lon):
    popt = [30.01595341, -1.1142313, 0.72453242, 4.34319137, 0.19605259]
    return two_sines_plus_const(lon, *popt)


def usno_lon_sun(JD):
    # http://aa.usno.navy.mil/faq/docs/SunApprox.php
    # as used in TESS GI software e.g.
    # https://github.com/tessgi/tobs/blob/master/tobs/tobs.py#L60
    D = JD - 2451545.0
    q = 280.459 + 0.98564736*D
    g = np.radians(357.529 + 0.98560028*D)
    return np.mod(q + 1.915*np.sin(g) + 0.020*np.sin(2.*g), 360.0)

def sin_plus_const(t, const, amp, freq, phase):
    return const + amp*np.sin(2.*np.pi*(freq*t + phase))

def sector_duration(i):
    # from best fit of _sin_plus_const to data from Tom Barclay
    # https://figshare.com/articles/TESS_Perigee_Times/6875525
    popt = np.array([2.72758623e+01, 1.49295494e+00, 3.45331787e-03, 1.70842155e-01])
    popt[:2] = popt[:2]*720
    return int(sin_plus_const(i/720, *popt))

def sin_plus_line(n, c0, c1, a1, f1, p1):
    return c0 + c1*n + a1*np.sin(2.*np.pi*(f1*n + p1))

sector_starts = [0]  # in cadences
sector_lengths = []  # in cadences

for i in range(26):
    sector_lengths.append(sector_duration(sector_starts[i]))
    sector_starts.append(sector_starts[i] + sector_lengths[i])

# assume pointing advances by (360 degrees)*(sector length/year)
# sector_lons = 315.8 + np.array(sector_starts)/720*360/365.5
# sector_lons = np.mod(sector_lons, 360.0)

# use fit of sine + line to Barclay data to estimate pointings
start_JD = 2458324.701388889
popt = np.array([15.22068518, 27.26620879,  2.38632187,  0.09370191, -0.06827264])
sector_lons = np.mod(usno_lon_sun(sin_plus_line(np.arange(26.0), *popt) + start_JD) + 180.0, 360.0)


def sigma_stassun(Imag):
    # see Stassun et al. (2017)
    # https://arxiv.org/abs/1706.00495v3

    # if Imag is less than local minimum aroung 5,
    # use flat line that roughly corresponds to systematic noise
    Imag = np.maximum(Imag, 4.97756051)
    F = 4.73508403525e-5
    E = -0.0022308015894
    D = 0.0395908321369
    C = -0.285041632435
    B = 0.850021465753
    lnA = 3.29685004771
    lnsigma = lnA + B*Imag + C*Imag**2 + D*Imag**3 + E*Imag**4 + F*Imag**5
    return np.exp(lnsigma)

def sigma_schofield(Imag, Teff, ELon=0.0, ELat=30.0, GLon=96.0, GLat=-30.0,
                    cadence=120.0, subcadence=2.0, npix_aper=4, frac_aper=0.76,
                    e_pix_ro=10, geom_area=60.0, pix_scale=21.1, sys_limit=0.0):
    # from Mat Schofield's ATL code
    # "from the TESS Wiki"...

    omega_pix = pix_scale**2
    n_exposures = cadence/subcadence

    # electrons from the star
    megaph_s_cm2_0mag = 1.6301336 + 0.14733937*(Teff-5000.0)/5000.0
    e_star = 10.0**(-0.4*Imag)*1e6*megaph_s_cm2_0mag*geom_area*cadence*frac_aper
    e_star_sub = e_star*subcadence/cadence

    # e/pix from zodi
    dlat = (abs(ELat) - 90.0)/90.0
    vmag_zodi = 23.345 - 1.148*dlat**2
    e_pix_zodi = 10.0**(-0.4*(vmag_zodi-22.8))*2.39e-3*geom_area*omega_pix*cadence

    # e/pix from background stars
    dlat = abs(GLat)/40.0

    dlon = GLon
    # q = np.where(dlon > 180.0)
    # if len(q[0]) > 0:
    # 	dlon[q] = 360.0 - dlon[q]
    # if dlon > 180.0:
    #     dlon = 360.0 - dlon
    dlon = np.minimum(dlon, 360.0 - dlon)

    dlon = abs(dlon)/180.0
    p = [18.97338, 8.833, 4.007, 0.805]
    Imag_bgstars = p[0] + p[1]*dlat + p[2]*dlon**p[3]
    e_pix_bgstars = 10.0**(-0.4*Imag_bgstars)*1.7e6*geom_area*omega_pix*cadence

    # compute noise sources
    noise_star = np.sqrt(e_star)/e_star
    noise_sky  = np.sqrt(npix_aper*(e_pix_zodi + e_pix_bgstars))/e_star
    noise_ro   = np.sqrt(npix_aper*n_exposures)*e_pix_ro/e_star
    noise_sys  = 0.0*noise_star + sys_limit/1e6/np.sqrt(cadence/3600.0)

    noise = np.sqrt(noise_star**2 + noise_sky**2 + noise_ro**2 + noise_sys**2)
    return noise*1e6*np.sqrt(cadence/3600.0)  # parts x sqrt(hr) -> ppm * sqrt(cadence)

def tess_sectors(lon, lat):
    """Returns a NumPy boolean array of which sectors the given stars
    appear.

    Parameters
    ----------
    lon: float or array
        Target stars' ecliptic longitudes.
    lat: float or array
        Target stars' ecliptic latitudes.

    Returns
    -------
    sectors: array
        Boolean NumPy array whose ``(i,j)``th element is ``True`` if
        the ``i``th star appears in the ``j``th sector.

    """
    Nsectors = 26

    one = False
    try:
        sectors = np.full((len(lon), Nsectors), False)
    except TypeError:
        one = True
        sectors = np.full((1, Nsectors), False)
        
    lon_range = np.minimum(12.0/abs(np.cos(np.radians(lat))), 180.0)

    for i in range(Nsectors):
        a = sector_lons[i]
        d1 = abs(lon - a)
        d1 = np.minimum(d1, 360.0 - d1)
        d2 = abs(lon - (a + 180.0)%360.0)
        d2 = np.minimum(d2, 360.0 - d2)
        h = (-1)**(i//13 + 1)  # sign for latitude of hemisphere
        sectors[:,i] = ((d1 <= lon_range) & (6.0 <= h*lat)) | ((d2 <= lon_range) & (78.0 <= h*lat))

    if one:
        return sectors[0]
    else:
        return sectors

    
def load_txt(filename):
    """Return data from `filename` with my simple text format (each line
    is `key = value`) in a dictionary.
    """
    d = {}
    with open(filename, 'r') as f:
        for line in f.readlines():
            k, v = line.split('=')
            d[k.strip()] = float(v)

    return d


def save_txt(filename, d):
    """Save a dictionary `d` to `filename` with my simple text format
    (each line is `key = value`).
    """
    with open(filename, 'w') as f:
        if type(d) is dict:
            for k, v in d.items():
                f.write('{:>16s} = {:.16g}\n'.format(k, v))
        elif type(d) is np.ndarray:
            for k in d.dtype.names:
                f.write('{:>16s} = {:.16g}\n'.format(k, float(d[k])))


def compare_sequence(a, b):
    for ai, bi in zip(a, b):
        assert(ai == bi)

def expand_bools(s):
    return [True if si == 'T' else False for si in s]
               

if __name__ == "__main__":
    print('Testing tess_sectors...')
    # test poles and ecliptic
    compare_sequence(tess_sectors(0.0, 90.0),
                     expand_bools('FFFFFFFFFFFFFTTTTTTTTTTTTT'))
    compare_sequence(tess_sectors(0.0, -90.0),
                     expand_bools('TTTTTTTTTTTTTFFFFFFFFFFFFF'))
    compare_sequence(tess_sectors(315.0, 0.0),
                     expand_bools('FFFFFFFFFFFFFFFFFFFFFFFFFF'))

    # test multiple latitudes
    truths = tess_sectors(sector_lons[-1]*np.ones(5),
                          np.array([12., 24., 36., 48., 60.]))
    for i, truth in enumerate(truths):
        compare_sequence(truth, expand_bools('FFFFFFFFFFFFFFFFFFFFFFFFFT'))

    # test multiple longitudes
    truths = tess_sectors(sector_lons[:13], -12.0*np.ones(13))
    for i, truth in enumerate(truths):
        reference = [False]*26
        reference[i] = True
        compare_sequence(truth, reference)

    print('Sector longitudes are:')
    for i, lon in enumerate(sector_lons):
        print('%2i %7.2f' % (i, lon))
        
    print('Done.')
