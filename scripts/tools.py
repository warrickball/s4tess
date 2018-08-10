import numpy as np

def _sin_plus_const(t, const, amp, freq, phase):
    return const + amp*np.sin(2.*np.pi*(freq*t + phase))

def _sector_duration(i):
    # from best fit of _sin_plus_const to data from Tom Barclay
    # https://figshare.com/articles/TESS_Perigee_Times/6875525
    popt = np.array([2.72758623e+01, 1.49295494e+00, 3.45331787e-03, 1.70842155e-01])
    popt[:2] = popt[:2]*720
    return int(_sin_plus_const(i/720, *popt))
    
sector_starts = [0]  # in cadences
sector_lengths = []  # in cadences

for i in range(39):
    sector_lengths.append(_sector_duration(sector_starts[i]))
    sector_starts.append(sector_starts[i] + sector_lengths[i])
    

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
            for k in d.type.names:
                f.write('{:>16s} = {:.16g}\n'.format(k, d[k]))
