import numpy as np

def load_txt(filename):
    d = {}
    with open(filename, 'r') as f:
        for line in f.readlines():
            k, v = line.split('=')
            d[k.strip()] = float(v)

    return d

def save_txt(filename, d):
    with open(filename, 'w') as f:
        if type(d) is dict:
            for k, v in d.items():
                f.write('{:>16s} = {:.16g}\n'.format(k, v))
        elif type(d) is np.ndarray:
            for k in d.type.names:
                f.write('{:>16s} = {:.16g}\n'.format(k, d[k]))
        
