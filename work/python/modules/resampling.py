import numpy as np
from scipy.signal.windows import hann
from tqdm import tqdm

def resampler_sinc_dyn(x, windowSize, SRO_ppm, outRange):

    y = np.zeros_like(x)
    x_len_min1 = np.size(x)-1
    ratio = 1+SRO_ppm/1e+6

    hh = hann(2*windowSize+1)

    for k in tqdm(range(outRange[0]-1, outRange[1])):

        k_ratio_rnd = np.round(ratio[k]*k)
        idx_start = int(max(0, k_ratio_rnd - windowSize))
        idx_end = int(min(x_len_min1, k_ratio_rnd + windowSize))
        m = np.arange(idx_start, idx_end+1)

        t1 = x[idx_start:idx_end+1] * np.sinc(k*ratio[k]-m)
        y[k] = np.sum( t1 * hh[-np.size(m):] )

    return y
