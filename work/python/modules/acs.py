import numpy as np

class ACS:

    def __init__(self):
        self.preACS = 0
        self.ACS = 0
        self.mapping = 'EXP'

    def process(self, GCSD2_avg):
        # GCSD2_avg: exp. averaged GCSD2 up until nyquist-freq.
        self.preACS = np.mean(np.abs(GCSD2_avg), axis=0)
        self.ACS = self._map_EXP(self.preACS)
        return self.ACS

    def _map_EXP(self, preACS, binary=True):
        # exponential mapping in log-domain.
        a, b, c = -0.636, -3.362, 23.4
        basic_map = lambda x : a - b*np.exp(-c*x)
        chi_max = 0.25 # mapped to ACS=1
        maxLogSROEstErr = 0 # 1ppm
        out = basic_map(preACS) - maxLogSROEstErr
        out = out / np.abs(basic_map(chi_max))
        out = -out
        out[out < 0] = 0
        out[out > 1] = 1
        if binary:
            out[out > 0] = 1
        
        return out
