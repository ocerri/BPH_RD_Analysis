from scipy.interpolate import interp1d
import pickle
import numpy as np

class pTCalReader:
    def __init__(self, calibration_file):
        self.kind = 'poly' if ('polyCoeff' in calibration_file) else 'ratio'

        if self.kind == 'poly':
            d = pickle.load(open( calibration_file, 'rb' ))
            self.beta = d['beta']
            self.betaVar = d['betaVar']
            self.nVar = len(d['betaVar'])

        elif self.kind == 'ratio':
            d = {}
            lines = open(calibration_file, 'r').readlines()
            keys = lines[0][1:-1].split('\t')
            for k in keys: d[k] = []

            for l in lines[1:]:
                l = l[:-1]
                v = l.split('\t')
                for i in range(len(v)):
                    d[keys[i]].append(float(v[i]))

            self.calibration_dic = d

            self.f = {}
            self.f['C'] = interp1d(d['pt'], d['w'],
                                   fill_value=(d['w'][0], d['w'][-1]),
                                   bounds_error=False,
                                   kind='cubic'
                                                    )
            self.f['Up'] = interp1d(d['pt'], d['wUp'],
                                   fill_value=(d['wUp'][0], d['wUp'][-1]),
                                   bounds_error=False,
                                   kind='cubic'
                                                    )
            self.f['Down'] = interp1d(d['pt'], d['wDown'],
                                   fill_value=(d['wDown'][0], d['wDown'][-1]),
                                   bounds_error=False,
                                   kind='cubic'
                                                    )

    def getWeights(self, B_pt, shape=0, scale=1., clipLimits=[0.7,1.7]):
        if self.kind == 'ratio':
            # Not implemented
            pass

        sign = np.sign(shape)
        idx = np.abs(shape) - 1
        delta_p = sign*scale*self.betaVar[idx]
        
        w = np.polyval(self.beta + delta_p, B_pt)
        return np.clip(w, clipLimits[0], clipLimits[1])
