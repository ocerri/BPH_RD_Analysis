import numpy as np
import scipy as sp

def doubleCrystalball(x, mu, sigma, beta, m):
    left = sp.stats.crystalball.pdf(x-mu, beta, m, loc=0, scale=sigma)
    right = sp.stats.crystalball.pdf(mu-x, beta, m, loc=0, scale=sigma)
    norm = 2*sp.stats.crystalball.cdf(0,beta,m,scale=sigma)
    return np.where(x < mu, left, right)/norm

def getBeamSpotWeights(ds, axis, parNew, parOld, ref):
    if ref == 'mean':
        x = 1e4*(ds['vtx_PV_'+axis] - np.mean(ds['vtx_PV_'+axis]))
    elif ref == 'bs':
        x = 1e4*(ds['vtx_PV_'+axis] - ds['beamSpot_'+axis])
    else:
        print 'Ref', ref, 'not accepted'
        raise
    return doubleCrystalball(x, *parNew)/doubleCrystalball(x, *parOld)

def getBeamSpotCorrectionWeights(ds, param, ref, dmu_x=0, dmu_y=0, clip_range=(0,10)):
    mu, sigma, beta, m = param['x']['data']
    mu += dmu_x
    out = getBeamSpotWeights(ds, 'x', [mu, sigma, beta, m], param['x']['MC'], ref)

    mu, sigma, beta, m = param['y']['data']
    mu += dmu_y
    out *= getBeamSpotWeights(ds, 'y', [mu, sigma, beta, m], param['y']['MC'], ref)

    return np.clip(out, clip_range[0], clip_range[1])
