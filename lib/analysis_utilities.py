import numpy as np
import uproot as ur
import ROOT as rt
from glob import glob
import yaml
import os

import operator
ops = {'>': operator.gt, '<': operator.lt, }

def drawOnCMSCanvas(CMS_lumi, dobj, opt = None, tag='', size=[800,600]):
    c = rt.TCanvas('c'+tag, 'c'+tag, 50, 50, size[0], size[1])
    c.SetTickx(0)
    c.SetTicky(0)

    if dobj.__class__ == rt.RooPlot:
        dobj.Draw()
    elif dobj[0].__class__ in [rt.TH1F, rt.TH1D, rt.TH2D, rt.TGraph, rt.TGraphErrors, rt.TGraphAsymmErrors, rt.TProfile, rt.TEfficiency]:
        for i, o in enumerate(dobj):
            do = ''
            if not (opt is None):
                if opt == 'same':
                    if i>0:
                        do = 'SAME'
                else:
                    do = opt[i]
            o.Draw(do)
    else:
        print dobj.__class__
        print dobj[0].__class__
        print 'Class not recognized'
        raise


    CMS_lumi.CMS_lumi(c, -1, 0)
    c.obj = dobj
    c.Draw()
    return c

def extarct(t, branches = []):
    if len(branches) == 0:
        branches = t.keys()
    l = {}
    for k in branches:
#         print 'Loading branch', k
        m = []
        for i, e in enumerate(t.array(k)):
            if e.__class__ == np.float32:
                m += [e]
            else:
                m += list(e)
        l[k] = np.array(m)

    return l

def extarct_multiple(fname, branches = [], flag=''):
    if len(branches) == 0:
        print 'Must give a branches list'
    l = {}
    for b in branches:
        l[b] = []

    if not isinstance(fname, list):
        flist = glob(fname)
    else:
        flist = fname

    for i,f in enumerate(flist):
        try:
            t = ur.open(f)
            if 'outA;1' in t.keys():
                t=t['outA']['Tevts']
                for k in branches:
                    if flag=='data' and k[:2] == 'MC':
                        continue
                    if not (k in t.keys()):
                        continue
                    for i, e in enumerate(t.array(k)):
                        try:
                            l[k] += list(e)
                        except:
                            l[k] += [e]
        except:
            print 'Error in file:', f

    for b in branches:
        l[b] = np.array(l[b])
    return l

def createSel(d, cut):
    sel = np.ones_like(d[cut.keys()[0]], dtype=bool)
    for k, v in cut.iteritems():
        sel = np.logical_and(sel, ops[v[0]](d[k], v[1]) )
    return sel

def getEff(k,N):
    e = k/float(N)
    de = np.sqrt(e*(1-e)/N)
    return [e, de]

class DSetLoader(object):
    def __init__(self, in_sample,
                 candLoc='/storage/user/ocerri/BPhysics/data/cmsMC_private/',
                 candDir='ntuples_B2DstMu',
                 site_loc_conf = '/mnt/hadoop/store/user/ocerri',
                 sampleFile = '/storage/user/ocerri/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/production/samples.yml'
                 ):
        samples = yaml.load(open(sampleFile))['samples']
        if not in_sample in samples.keys():
            raise
        self.sample = in_sample
        self.candLoc = candLoc
        self.candDir = candDir

        self.MINIAOD_dirs = []
        for part in samples[self.sample]['parts']:
            aux = glob(part)
            if len(aux) > 0:
                aux = [os.path.dirname(part)]
            else:
                aux = glob(site_loc_conf + part[:-38].replace('ocerri-','') + '/*/*')
            self.MINIAOD_dirs += aux

        self.full_name = samples[self.sample]['dataset']

        res = glob(os.path.join(self.candLoc, self.full_name, self.candDir))
        if res:
            self.ntuples_dir = res[0]
            self.skimmed_dir = os.path.join(self.ntuples_dir, 'skimmed')
        else:
            self.ntuples_dir = ''
            self.skimmed_dir = ''

        effMCgenFile = os.path.join(self.candLoc, self.full_name, 'effMCgenerator.yaml')
        if os.path.isfile(effMCgenFile):
            self.effMCgen = yaml.load(open(effMCgenFile, 'r'))

        effCandFile = os.path.join(self.ntuples_dir, 'effCAND.yaml')
        if os.path.isfile(effCandFile):
            self.effCand = yaml.load(open(effCandFile, 'r'))

    def getSkimEff(self, catName='probe'):
        if catName is 'probe':
            with open(self.skimmed_dir + '/selTree.log', 'r') as f:
                aux = f.readlines()[-1][:-1].split(' ')
                return [float(aux[1])*1e-2, float(aux[3])*1e-2]
        else:
            with open(self.skimmed_dir + '/{}.log'.format(catName), 'r') as f:
                aux = f.readlines()[-1][:-1].split(' ')
                return [float(aux[1])*1e-2, float(aux[3])*1e-2]
