import numpy as np
import uproot as ur
import ROOT as rt
from glob import glob

import operator
ops = {'>': operator.gt, '<': operator.lt, }

def drawOnCMSCanvas(CMS_lumi, dobj, opt = None, tag=''):
    c = rt.TCanvas('c'+tag, 'c'+tag, 50, 50, 800, 600)
    c.SetTickx(0)
    c.SetTicky(0)

    if dobj.__class__ == rt.RooPlot:
        dobj.Draw()
    elif dobj[0].__class__ in [rt.TH1F, rt.TH1D, rt.TH2D]:
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
        print 'Loading branch', k
        m = []
        for i, e in enumerate(t.array(k)):
            m += list(e)
        l[k] = np.array(m)

    return l

def extarct_multiple(fname, branches = [], flag=''):
    if len(branches) == 0:
        print 'Must give a branches list'
    l = {}
    for b in branches:
        l[b] = []

    flist = glob(fname)

    for i,f in enumerate(flist):
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

    for b in branches:
        l[b] = np.array(l[b])
    return l

def createSel(d, cut):
    sel = np.ones_like(d[cut.keys()[0]], dtype=bool)
    for k, v in cut.iteritems():
        sel = np.logical_and(sel, ops[v[0]](d[k], v[1]) )
    return sel
