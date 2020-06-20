import sys, os, pickle
from glob import glob
sys.path.append('../lib')
import commands
from prettytable import PrettyTable
import json, yaml
import numpy as np
from scipy.interpolate import interp1d
from array import array

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

def getUncertaintyFromLimitTree(name, verbose=True):
    f = ur.open(name)
    r_arr = f['limit']['r'].array()
    nll_arr = f['limit']['deltaNLL'].array()
    c = r_arr[0]
    r_u = r_arr[r_arr>r_arr[0]]
    nll_u = nll_arr[r_arr>r_arr[0]]
    f_u = interp1d(nll_u, r_u, 'quadratic')
    u = f_u(0.5)
    r_l = r_arr[r_arr<r_arr[0]]
    nll_l = nll_arr[r_arr<r_arr[0]]
    f_l = interp1d(nll_l, r_l, 'quadratic')
    l = f_l(0.5)
    if verbose:
        print '----------------------------------'
        print 'R(D*) = {:.3f} +{:.3f}/-{:.3f} [{:.1f} %]'.format(c, u-c, c-l, 100*(u-l)*0.5/c)
        print 'Sigma = {:.3f}'.format((u-l)*0.5)
        print '----------------------------------\n'
    return c, c-l, u-c, (u-l)*0.5

def dumpDiffNuisances(output, outdir):
    name = []
    inVal = []
    outVal = []
    outDipls = []
    for line in output.split('\n')[3:]:
        aux = [i for i in line.split('  ') if i]
        if aux[0] == 'r': continue
        name.append(aux[0])
        inVal.append(aux[1])
        outVal.append(aux[3])

        xIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[0])
        xOut = float(aux[3].replace('!','').replace('*','').split(' +/- ')[0])
        sigIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[1])
        outDipls.append((xIn-xOut)/sigIn)

    outDipls = np.argsort(np.abs(np.array(outDipls)))
    for st_idx in [-15, 0]:
        idxList = outDipls[st_idx:]
        t = PrettyTable()
        t.field_names = ['Parameter', 'pre-fit', 'post-fit']
        for i in reversed(list(idxList)):
            t.add_row([name[i], inVal[i], outVal[i]])
        if st_idx == 0:
            with open(outdir+'/nuisance_difference.txt', 'w') as dumpfile:
                dumpfile.write('{}\n'.format(t))
        else: print t

stringJubCustomizationCaltechT2 = '''
+RunAsOwner = True
+InteractiveUser = True
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7"
+SingularityBindCVMFS = True
+MaxRuntime = 3600
RequestDisk = 200000
RequestMemory = 2500
RequestCpus = 1
x509userproxy = $ENV(X509_USER_PROXY)
'''
