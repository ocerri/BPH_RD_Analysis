import sys, os, pickle, copy
from glob import glob
sys.path.append('../lib')
import commands
from prettytable import PrettyTable
import json, yaml
import numpy as np
from scipy.interpolate import interp1d
from array import array
import matplotlib.pyplot as plt

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

def loadHisto4CombineFromRoot(histo_file_dir, card_name, loadShapeVar=False):
    histo = {}
    for fname in glob(histo_file_dir+'{}_[EAMU]*.root'.format(card_name)):
        regionName = os.path.basename(fname)[len(card_name)+1:-5]
        tfReg = rt.TFile.Open(histo_file_dir+'{}_{}.root'.format(card_name, regionName), 'READ')
        histo[regionName] = {}
        for n in [k.GetName() for k in tfReg.GetListOfKeys()]:
            if not loadShapeVar and '__' in n:
                continue
            key = 'data' if 'data' in n else n
            histo[regionName][key] = copy.deepcopy(tfReg.Get(n))
    return histo

def getUncertaintyFromLimitTree(name, verbose=True, drawPlot=False):
    f = ur.open(name)
    iToy_arr = f['limit']['iToy'].array()
    r_arr_raw = f['limit']['r'].array()
    nll_arr_raw = f['limit']['deltaNLL'].array()
    res = []
    for iToy in np.unique(iToy_arr):
        sel = iToy_arr == iToy
        r_arr = r_arr_raw[sel]
        nll_arr = nll_arr_raw[sel]
        c = r_arr[0]

        r_u = r_arr[r_arr>c]
        nll_u = nll_arr[r_arr>c]

        r_l = r_arr[r_arr<c]
        nll_l = nll_arr[r_arr<c]

        if drawPlot:
            color = ['b', 'g', 'r', 'c', 'm', 'y', 'k'][iToy]
            plt.plot(2*nll_u, r_u, '.--', color=color)
            plt.plot(2*nll_l, r_l, '.--', color=color)
            plt.plot(nll_arr[0], c, 'o', color='gray')
            idxL = np.argmax(nll_l < 4) if np.max(nll_l) >= 4 else -1
            idxU = np.argmax(nll_u < 4) if np.max(nll_u) >= 4 else -1
            plt.ylim(r_l[idxL], r_u[idxU])
            plt.xlabel('$-2\Delta\log(L)$')
            plt.ylabel('POI')
            plt.grid()
        if np.all(nll_l[:-1] >= nll_l[1:]) and np.all(nll_u[:-1] <= nll_u[1:]):
            f_l = interp1d(nll_l, r_l, 'quadratic', fill_value='extrapolate')
            l = f_l(0.5)
            l2 = f_l(2)
            f_u = interp1d(nll_u, r_u, 'quadratic', fill_value='extrapolate')
            u = f_u(0.5)
            u2 = f_u(2.0)
        else:
            if not np.all(nll_l[:-1] >= nll_l[1:]):
                print 'Low error'
                print nll_l
            if not np.all(nll_u[:-1] <= nll_u[1:]):
                print 'High error'
                print nll_u
            print 'ERROR: X array not sorted'
            u, l = 0, 0
            u2, l2 = 0, 0
        if verbose:
            print '----------------------------------'
            if iToy:
                print 'Toy', iToy
            print 'R(D*) = {:.3f} +{:.3f}/-{:.3f} [{:.1f} %]'.format(c, u-c, c-l, 100*(u-l)*0.5/c)
            print 'Sigma = {:.3f}'.format((u-l)*0.5)
        res.append([c, c-l, u-c, (u-l)*0.5, l2, l, u, u2])
    if verbose:
        print '----------------------------------\n'
    return np.array(res) if len(res) > 1 else res[0]

def dumpDiffNuisances(output, outdir, tag='', useBonlyResults=False, parsToPrint=15):
    name = []
    inVal = []
    outVal = []
    outDipls = []
    for line in output.split('\n')[3:]:
        aux = [i for i in line.split('  ') if i]
        if aux[0] == 'r': continue
        name.append(aux[0])
        inVal.append(aux[1])
        idxOutVal = 2 if useBonlyResults else 3
        outVal.append(aux[idxOutVal])

        xIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[0])
        xOut = float(aux[idxOutVal].replace('!','').replace('*','').split(' +/- ')[0])
        sigIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[1])
        outDipls.append((xIn-xOut)/sigIn)

    outDipls = np.argsort(np.abs(np.array(outDipls)))
    for st_idx in [-parsToPrint, 0]:
        idxList = outDipls[st_idx:]
        t = PrettyTable()
        t.field_names = ['Parameter', 'pre-fit', 'post-fit']
        for i in reversed(list(idxList)):
            t.add_row([name[i], inVal[i], outVal[i]])
        if st_idx == 0:
            with open(outdir+'/nuisance_difference'+('_'+tag if tag else '')+'.txt', 'w') as dumpfile:
                dumpfile.write('{}\n'.format(t))

            print 'Dumping text table'
            with open(outdir+'/nuisance_difference'+('_'+tag if tag else '')+'_texTable.txt', 'w') as dumpfile:
                cols = 2
                for i in range(1, idxList.shape[0]):
                    if not i%cols == 1:
                        continue

                    s = ''
                    for j in range(cols):
                        idx = idxList[-i-j]
                        s += name[idx].replace('prop_bin', 'stat_').replace('_AddTk', '').replace('_', '\_')
                        s += ' & $'
                        s += outVal[idx].replace('!', '').split(' (')[0].replace('+/-', '\pm') + '$'
                        if j == cols -1 :
                            s += ' \\\\'
                        else:
                            s += ' & '
                    dumpfile.write(s + '\n')
                    lastVal = outVal[idx].replace('!', '').split(' +/-')[0]
                    if np.abs(float(lastVal)) < 1.6:
                        break

        else:
            print t

stringJubCustomizationCaltechT2 = '''
+RunAsOwner = True
+InteractiveUser = True
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m20200605"
+SingularityBindCVMFS = True
+MaxRuntime = 3600
RequestDisk = 200000
RequestMemory = 2500
RequestCpus = 1
x509userproxy = $ENV(X509_USER_PROXY)
'''
