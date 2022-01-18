import sys, os, pickle, copy, re
from glob import glob
sys.path.append('../lib')
import commands
from prettytable import PrettyTable
import json, yaml
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from array import array
import matplotlib.pyplot as plt

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

def loadHisto4CombineFromRoot(histo_file_dir, card_name, loadShapeVar=False, verbose=False):
    if not histo_file_dir[-1] == '/':
        histo_file_dir += '/'
    histo = {}
    for fname in glob(histo_file_dir+'{}_*.root'.format(card_name)):
        regionName = os.path.basename(fname)[len(card_name)+1:-5]
        if regionName.startswith('NoMCstats_'): continue
        if regionName.startswith('blinded_'): continue
        if regionName.startswith('Asimov_'): continue
        if regionName.startswith('_'): continue
        if verbose:
            print regionName
            print 'Loading histos from:', fname
        tfReg = rt.TFile.Open(fname, 'READ')
        histo[regionName] = {}
        for n in [k.GetName() for k in tfReg.GetListOfKeys()]:
            if not loadShapeVar and '__' in n:
                continue
            key = 'data' if n == 'data_obs' else n
            histo[regionName][key] = copy.deepcopy(tfReg.Get(n))
    return histo

def getUncertaintyFromLimitTree(name, verbose=True, drawPlot=False):
    def y_lineByPoints(x, x1, x2, y1, y2):
        # (x - x1) / (x2 -x1) = (y - y1)/(y2 - y1)
        if x1 != x2:
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        else:
            return (y1 + y2)*0.5

    if not os.path.isfile(name):
        print 'File not found:', name
        raise
    arr = rtnp.root2array(name, treename='limit')
    iToy_arr = arr['iToy']
    r_arr_raw = arr['r']
    nll_arr_raw = arr['deltaNLL']
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
            idxU = np.argmax(nll_u < 4) if np.max(nll_u) >= 4 else 0
            plt.ylim(r_l[idxL], r_u[idxU])
            plt.xlabel('$-2\Delta\log(L)$')
            plt.ylabel('POI')
            plt.grid()

        if np.all(nll_l[:-1] >= nll_l[1:]):
            if len(r_l) > 2:
                kind = 'linear' if 1 > np.max(nll_l) else 'quadratic'
                f_l = interp1d(nll_l, r_l, kind, fill_value='extrapolate')
                l = f_l(0.5)
                l2 = f_l(2)
            else:
                print '[WARNING] Minimum at lower bound'
                if len(r_l):
                    l, l2 = r_l[0], r_l[0]
                else:
                    l, l2 = c, c
        else:
            print 'Low error'
            print nll_l
            print '[ERROR] Likelihood not convex (X array not sorted)'
            l, l2 = 0, 0
            raise


        if np.all(nll_u[:-1] <= nll_u[1:]):
            if len(r_u) > 2:
                kind = 'linear' if 1 > np.max(nll_u) else 'quadratic'
                f_u = interp1d(nll_u, r_u, kind, fill_value='extrapolate')
                u = f_u(0.5)
                u2 = f_u(2.0)
            else:
                print '[WARNING] Minimum at upper bound'
                if len(r_u):
                    u, u2 = r_u[-1], r_u[-1]
                else:
                    u, u2 = c, c
        else:
            print 'High error'
            print nll_u
            print '[ERROR] Likelihood not convex (X array not sorted)'
            i = 0
            u, u2 = None, None
            while i < len(nll_u):
                if nll_u[i] < 0.5 or (u is not None and nll_u[i] < 2.):
                    i += 1
                    continue
                elif u is None:
                    aux_r_u = [c] + list(r_u)
                    aux_nll_u = [0] + list(nll_u)

                    x1, x2 = aux_nll_u[i], aux_nll_u[i+1]
                    y1, y2 = aux_r_u[i], aux_r_u[i+1]
                    u = y_lineByPoints(0.5, x1, x2, y1, y2)
                elif u2 is None:
                    x1, x2 = nll_u[i-1], nll_u[i]
                    y1, y2 = r_u[i-1], r_u[i]
                    u2 = y_lineByPoints(2., x1, x2, y1, y2)
                else:
                    break
            if u is None:
                imax = np.argmax(nll_u)
                x1, x2 = nll_u[imax], 0
                y1, y2 = r_u[imax], c
                u = y_lineByPoints(0.2, x1, x2, y1, y2)
            if u2 is None:
                u2 = 2*u


        if verbose:
            print '----------------------------------'
            if iToy:
                print 'Toy', iToy
            outStr = 'R(D*) = {:.3f} +{:.3f}/-{:.3f}'.format(c, u-c, c-l)
            if c > 9e-4:
                outStr += ' [{:.1f} %]'.format(100*(u-l)*0.5/c)
                outStr += '\nSigma = {:.3f}'.format((u-l)*0.5)
            print outStr
        res.append([c, c-l, u-c, (u-l)*0.5, l2, l, u, u2])
    if verbose:
        print '----------------------------------\n'
    return np.array(res) if len(res) > 1 else res[0]


def getResultsFromMultiDimFitSingles(name, verbose=True, getTrackedParam=False, getMCstatNuisance=False):
    if not os.path.isfile(name):
        print 'File not found:', name
        raise
    arr = pd.DataFrame(rtnp.root2array(name, treename='limit'))

    res = []
    trackedParam = {}
    for varName in arr.columns:
        if varName.startswith('trackedParam_'):
            n = varName[13:]
            if n.startswith('prop_bin') and not getMCstatNuisance:
                continue
            trackedParam[n] = []
    for iToy in np.unique(arr['iToy']):
        sel = arr['iToy'] == iToy

        c, l, u = arr['r'][sel]

        for n in trackedParam.keys():
            trackedParam[n].append(arr['trackedParam_'+n][0])

        if verbose:
            print '----------------------------------'
            if iToy:
                print 'Toy', iToy
            print 'R(D*) = {:.3f} +{:.3f}/-{:.3f} [{:.1f} %]'.format(c, u-c, c-l, 100*(u-l)*0.5/c)
            print 'Sigma = {:.3f}'.format((u-l)*0.5)
        res.append([c, (u-l)*0.5, c-l, u-c, l, u])
    if verbose:
        print '----------------------------------\n'

    res = np.array(res) if len(res) > 1 else res[0]
    if getTrackedParam:
        return res, trackedParam
    else:
        return res

def dumpDiffNuisances(output, outdir, tag='', useBonlyResults=False, parsToPrint=15):
    name = []
    inVal = []
    outVal = []
    outDipls = []
    for line in output.split('\n')[3:]:
        aux = [i for i in line.split('  ') if i]
        if aux[0] == 'r': continue
        name.append(aux[0])

        if ' +/- ' in aux[1]:
            inVal.append(aux[1])
            xIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[0])
            sigIn = float(aux[1].replace('!','').replace('*','').split(' +/- ')[1])
        elif '[' in aux[1][:3] and ']' == aux[1][-1]: #is a rate param
            name[-1] += ' (rate param.)'
            aux[1] = aux[1].replace('[', '').replace(']', '').split(',')
            inVal.append('[{:+1.2e}, {:+1.2e}]'.format(float(aux[1][0]), float(aux[1][1])))
            xIn = 0.5*(float(aux[1][1]) + float(aux[1][0]))
            sigIn = 1e12

        idxOutVal = 2 if useBonlyResults else 3
        outVal.append(aux[idxOutVal])
        xOut = float(aux[idxOutVal].replace('!','').replace('*','').split(' +/- ')[0])
        outDipls.append((xIn-xOut)/sigIn)

    outputSigmas = outDipls
    outDipls = np.argsort(np.abs(np.array(outDipls)))

    # Print top parameters
    t = PrettyTable()
    t.field_names = ['Parameter', 'pre-fit', 'post-fit']
    i = outDipls.shape[0] - 1
    nPrinted = 0
    while i>=0 and nPrinted < parsToPrint:
        idx = outDipls[i]
        if name[idx].startswith('prop_bin'):
            i -= 1
        else:
            t.add_row([name[idx], inVal[idx], outVal[idx]])
            nPrinted += 1
            i -= 1
    print t

    # Saving the whole diff
    t = PrettyTable()
    t.field_names = ['Parameter', 'pre-fit', 'post-fit']
    for i in reversed(list(outDipls)): t.add_row([name[i], inVal[i], outVal[i]])

    with open(outdir+'/nuisance_difference'+('_'+tag if tag else '')+'.txt', 'w') as dumpfile:
        dumpfile.write('{}\n'.format(t))

    with open(outdir+'/nuisance_difference'+('_'+tag if tag else '')+'_texTable.txt', 'w') as dumpfile:
        cols = 2
        # s = '\\hline\n'
        # s += r' & '.join(cols*[r'Parameter & Postfit [$\sigma$]']) + r'\\' + '\n'
        # s += '\\hline\n'
        # dumpfile.write(s)
        for i in range(1, outDipls.shape[0]):
            if not i%cols == 1:
                continue

            s = ''
            for j in range(cols):
                idx = outDipls[-i-j]
                bName = name[idx].replace('prop_bin', '').replace('_AddTk', '')
                if 'bin' in bName:
                    bName = bName.replace('_bin', ' (') + ')'
                bName = bName.replace('_', ' ')
                s += bName
                s += ' & $'
                s += outVal[idx].replace('*', '').replace('!', '').split(' (')[0].replace('+/-', '\pm') + '$'
                if j == cols -1 :
                    s += ' \\\\'
                else:
                    s += ' & '
            dumpfile.write(s + '\n')
            lastVal = outVal[idx].replace('*', '').replace('!', '').split(' +/-')[0]
            # if np.abs(float(lastVal)) < 1.6:
            #     break
        # dumpfile.write('\\hline\n')

    return name, outVal, np.array(outputSigmas)


stringJubCustomizationCaltechT2 = '''
+RunAsOwner = True
+InteractiveUser = True
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"
+SingularityBindCVMFS = True
+MaxRuntime = 3600
+JobQueue   = "Short"
RequestDisk = 200000
RequestMemory = 2500
RequestCpus = 1
x509userproxy = $ENV(X509_USER_PROXY)
'''
