#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.stats as sps
from glob import glob
import pickle, re
import time
from array import array
from scipy.interpolate import interp1d
import multiprocessing
import matplotlib.pyplot as plt
import uproot as ur
import ROOT as rt
import root_numpy as rtnp
import ROOT.RooFit as rf
from scipy.special import erf
import sys, os
import itertools
sys.path.append('../lib')
if os.environ['CMSSW_VERSION'] != 'CMSSW_10_2_3':
    raise
from histo_utilities import create_TH1D, create_TH2D, std_color_list, SetMaxToMaxHist, make_ratio_plot
from cebefo_style import Set_2D_colz_graphics
from progressBar import ProgressBar

from analysis_utilities import drawOnCMSCanvas, extarct, extarct_multiple, createSel
from lumi_utilities import getLumiByTrigger

import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1

donotdelete = []

import argparse
parser = argparse.ArgumentParser(description='Script used to run trigger efficiencies.',
                                 epilog='Test example: ./triggerEfficiencies.py',
                                 add_help=True
                                 )
parser.add_argument ('--HELP', '-H', default=False, action='store_true', help='Print help message.')
parser.add_argument ('--version', '-v', default='test', help='Version name.')

parser.add_argument ('--dataset', '-d', type=str, default='MC', choices=['RD', 'MC'], help='Dataset to use.')
parser.add_argument ('--trigger', '-t', type=str, default='Mu7_IP4', choices=['Mu7_IP4', 'Mu9_IP6', 'Mu12_IP6'], help='Trigger to probe.')
parser.add_argument ('--tagTrigger', type=str, default='', choices=['Mu7_IP4', 'Mu9_IP6', 'Mu12_IP6'], help='Trigger of the tag muon.')

parser.add_argument ('--parallel', '-p', type=int, default=10, help='Number of parallel CPU to use.')

parser.add_argument ('--verbose', default=False, action='store_true', help='Verbose switch.')
parser.add_argument ('--submit', default=False, action='store_true', help='Submit a job instead of running the call interactively.')


args = parser.parse_args()
if args.HELP:
    parser.print_help()
    exit()

rt.gROOT.SetBatch(True)
plt.ioff()
plt.switch_backend('Agg')



rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
webFolder = '/storage/user/ocerri/public_html/BPH_RDst/triggerScaleFactors/'+args.version
if not os.path.isdir(webFolder):
    print 'Creating', webFolder
    os.makedirs(webFolder)
    os.system('cp '+webFolder+'/../index.php '+webFolder+'/')




cl = rt.TLine()
cl.SetLineColor(6)
cl.SetLineStyle(9)
cl.SetLineWidth(2)



colors = [rt.kBlack, rt.kAzure+1, rt.kRed-4, rt.kGreen+1, rt.kViolet-7]



branchesToLoad = ['mTag_pt', 'mTag_eta', 'mTag_phi', 'mTag_sigdxy',
                  'mTag_softID', 'mTag_tightID',
                  'mTag_HLT_Mu7_IP4', 'mTag_HLT_Mu9_IP6', 'mTag_HLT_Mu12_IP6',
                  'mProbe_pt', 'mProbe_eta', 'mProbe_phi', 'mProbe_sigdxy',
                  'mProbe_softID', 'mProbe_tightID',
                  'mProbe_HLT_Mu7_IP4', 'mProbe_HLT_Mu9_IP6', 'mProbe_HLT_Mu12_IP6',
                  'deltaR_tagProbe', 'massMuMu', 'vtx_isGood', 'massMuMu_refit',
                  'prescaleMu7_IP4', 'prescaleMu9_IP6', 'prescaleMu12_IP6', 'nVtx',
                 ]


def loadDF(loc, branches):
    dfL = []
    for l in loc:
        print l
        dfL.append(pd.DataFrame(rtnp.root2array(l, branches=branches)))
    if len(dfL) == 1:
        return dfL[0]
    else:
        return pd.concat(dfL)


class pileupReweighter(object):
    def __init__(self, mcSkimFile, cat, histoName='hAllNvtx', dataDate='200515'):
        loc = '/storage/user/ocerri/BPhysics/data/cmsRD/ParkingBPH{}/'+'Run2018D-05May2019promptD-v1_RDntuplizer_PrescaleVertices_{}_CAND.root'.format(dataDate)
        fAuxPileupRD = []

        hPileupTarget = None

        for i in range(1, 6):
            fAuxPileupRD.append(rt.TFile.Open(loc.format(i), 'READ'))
            if hPileupTarget is None:
                hPileupTarget = fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg).Clone()
            else:
                hPileupTarget.Add(fAuxPileupRD[-1].Get('nVtx/hNvtxPassed'+cat.trg))

        hPileupTarget.Scale(1./hPileupTarget.Integral())

        fAuxPileupMC = rt.TFile.Open(mcSkimFile, 'READ')
        hPileupGen = fAuxPileupMC.Get(histoName)

        weights = np.ones(hPileupGen.GetNbinsX())
        s = 0
        for i in range(weights.shape[0]):
            if hPileupGen.GetBinContent(i+1) == 0:
                continue
            weights[i] = hPileupTarget.GetBinContent(i+1)/(hPileupGen.GetBinContent(i+1)/hPileupGen.Integral())
            s += (hPileupGen.GetBinContent(i+1)/hPileupGen.Integral()) * weights[i]

        self.weightsPileupMC = weights/s

        for f in fAuxPileupRD + [fAuxPileupMC]:
            f.Close()

    def getPileupWeights(self, arrNvtx, selection=None):
        x = arrNvtx
        if not selection is None:
            x = x[selection]
        return self.weightsPileupMC[x.astype(np.int)]


class Bauble(object):
    pass

if args.dataset == 'RD':
    dataDir = '../data/cmsRD'
    RDdsLoc = glob(dataDir + '/ParkingBPH*/Run2018D-05May2019promptD-v1_RDntuplizer_TagAndProbeTrigger_210209_CAND.root')
    df = loadDF(RDdsLoc, branchesToLoad)
    print 'Data probe muons:', df.shape[0]
    CMS_lumi.extraText = "     Internal"
elif args.dataset == 'MC':
    mcDir = '../data/cmsMC_private/BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3'
    MCdsLoc = glob(mcDir + '/ntuples_TagAndProbeTrigger_Jpsi/merged/out_CAND.root')
    df = loadDF(MCdsLoc, branchesToLoad + ['sfMuonID'])
    print 'MC probe muons:', df.shape[0]

    aux = Bauble()
    aux.trg = args.trigger
    puRew = pileupReweighter(MCdsLoc[0], aux, histoName='TnP/hAllNvts')
    nMax = np.max(df['nVtx'])
    while nMax > (puRew.weightsPileupMC.shape[0] - 1):
        puRew.weightsPileupMC = np.append(puRew.weightsPileupMC, puRew.weightsPileupMC[-1])
        print args.trigger, puRew.weightsPileupMC.shape
    df['w'+args.trigger] = puRew.weightsPileupMC[df['nVtx'].astype(np.int)]
    df['w'] = df['sfMuonID']*df['w'+args.trigger]
    CMS_lumi.extraText = "     Simulation Internal"



def analyzeBin(idx, reloadDf=False, verbose=False):
    print idx, 'started'
    if reloadDf:
        raise
    psCut = []
    lim = {}
    selTot = None
    st = time.time()
    for n, i in idx.iteritems():
        lim[n] = [binning[n][i], binning[n][i+1]]
        if n=='eta':
            aux = np.abs(df['mProbe_'+n])
            sel = np.logical_and(aux > lim[n][0], aux < lim[n][1])
        else:
            sel = np.logical_and(df['mProbe_'+n] > lim[n][0], df['mProbe_'+n] < lim[n][1])
        if selTot is None:
            selTot = sel
        else:
            selTot = np.logical_and(sel, selTot)
    selTot = np.logical_and(selTot, df['prescale'+probeTrigger[4:]] > 0)
    selTot = np.logical_and(selTot, df['mProbe_softID'] > 0.5)
    selTot = np.logical_and(selTot, df['deltaR_tagProbe'] > 0.3)
    selTot = np.logical_and(selTot, df['vtx_isGood'] > 0.5)
    selTot = np.logical_and(selTot, np.abs(df['massMuMu_refit'] - 3.09691) < 0.1)
    if args.tagTrigger:
        selTot = np.logical_and(selTot, df['mTag_HLT_' + args.tagTrigger] == 1)


    if verbose:
        print ' --- Total ---'
    st = time.time()
    if args.dataset == 'RD':
        nSigTot = np.sum(selTot)
    else:
        nSigTot = np.sum(df['w'][selTot])
    if verbose:
        print 'Time: {:.1f} s'.format(time.time()-st)
        print ' --- Passed ---'
    st = time.time()
    selTot = np.logical_and(selTot, df['mProbe_' + probeTrigger] == 1)
    if args.dataset == 'RD':
        nSigPass = np.sum(selTot)
    else:
        nSigPass = np.sum(df['w'][selTot])

    if verbose:
        print 'Time: {:.1f} s'.format(time.time()-st)
    print idx, 'done'
    return idx, nSigTot, nSigPass


# # Run the fit in each bin
probeTrigger = 'HLT_'+args.trigger

if args.trigger == 'Mu7_IP4':
    binning = {'pt': array('d', [5.5, 6.5, 7, 7.1, 7.2, 7.3, 7.6, 8, 9, 9.2, 10, 12, 14]),
               'eta': array('d', [0, 0.4, 0.8, 1.5]),
               'sigdxy': array('d', [4, 5, 5.5, 6, 10, 20, 200])
              }
elif args.trigger == 'Mu9_IP6':
    binning = {'pt': array('d', [8.5, 9, 9.1, 9.2, 9.3, 9.6, 10.2, 11, 12, 12.2, 14]),
               'eta': array('d', [0, 0.4, 0.8, 1.5]),
               'sigdxy': array('d', [4, 6, 7, 7.5, 8, 10, 20, 200])
              }
elif args.trigger == 'Mu12_IP6':
    binning = {'pt': array('d', [11, 12, 12.2, 13, 14, 16, 18, 20, 22, 25, 28, 35]),
               'eta': array('d', [0, 0.4, 0.8, 1.5]),
               'sigdxy': array('d', [4, 6, 7, 8, 10, 20, 200])
              }


h2 = {}
for var, cat in itertools.product(['N', 'Chi2'], ['tot', 'pass']):
    h2[var+cat] = rt.TH3D('h2'+var+cat, '',
                          len(binning['pt'])-1, binning['pt'],
                          len(binning['sigdxy'])-1, binning['sigdxy'],
                          len(binning['eta'])-1, binning['eta'],)


start = time.time()
testOutput = analyzeBin({'pt': 2, 'sigdxy':2, 'eta':0}, verbose=True)
print testOutput
print 'Total time: {:.1f} mins'.format((time.time() - start)/60.)



inputs = []
for ipt in range(len(binning['pt'])-1):
    for iip in range(len(binning['sigdxy'])-1):
        for ieta in range(len(binning['eta'])-1):
            idx = {'pt': ipt, 'sigdxy':iip, 'eta': ieta}
            inputs.append(idx)
print 'Total bins:', len(inputs)


if args.parallel:
    N_max = min(args.parallel, max(1, multiprocessing.cpu_count() - 10))
    N_request = min(len(inputs), N_max)
    print 'Parallelization factor:', N_request
    p = multiprocessing.Pool(N_request)
    output = p.map(analyzeBin, inputs)
else:
    output = []
    for counter, i in enumerate(inputs):
        if not args.verbose:
            print '{}/{}'.format(counter+1, len(inputs))
        output.append(analyzeBin(i, verbose=args.verbose))



for idx, nSigTot, nSigPass in output:
    ip = idx['pt']+1
    ii = idx['sigdxy']+1
    ie = idx['eta']+1
    h2['Ntot'].SetBinContent(ip, ii, ie, nSigTot)
    h2['Ntot'].SetBinError(h2['Ntot'].GetBin(ip, ii, ie), np.sqrt(nSigTot))
    h2['Npass'].SetBinContent(ip, ii, ie, nSigPass)
    h2['Npass'].SetBinError(h2['Npass'].GetBin(ip, ii, ie), np.sqrt(nSigPass))



if not rt.TEfficiency.CheckConsistency(h2['Npass'], h2['Ntot']): raise
pEff = rt.TEfficiency(h2['Npass'], h2['Ntot'])
pEff.SetStatisticOption(rt.TEfficiency.kFCP)
pEff.SetNameTitle('eff_'+probeTrigger, 'Efficience for '+probeTrigger)

tf = rt.TFile('../data/calibration/triggerScaleFactors/{}_{}_{}.root'.format(probeTrigger, args.dataset, args.version), 'RECREATE')
pEff.Write()
for h in h2.values():
    h.Write()
tf.Close()


outCanvases = []
tdrstyle.setTDRStyle()
BRY_colors = [rt.kBlack, rt.kGray+1,
              rt.kBlue, rt.kAzure+1,
              rt.kViolet-7, rt.kMagenta-9, rt.kRed-4,
              rt.kOrange-3, rt.kYellow+7, rt.kGreen+1]
hRef = h2['Npass']
for iz in range(1, hRef.GetNbinsZ()+1):
    title = 'Efficiency {} {} ({:.1f} < |#eta| < {:.1f})'.format(probeTrigger, args.dataset, binning['eta'][iz-1], binning['eta'][iz])

    leg = rt.TLegend(0.7, 0.2, 0.98, 0.5)
    leg.SetLineWidth(0)
    leg.SetBorderSize(0)
    gr2draw = []

    for iy in range(1, hRef.GetNbinsY()+1):
        gr = rt.TGraphAsymmErrors()
        gr.SetName('gr_{}_{}'.format(iy,iz))
        for ix in range(1, hRef.GetNbinsX()+1):
            idx = pEff.GetGlobalBin(ix, iy, iz)
            # if iz == 3 and iy==1:
            #     print ix, iy, pEff.GetEfficiency(idx)
            x = binning['pt'][ix-1] + 0.5*(binning['pt'][ix] - binning['pt'][ix-1])
            gr.SetPoint(ix-1, x, pEff.GetEfficiency(idx))
            gr.SetPointError(ix-1, x-binning['pt'][ix-1], binning['pt'][ix]-x,
                             pEff.GetEfficiencyErrorLow(idx), pEff.GetEfficiencyErrorUp(idx)
                            )
        gr.SetLineColor(BRY_colors[iy-1])
        gr.SetMarkerColor(BRY_colors[iy-1])
        leg.AddEntry(gr, '{:.1f} < IP < {:.1f}'.format(binning['sigdxy'][iy-1], binning['sigdxy'][iy]), 'lep')
        gr2draw.append(gr)

    M = 1.2
    m = 0
    gr2draw[0].GetYaxis().SetRangeUser(m ,M)
    c = drawOnCMSCanvas(CMS_lumi, gr2draw, ['AP'] + (len(gr2draw)-1)*['P'], tag='_eff'+str(iz))
    gr2draw[0].GetYaxis().SetTitle('Efficiency')
    gr2draw[0].GetXaxis().SetTitle('Muon p_{T} [GeV]')
    leg.Draw()
#     c.SetLogx()
#     c.SetLogy()
#     m = 1e-3
#     M = 5
    gr2draw[0].GetXaxis().SetTitleOffset(1.1)

    trgThr = float(re.search(r'_Mu[0-9]+_', probeTrigger).group(0)[3:-1])
    cl.DrawLine(trgThr, m, trgThr, 1)

    l = rt.TLine()
    l.SetLineWidth(1)
    l.SetLineColor(rt.kGray)
    l.SetLineStyle(7)
    l.DrawLine(binning['pt'][0], 1, binning['pt'][-1], 1)

    rt.TLatex()
    text = rt.TLatex()
    text.SetTextAlign(22)
    text.SetTextSize(0.04)
    text.SetTextFont(42)
    text.DrawLatexNDC(0.6, 0.9, title);

    imgLoc = '../data/calibration/triggerScaleFactors/figEff/'
    c.SaveAs(imgLoc + probeTrigger+ '_' + args.dataset + '_eta{}_{}.png'.format(iz-1, args.version))
    c.SaveAs(webFolder+'/eff_'+probeTrigger+ '_' + args.dataset + '_eta{}_{}.png'.format(iz-1, args.version))
    outCanvases.append([c, gr2draw, leg])