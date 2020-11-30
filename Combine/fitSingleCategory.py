#!/usr/bin/env python

"""
######### Notes #########
To activate the environment run: cd ~/work/CMSSW_10_2_13/src/; cmsenv; cd ~/BPhysics/Combine/

######## Release notes #########
New from v16:
- Only positive q2 and M2_miss
- Revert to old B pT calibration
- Remove Muon pt correction

To do:
- Add combinatorial and fake from data template
"""


import sys, os, pickle, time, json, yaml, itertools, commands, re
from glob import glob
sys.path.append('../lib')
sys.path.append('../analysis')
from multiprocessing import Pool
from prettytable import PrettyTable
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from array import array

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

from progressBar import ProgressBar
from categoriesDef import categories as categoriesDef
from analysis_utilities import drawOnCMSCanvas, getEff, DSetLoader
from pT_calibration_reader import pTCalReader
from histo_utilities import create_TH1D, create_TH2D, std_color_list
from gridVarQ2Plot import plot_gridVarQ2, plot_SingleCategory, getControlXtitle, getControlSideText
from lumi_utilities import getLumiByTrigger
from combine_utilities import getUncertaintyFromLimitTree, dumpDiffNuisances, stringJubCustomizationCaltechT2, loadHisto4CombineFromRoot

import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "     Preliminary"
donotdelete = []

import argparse
parser = argparse.ArgumentParser()

parser.add_argument ('--category', '-c', type=str, default='low', choices=['single', 'low', 'mid', 'high'], help='Category')
parser.add_argument ('--useMVA', default=False, choices=[False, 'v0', 'v1'], help='Use MVA in the fit')
parser.add_argument ('--schemeFF', default='CLN', choices=['CLN', 'BLPR', 'NoFF'], help='Form factor scheme to use')
parser.add_argument ('--cardTag', '-v', default='vTest', help='Card name initial tag')

parser.add_argument ('--unblinded', default=False, action='store_true', help='Unblind the fit regions')
parser.add_argument ('--asimov', default=False, action='store_true', help='Use Asimov dataset insted of real data')
parser.add_argument ('--dataType', default=0, choices=[0,1,2], help='0: both, 1: only B0, 2: only anti-B0')
parser.add_argument ('--noMCstats', default=False, action='store_true', help='Do not include MC stat systematic')


availableSteps = ['histos', 'preFitPlots', 'card', 'workspace', 'bias', 'scan', 'fitDiag', 'postFitPlots', 'uncBreakdown', 'impacts', 'GoF']
defaultPipeline = ['histos', 'preFitPlots', 'card', 'workspace', 'scan', 'fitDiag', 'postFitPlots', 'uncBreakdown', 'GoF']
parser.add_argument ('--step', '-s', type=str, default=defaultPipeline, choices=availableSteps, help='Category(ies)', nargs='+')

parser.add_argument ('--forceRDst', default=False, action='store_true', help='Perform fit fixing R(D*) to 0.295')
# parser.add_argument ('--fitRegionsOnly', default=False, action='store_true', help='Include only regions used for fitting')
parser.add_argument ('--seed', default=6741, type=int, help='Seed used by Combine')
parser.add_argument ('--RDstLims', default=[0.1, 0.5], type=int, help='Initial boundaries of R(D*)', nargs='+')

# Scan options
parser.add_argument ('--maskScan', type=str, default=[], nargs='+', help='Channels to mask during likelyhood scan. If this list is non empty, the full card is used (default is fitregionsOnly)')
parser.add_argument ('--tagScan', type=str, default='')

# Impacts options
parser.add_argument ('--collectImpacts', default=False, action='store_true', help='Only collect impact fits which have been previously run')
parser.add_argument ('--subOnlyImpacts', default=False, action='store_true', help='Only submit impact fits, do not collect results')

# Goodness Of Fit options
parser.add_argument ('--algoGoF', type=str, default=['Sat', 'AD', 'KS'], choices=['Sat', 'AD', 'KS'], help='Algorithm to be used for the goodness of fit test', nargs='+')
parser.add_argument ('--maskEvalGoF', type=str, default=[], nargs='+', help='Additional channels to mask during GoF evaluation')
parser.add_argument ('--tagGoF', type=str, default='all')

parser.add_argument ('--showPlots', default=False, action='store_true', help='Show plots by setting ROOT batch mode OFF (default ON)')
parser.add_argument ('--showCard', default=False, action='store_true', help='Dump card on std outoput')

args = parser.parse_args()

category = args.category
schemeFF = args.schemeFF
if not args.showPlots:
    rt.gROOT.SetBatch(True)
    plt.ioff()

if len(args.RDstLims) > 2:
    raise
elif args.RDstLims[1] <= args.RDstLims[0]:
    raise

########################### ---- Define standards ----------- ########################

cat = categoriesDef[category]
binning = {'q2': array('d', [0, 3.5, 6, 9.4, 12])}
# binning['q2'] = array('d', [-2, 2.5, 6, 9.4, 12])

SM_RDst = 0.295
CMS_lumi.integrated_lumi = {'low':6.4, 'mid':20.7, 'high':26.4, 'single':20.7}[category] #fb^-1

FreeParFF = {
   'CLN': ['R0', 'eig1', 'eig2', 'eig3'],
   'BLPR': ['eig1', 'eig2', 'eig3', 'eig4', 'eig5', 'eig6'],
   'NoFF': []
}[schemeFF]

processOrder = ['tau', 'mu',
                'DstPip','DstPi0',
                'DstPipPi0','DstPipPim','DstPi0Pi0',
                'TauDstPi0', 'TauDstPip',
                'DstmDsp','DstmD0','DstmDp',
                'BpDstmHc','BmDstmHc','antiB0DstmHc']


SamplesB0 = ['mu', 'tau',
             'DstmD0', 'DstmDp', 'DstmDsp',
             'antiB0DstmHc',
             'DstPi0',
             'DstPipPim', 'DstPi0Pi0',
             'TauDstPi0'
            ]

SamplesBp = ['BpDstmHc', 'BmDstmHc',
             'DstPip',
             'DstPipPi0',
             'TauDstPip'
            ]


########################### --------------------------------- #########################
if args.asimov:
    CMS_lumi.extraText = "     Simulation Preliminary"

def createCardName(a):
    c = a.cardTag + a.category + '_' + a.schemeFF
    if a.useMVA:
        c += '_MVA'+useMVA
    if a.asimov:
        c += '_Asimov'
    elif a.dataType:
        c += '_onlyB0'
        if a.dataType == 2:
            c += 'bar'
    if not a.unblinded:
        c += '_blinded'
    if a.noMCstats:
        c += '_NoMCstats'
    return c

card_name = createCardName(args)
print 'Card name:', card_name

outdir = '/storage/user/ocerri/BPhysics/Combine/results/' + card_name
if not os.path.isdir(outdir):
    os.system('mkdir -p ' + outdir + '/fig')
card_location = '/storage/user/ocerri/BPhysics/Combine/cards/{}.txt'.format(card_name)
histo_file_dir = '/storage/user/ocerri/BPhysics/data/_root/histos4combine/'

webFolder = '/storage/user/ocerri/public_html/BPH_RDst/Combine/' + card_name
if not os.path.isdir(webFolder):
    os.makedirs(webFolder)
    os.system('cp '+webFolder+'/../index.php '+webFolder)


########################### -------- Create histrograms ------------------ #########################

def loadDatasets(loadRD, dataType):
    print 'Loading MC datasets'
    #They all have to be produced with the same pileup
    MCsample = {
    'tau' : DSetLoader('B0_TauNuDmst_PUc0'),
    'mu' : DSetLoader('B0_MuNuDmst_PUc0'),
    'DstmD0' : DSetLoader('B0_DstmD0_PUc0'),
    'DstmDp' : DSetLoader('B0_DstmDp_PUc0'),
    'DstmDsp' : DSetLoader('B0_DstmDsp_PUc0'),
    'BpDstmHc' : DSetLoader('Bp_DstmHc_PUc0'),
    'BmDstmHc' : DSetLoader('Bm_DstmHc_PUc0'),
    'antiB0DstmHc' : DSetLoader('antiB0_DstmHc_PUc0'),
    'DstPip' : DSetLoader('Bp_MuNuDstst_Pip_PUc0'),
    'DstPi0' : DSetLoader('B0_MuNuDstst_Pi0_PUc0'),
    'DstPipPi0' : DSetLoader('Bp_MuNuDstst_PipPi0_PUc0'),
    'DstPipPim' : DSetLoader('B0_MuNuDstst_PipPim_PUc0'),
    # 'DstPipPi0' : DSetLoader('Bp_MuNuDstPipPi0_PUc0'),
    # 'DstPipPim' : DSetLoader('B0_MuNuDstPipPim_PUc0'),
    'DstPi0Pi0' : DSetLoader('B0_MuNuDstst_Pi0Pi0_PUc0'),
    'TauDstPi0' : DSetLoader('B0_TauNuDstst_Pi0_PUc0'),
    'TauDstPip' : DSetLoader('Bp_TauNuDstst_Pip_PUc0')
    }

    dSet = {}
    dSetTkSide = {}
    for n, s in MCsample.iteritems():
        if not n in processOrder: raise
        dSet[n] = pd.DataFrame(rtnp.root2array(s.skimmed_dir + '/{}_corr.root'.format(cat.name)))
        dSetTkSide[n] = rtnp.root2array(s.skimmed_dir + '/{}_trkCtrl_corr.root'.format(cat.name))

    if loadRD:
        print 'Loading data datasets'
        dataDir = '/storage/user/ocerri/BPhysics/data/cmsRD'
        lumi_tot = 0

        #B0 data
        if dataType == 0 or dataType == 1:
            creation_date = '201101'
            locRD = dataDir+'/skimmed/B2DstMu_B0_{}_{}'.format(creation_date, cat.name)
            dSet['data'] = pd.DataFrame(rtnp.root2array(locRD + '_corr.root'))
            dSetTkSide['data'] = pd.DataFrame(rtnp.root2array(locRD + '_trkCtrl_corr.root'))
            datasets_loc = glob(dataDir + '/ParkingBPH*/*RDntuplizer_B2DstMu_{}_CAND.root'.format(creation_date))
            lumi_tot = getLumiByTrigger(datasets_loc, cat.trg, verbose=True)

        #Anti-B0 data
        if dataType == 0 or dataType == 2:
            creation_date = '201122'
            locRD = dataDir+'/skimmed/B2DstMu_antiB0_{}_{}'.format(creation_date, cat.name)

            aux = pd.DataFrame(rtnp.root2array(locRD + '_corr.root'))
            auxTk = pd.DataFrame(rtnp.root2array(locRD + '_trkCtrl_corr.root'))
            if dataType == 0:
                dSet['data'] = pd.concat([dSet['data'], aux])
                dSetTkSide['data'] = pd.concat([dSetTkSide['data'], auxTk])
            else:
                dSet['data'] = dSet['data']
                dSetTkSide['data'] = dSetTkSide['data']

        if not lumi_tot == 0:
            CMS_lumi.integrated_lumi = lumi_tot

    return MCsample, dSet, dSetTkSide

def computeBrVarWeights(ds, selItems={}, relScale=0.2, keepNorm=False):
    sel = np.ones_like(ds['mu_pt']).astype(np.bool)
    for var, val in selItems.iteritems():
        sel = np.logical_and(ds[var].astype(np.int) == val, sel)
    w = np.ones_like(sel)
    up = np.where(sel, 1.+relScale, 1.)
    down = np.where(sel, max(0, 1.-relScale), 1.)
    if keepNorm:
        up = float(up.shape[0])/np.sum(up)
        down = float(down.shape[0])/np.sum(down)
    return w, up, down

def computeWidthVarWeights(ds, selItems=[], relScale=0.1): #Gamma modification factor
    # selItems=[ [pdgId, mass, Gamma] ]
    up = np.ones_like(ds['mu_pt'])
    down = np.ones_like(ds['mu_pt'])
    for pdgId, mass, gamma in selItems:
        # print pdgId, mass, gamma
        gUp = gamma*(1+relScale)
        gDown = gamma*(1-relScale)

        dx2 = np.clip(np.square(ds['MC_MassCharmedBDaughter'] - mass), 0, 9*(gamma**2))
        wUp = ((dx2 + gamma**2)*gUp)/(gamma*(dx2 + gUp**2))
        wDown = ((dx2 + gamma**2)*gDown)/(gamma*(dx2 + gDown**2))

        sel = np.abs(ds['MC_DstMotherPdgId'].astype(np.int)) == np.abs(pdgId)
        up = np.where(sel, wUp, up)
        down = np.where(sel, wDown, down)

    w = np.ones_like(sel)
    return w, up, down

def computeTksPVweights(ds, relScale=0.05, centralVal=0.39/0.10):
    selPdgID0 = np.logical_and(np.abs(ds['MC_tkMotherPdgId_0']) < 6, ds['MC_tkMotherPdgId_0'] != 0)
    selPdgID0 = np.logical_or(selPdgID0, ds['MC_tkMotherPdgId_0']==2212)
    selPdgID0 = np.logical_and(selPdgID0, ds['MC_tkFlag_0'] == 1)
    selPdgID1 = np.logical_and(np.abs(ds['MC_tkMotherPdgId_1']) < 6, ds['MC_tkMotherPdgId_1'] != 0)
    selPdgID1 = np.logical_or(selPdgID1, ds['MC_tkMotherPdgId_1']==2212)
    selPdgID1 = np.logical_and( selPdgID1, ds['MC_tkFlag_1'] == 1)
    exponent = selPdgID0.astype(np.int) + selPdgID1.astype(np.int)
    w = np.power(centralVal, exponent)
    up = np.power(centralVal*(1+relScale), exponent)/w
    down = np.power(centralVal*(1-relScale), exponent)/w
    return w, up, down

def createHistograms():
    MCsample, dSet, dSetTkSide = loadDatasets(not args.asimov, args.dataType)

    ######################################################
    ########## Load calibrations
    ######################################################
    from pileup_utilities import pileupReweighter
    puReweighter = pileupReweighter(MCsample['mu'].skimmed_dir + '/{}_corr.root'.format(cat.name), cat)

    dataDir = '/storage/user/ocerri/BPhysics/data'
    decayBR = pickle.load(open(dataDir+'/forcedDecayChannelsFactors.pickle', 'rb'))

    loc = dataDir+'/calibration/triggerScaleFactors/'
    fTriggerSF = rt.TFile.Open(loc + 'HLT_' + cat.trg + '_SF_v3.root', 'READ')
    hTriggerSF = fTriggerSF.Get('hSF_HLT_' + cat.trg)
    def computeTrgSF(ds, selection=None):
        trgSF = np.ones_like(ds['q2'])
        trgSFUnc = np.zeros_like(ds['q2'])
        ptmax = hTriggerSF.GetXaxis().GetXmax() - 0.01
        ipmax = hTriggerSF.GetYaxis().GetXmax() - 0.01
        etamax = hTriggerSF.GetZaxis().GetXmax() - 0.01
        x = np.column_stack((ds['mu_pt'], ds['mu_eta'], ds['mu_sigdxy']))
        if not selection is None:
            x = x[selection]
        for i, (pt, eta, ip) in enumerate(x):
            ix = hTriggerSF.GetXaxis().FindBin(min(ptmax, pt))
            iy = hTriggerSF.GetYaxis().FindBin(min(ipmax, ip))
            iz = hTriggerSF.GetZaxis().FindBin(min(etamax, np.abs(eta)))
            trgSF[i] = hTriggerSF.GetBinContent(ix, iy, iz)
            ib = hTriggerSF.GetBin(ix, iy, iz)
            trgSFUnc[i] = hTriggerSF.GetBinError(ib)
            if trgSF[i] == 0:
                print pt, ip, np.abs(eta)
                raise
        # Divide them for the weight so later you can simply multiply back to get the value
        up = 1 + trgSFUnc/trgSF
        down = 1 - trgSFUnc/trgSF
        return trgSF, up, down

    fMuonIDSF = rt.TFile.Open(dataDir+'/calibration/muonIDscaleFactors/Run2018ABCD_SF_MuonID_Jpsi.root', 'READ')
    hMuonIDSF = fMuonIDSF.Get('NUM_SoftID_DEN_genTracks_pt_abseta')
    def computeMuonIDSF(ds, selection=None):
        muonSF = np.ones_like(ds['q2'])
        muonSFUnc = np.zeros_like(ds['q2'])
        ptmax = hMuonIDSF.GetXaxis().GetXmax() - 0.01
        etamax = hMuonIDSF.GetYaxis().GetXmax() - 0.01
        x = np.column_stack((ds['MC_mu_pt'], ds['MC_mu_eta']))
        if not selection is None:
            x = x[selection]
        for i, (pt, eta) in enumerate(x):
            ix = hMuonIDSF.GetXaxis().FindBin(min(pt, ptmax))
            if ix == 0: ix = 1 #Remove underflows (Meaning that the MC matching failed)
            iy = hMuonIDSF.GetYaxis().FindBin(min(np.abs(eta), etamax))
            muonSF[i] = hMuonIDSF.GetBinContent(ix, iy)
            muonSFUnc[i] = hMuonIDSF.GetBinError(hMuonIDSF.GetBin(ix, iy))
            if muonSF[i] == 0:
                print pt, eta
                print ix, iy
                raise
        up = 1 + muonSFUnc/muonSF
        down = 1 - muonSFUnc/muonSF
    #     print np.column_stack((muonSF, up, down, muonSF*up, muonSF*down))
    #     raise
        return muonSF, up, down


    aux = 'Low' if category == 'single' else cat.name
    cal_pT_B0 = pTCalReader(calibration_file=dataDir+'/calibration/B0pTspectrum/pwWeights_{}.txt'.format(aux))
    # cal_pT_B0 = pTCalReader(calibration_file=dataDir+'/calibration/B0pTspectrum/polyCoeffWeights_{}.pkl'.format(aux))
    cal_pT_Bp = pTCalReader(calibration_file=dataDir+'/calibration/Bcharged_pTspectrum/pwWeights_{}.txt'.format(aux))
    # cal_pT_mu = pTCalReader(calibration_file=dataDir+'/calibration/MuonPtSpectrum/polyCoeffWeights_{}.pkl'.format(aux))
    def computePtWeights(ds, var, tag, cal_pT):
        if cal_pT.kind == 'poly':
            # The denominator (sum of weights) for this weights is not known but it cancel out in the ratio
            w = cal_pT.getWeights(ds[var], shape=0)
            if np.sum(w==0):
                print np.sum(w==0)
                raise

            varDic = {}
            for iShape in range(1, cal_pT.nVar+1):
                varDic[tag+'_lam{}Down'.format(iShape)] = cal_pT.getWeights(ds[var], shape=-iShape)/w
                varDic[tag+'_lam{}Up'.format(iShape)] = cal_pT.getWeights(ds[var], shape=+iShape)/w
            return w, varDic
        elif cal_pT.kind == 'ratio':
            w = cal_pT.f['C'](ds[var])
            if np.sum(w==0):
                print np.sum(w==0)
                raise
            up = cal_pT.f['Up'](ds[var])/w
            down = cal_pT.f['Down'](ds[var])/w
            return w, up, down
        else:
            print 'Unknown calibration'
            raise

    pWeightsEta = pickle.load(open(dataDir+'/calibration/B0pTspectrum/etaWeights_poly_{}.p'.format(cat.name), 'rb'))
    def computeB0etaWeights(ds):
        w = np.polyval(p=pWeightsEta, x=ds['B_eta'])
        return np.clip(w, a_min=0.5, a_max=1.5)

    if args.useMVA:
        fname = dataDir+'../plot_scripts/kinObsMVA/clfGBC_tauVall_{}{}.p'.format(args.useMVA, cat.name)
        clfGBC = pickle.load(open(fname, 'rb'))

        def computeVarMVA(ds):
            if args.useMVA == 'v0':
                aux = np.column_stack((ds['q2'], ds['Est_mu'], ds['M2_miss']))
            elif args.useMVA == 'v1':
                ds['pt_vis'] = ds['B_pt']*ds['mass_D0pismu']/5.27963
                aux = np.column_stack((ds['q2'], ds['Est_mu'], ds['M2_miss'],
                                       ds['pt_vis'], ds['mass_D0pismu'],
                                       ds['B_eta']
                                             ))
            else: raise
            p = clfGBC.predict_proba(aux)
            return p[:,1]

    histo = {}
    eventCountingStr = {}
    RDoMC_normRatio = 6 # Difference in Pythia prediction between Hard and Soft QCD production
    RDoMC_normRatio *= 2 #Both charge signs

    ######################################################
    ########## Signal region
    ######################################################
    n_q2bins = len(binning['q2'])-1
    binning['M2_miss'] = [
    #         array('d', [-2.5] + list(np.arange(-1.8, -0.2, 0.4)) + [-0.2, 0., 0.2, 0.6, 8] ),
            array('d', [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 4] ),
    #         array('d', [-2.5] + list(np.arange(-1.8, -0.1, 0.4)) + [-0.1, 0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.0, 0.4)) + [8] ),
            array('d', [0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.5, 0.2)) + [8] ),

    #         array('d', [-2.5] + list(np.arange(-1.8, 5.6, 0.4)) + [8] ),
            array('d', list(np.arange(0, 6, 0.2)) + [8] ),

    #         array('d', [-2.5] + list(np.arange(-1.8, 7.6, 0.4)) + [8] ),
            array('d', list(np.arange(0, 7.8, 0.2)) + [8] ),
        ]
    binning['Est_mu'] = [
            array('d', [0.3] + list(np.arange(0.5, 2.3, 0.05)) + [2.3] ),
            array('d', [0.3] + list(np.arange(0.5, 2.2, 0.05)) + [2.2] ),
            array('d', [0.3] + list(np.arange(0.5, 2.1, 0.05)) + [2.1] ),
            [24, 0.3, 2.0],
        ]
    binning['mu_pt'] = n_q2bins*[{'low': array('d', list(np.arange(7, 9, 0.05))+[9] ),
                                  'mid': array('d', list(np.arange(9, 12, 0.05)) +[12] ),
                                  'high': array('d', list(np.logspace(np.log10(12), np.log10(50), 40))
                                 }[category]]

    binning['Dst_pt'] = n_q2bins*[array('d', list(np.arange(5, 30, 1)) )]
    binning_2D = [
        [
    #         array('d', [-2.5] + list(np.arange(-1.8, -0.2, 0.4)) + [-0.2, 0., 0.2, 0.6, 8] ),
            array('d', [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 4] ),

            array('d', [0.3] + list(np.arange(0.7, 2.3, 0.3)) + [2.3] )
        ],
        [
    #         array('d', [-2.5] + list(np.arange(-1.8, 3.0, 0.4)) + [8] ),
            array('d', [0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.0, 0.2)) + [8] ),

            array('d', [0.3] + list(np.arange(0.7, 2.2, 0.3)) + [2.2] )
        ],
        [
    #         array('d', [-2.5] + list(np.arange(-1.8, 5.6, 0.4)) + [8] ),
            array('d', list(np.arange(0, 5.6, 0.4)) + [8] ),

            array('d', [0.3] + list(np.arange(0.5, 2.1, 0.3)) + [2.1] )
        ],
        [
    #         array('d', [-2.5] + list(np.arange(-1.8, 7.6, 0.4)) + [8] ),
            array('d', list(np.arange(0, 7.6, 0.4)) + [8] ),

            array('d', list(np.linspace(0.3, 2.0, 10)) )
        ]

    ]
    binning['B_pt'] = {'low': array('d', list(np.arange(10, 75, 2)) ), 'mid': array('d', list(np.arange(14, 90, 2)) ), 'high': array('d', list(np.arange(18, 110, 2)))}[category]
    binning['B_eta'] = array('d', list(np.arange(-1.9, 1.91, 0.05)) )

    if args.unblinded:
        binning['MVA'] = array('d', list(np.arange(0., 0.83, 0.02)))
    else:
        binning['MVA'] = array('d', list(np.arange(0., 0.51, 0.03)))


    totalCounting = [0,0]
    print '---------> Fill signal region histograms'
    for n in processOrder:
        ds = dSet[n]
        if n == 'data': continue
        print '\n----------->', n, '<-------------'
        sMC = MCsample[n]

        nTotSelected = ds['q2'].shape[0]
        print 'N tot selected: {:.1f}k'.format(1e-3*nTotSelected)
        totalCounting[1] += 1e-3*nTotSelected
        nGenExp = sMC.effMCgen['xsec'][0] * CMS_lumi.integrated_lumi * RDoMC_normRatio
        eff = [1, 0]
        for f, df in [sMC.effMCgen['effGEN'], decayBR[n], sMC.effCand['effCAND'], sMC.getSkimEff(cat.name+'_corr')]:
            eff[0] *= f
            eff[1] += np.square(df/f)
        eff[1] = eff[0] * np.sqrt(eff[1])
        nTotExp = nGenExp*eff[0]
        print 'N tot expected (before weights): {:.2f}k'.format(1e-3*nTotExp)

        wVar = {}
        weights = {}

        print 'Including pileup reweighting'
        weights['pileup'] = puReweighter.weightsPileupMC[ds['N_vtx'].astype(np.int)]

        print 'Including trigger corrections'
        weights['trg{}SF'.format(cat.trg)], wVar['trg{}SFUp'.format(cat.trg)], wVar['trg{}SFDown'.format(cat.trg)] = computeTrgSF(ds)
    #     weights['trg{}SF'.format(cat.trg)], _, _ = computeTrgSF(ds)

        print 'Including muon ID corrections'
    #     weights['muonIdSF'], wVar['muonIdSFUp'], wVar['muonIdSFDown'] = computeMuonIDSF(ds)
        weights['muonIdSF'], _, _ = computeMuonIDSF(ds)

    #     print 'Including muon pT corrections'
    #     weights['MuPt'], auxVarDic = computePtWeights(ds, 'mu_pt', 'MuPt', cal_pT_mu)
    #     wVar.update(auxVarDic)

        # B phase space corrections
        weights['etaB'] = computeB0etaWeights(ds)
        if n in SamplesB0:
            print 'Including B0 pT corrections'
            if cal_pT_B0.kind == 'ratio':
                weights['B0pT'], wVar['B0pTUp'], wVar['B0pTDown'] = computePtWeights(ds, 'MC_B_pt', None, cal_pT_Bp)
            else:
                weights['B0pT'], auxVarDic = computePtWeights(ds, 'MC_B_pt', 'B0pT', cal_pT_B0)
                wVar.update(auxVarDic)
        if n in SamplesBp:
            print 'Including B +/- pT corrections'
            weights['BpPt'], wVar['BpPtUp'], wVar['BpPtDown'] = computePtWeights(ds, 'MC_B_pt', None, cal_pT_Bp)
        # Hammer corrections to the FF
        if n in ['mu', 'tau'] and schemeFF != 'NoFF':
            print 'Including FF corrections (Hammer)'
            weights['B2DstFF'] = ds['wh_'+schemeFF+'Central']*sMC.effCand['rate_den']/sMC.effCand['rate_'+schemeFF+'Central']
            for nPar in FreeParFF:
                for var in ['Up', 'Down']:
                    tag = schemeFF + nPar + var
                    wVar['B2Dst'+tag] = ds['wh_'+tag]/ds['wh_'+schemeFF+'Central']
                    wVar['B2Dst'+tag] *= sMC.effCand['rate_'+schemeFF+'Central']/sMC.effCand['rate_' + tag]
        #Dstst resonance mix
        if n == 'DstPip' or n =='TauDstPip':
            print 'Including D**->D*Pi width variations'
            _, wVar['fDststWideUp'], wVar['fDststWideDown'] = computeBrVarWeights(ds, {'MC_munuSisterPdgId': -20423}, 0.6/2.7, keepNorm=True) #Gamma 14 pdg 2020

            _, wVar['D2420_10WidthUp'], wVar['D2420_10WidthDown'] = computeWidthVarWeights(ds, selItems=[[10423, 2.422, 0.020]], relScale=0.15)
            # _, wVar['D2430_10WidthUp'], wVar['D2430_10WidthDown'] = computeWidthVarWeights(ds, selItems=[[20423, 2.445, 0.250]], relScale=0.1)
            _, wVar['D2460_1StWidthUp'], wVar['D2460_1StWidthDown'] = computeWidthVarWeights(ds, selItems=[[425, 2.461, 0.043]], relScale=0.15)
        if n == 'DstPipPim' or n == 'DstPi0Pi0':
            print 'Including D**->D*PiPi width variations'
            widthMods = [[100413, 2.640, 0.200]]
            _, wVar['DstPiPiWidthUp'], wVar['DstPiPiWidthDown'] = computeWidthVarWeights(ds, selItems=widthMods, relScale=0.1)
        #Hc mix variations
        if n == 'DstmD0':
            _, wVar['BrB02DstD0KpUp'], wVar['BrB02DstD0KpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 421, 'MC_DstSisterPdgId_light': 321}, 0.21/2.47) #Gamma 169 pdg 2020
            _, wVar['BrB02DstD0KstpUp'], wVar['BrB02DstD0KstpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 421, 'MC_DstSisterPdgId_light': 323}, 0.5) # Guess
            _, wVar['BrB02DstDst0KpUp'], wVar['BrB02DstDst0KpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 423, 'MC_DstSisterPdgId_light': 321}, 0.09/1.06) #Gamma 170 pdg 2020
            _, wVar['BrB02DstDst0KstpUp'], wVar['BrB02DstDst0KstpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 423, 'MC_DstSisterPdgId_light': 323}, 0.5) # Guess
            _, wVar['BrB02DstDstpK0Up'], wVar['BrB02DstDstpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 311}, 0.5/5.3) #Gamma 173 pdg 2020
            _, wVar['BrB02DstDstpKst0Up'], wVar['BrB02DstDstpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
        if n == 'DstmDp':
            _, wVar['BrB02DstDpK0Up'], wVar['BrB02DstDpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 411, 'MC_DstSisterPdgId_light': 311}, 0.5/3.2) #Gamma 172 pdg 2020
            _, wVar['BrB02DstDpKst0Up'], wVar['BrB02DstDpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 411, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
            _, wVar['BrB02DstDstpK0Up'], wVar['BrB02DstDstpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 311}, 0.2/2.7) #Gamma 173 pdg 2020
            _, wVar['BrB02DstDstpKst0Up'], wVar['BrB02DstDstpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
        if n == 'DstmDsp':
            _, wVar['BrB02DstDsUp'], wVar['BrB02DstDsDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 431}, 1.1/8.0) #Gamma 83 pdg 2020
            _, wVar['BrB02DstDsstUp'], wVar['BrB02DstDsstDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 433}, .14/1.77) #Gamma 85 pdg 2020
            _, wVar['BrB02DstDs0stUp'], wVar['BrB02DstDs0stDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 10431}, .6/1.5) #Gamma 95 pdg 2020

        print 'Computing total weights'
        weightsCentral = np.ones_like(ds['q2'])
        for w in weights.values():
            weightsCentral *= w
        print 'N tot expected (after weights): {:.3f}k'.format(1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected)
        totalCounting[0] += 1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected
        evCountStr = '{:.2f} ({:.2f})'.format(1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected, 1e-3*nTotSelected)
        eventCountingStr[n] = evCountStr
        wVar[''] = np.ones_like(weightsCentral)


        if args.useMVA:
            print 'Evaluating MVA'
            if not 'MVA' in histo.keys():
                histo['MVA'] = {}
            sMVA = computeVarMVA(ds)
            for name_wVar, v_wVar in wVar.iteritems():
                    h_name = n
                    if not name_wVar == '':
                        h_name += '__' + name_wVar
                    w = weightsCentral*v_wVar
                    scale = nTotExp/nTotSelected
                    histo['MVA'][h_name] = create_TH1D(sMVA, name=h_name, weights=w, scale_histo=scale,
                                                       binning=binning['MVA'], opt='underflow,overflow')

        # Variables to be broken in q2 bins
        for i_q2 in range(n_q2bins):
            q2_l = binning['q2'][i_q2]
            q2_h = binning['q2'][i_q2 + 1]
            sel_q2 = np.logical_and(ds['q2'] > q2_l, ds['q2'] <= q2_h)
            name2D = 'h2D_q2bin'+str(i_q2)
            if not name2D in histo.keys():
                    histo[name2D] = {}
            for var in ['M2_miss', 'Est_mu', 'mu_pt', 'Dst_pt']:
                cat_name = var+'_q2bin'+str(i_q2)

                if not cat_name in histo.keys():
                    histo[cat_name] = {}

                for name_wVar, v_wVar in wVar.iteritems():
                    h_name = n
                    if not name_wVar == '':
                        h_name += '__' + name_wVar
                    w = weightsCentral*v_wVar
                    scale = nTotExp/nTotSelected
                    histo[cat_name][h_name] = create_TH1D(
                                                          ds[var][sel_q2],
                                                          name=h_name, title=h_name,
                                                          binning=binning[var][i_q2],
                                                          opt='underflow,overflow',
                                                          weights=w[sel_q2], scale_histo=scale,
                                                          )
                    if not var == 'M2_miss':
                        continue
                    auxS = np.column_stack((ds['M2_miss'][sel_q2], ds['Est_mu'][sel_q2]))
                    histo[name2D][h_name] = create_TH2D(
                                                          auxS,
                                                          name=h_name, title=h_name,
                                                          binning=binning_2D[i_q2],
                                                          weights=w[sel_q2], scale_histo=scale,
                                                       )
        # Variables in the whole spectrum
        for var in ['B_pt', 'B_eta']:
            if not var in histo.keys():
                histo[var] = {}
            for name_wVar, v_wVar in wVar.iteritems():
                h_name = n
                if not name_wVar == '':
                    h_name += '__' + name_wVar
                w = weightsCentral*v_wVar
                scale = nTotExp/nTotSelected
                histo[var][h_name] = create_TH1D(ds[var], name=h_name, weights=w, scale_histo=scale,
                                                    binning=binning[var], opt='underflow,overflow')

    evCountStr = '{:.1f} ({:.1f})'.format(*totalCounting)
    eventCountingStr['tot'] = evCountStr

    # Do the unrolling
    nonEmptyIdxs = []
    unrollingCutOff = 3
    for i_q2 in range(len(binning['q2'])-1):
        nonEmptyIdxs.append([])
        name2D = 'h2D_q2bin'+str(i_q2)
        hSum = None
        nDroppedBins = 0
        nExpectedDroppedEvents = 0
        for key, hN in histo[name2D].iteritems():
            if '__' in key:
                continue
            if hSum is None:
                hSum = hN.Clone('hSum_'+str(i_q2))
            else:
                scale = SM_RDst if 'tau' in n else 1.
                hSum.Add(hN, scale)
        for ix in range(1, hSum.GetNbinsX()+1):
            for iy in range(1, hSum.GetNbinsY()+1):
                if hSum.GetBinContent(ix, iy) > unrollingCutOff:
                    nonEmptyIdxs[i_q2].append([ix, iy])
                else:
                    nDroppedBins += 1
                    nExpectedDroppedEvents += hSum.GetBinContent(ix, iy)
        print 'Dropped bins:', nDroppedBins
        print 'Expected dropped candidates:', nExpectedDroppedEvents

        nameU = 'Unrolled_q2bin'+str(i_q2)
        histo[nameU] = {}
        validBins = nonEmptyIdxs[i_q2]
        for n, h in histo[name2D].iteritems():
            hUnrolled = rt.TH1D(h.GetName(), h.GetTitle(), len(validBins), 0.5, len(validBins)+0.5)
            for i, (ix, iy) in enumerate(validBins):
                hUnrolled.SetBinContent(i+1, h.GetBinContent(ix, iy))
                hUnrolled.SetBinError(i+1, h.GetBinError(ix, iy))
            histo[nameU][n] = hUnrolled


    ######################################################
    ########## Control region
    ######################################################


    sideSelecton = {}
    sideVar = {}

    def selfun__TkPlus(ds):
        sel = np.logical_and(ds['N_goodAddTks'] == 1, ds['tkCharge_0'] > 0)
        return sel
    sideSelecton['AddTk_p_mHad'] = selfun__TkPlus
    sideVar['AddTk_p_mHad'] = 'massHadTks'
    # sideVar['AddTk_p_mHad'] = 'massHadTks_DstMassConstraint'
    binning['AddTk_p_mHad'] = [35, 2.13, 2.83]

    def selfun__TkMinus(ds):
        sel = np.logical_and(ds['N_goodAddTks'] == 1, ds['tkCharge_0'] < 0)
        return sel
    sideSelecton['AddTk_m_mHad'] = selfun__TkMinus
    sideVar['AddTk_m_mHad'] = 'massHadTks'
    binning['AddTk_m_mHad'] = [30, 2.1, 3.3]


    def selfun__TkPlusMinus(ds):
        sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == 0, ds['N_goodAddTks'] == 2)
        sel = np.logical_and(ds['massVisTks'] < 5.3, sel)
        return sel
    sideSelecton['AddTk_pm_mVis'] = selfun__TkPlusMinus
    sideVar['AddTk_pm_mVis'] = 'massVisTks'
    binning['AddTk_pm_mVis'] = array('d', [2.8] + list(np.arange(3., 5.3, 0.1)) + [5.3] )
    # binning['AddTk_pm_mVis'] = array('d', list(np.arange(3., 5.3, 0.1)) + [5.3] )


    sideSelecton['AddTk_pm_mHad'] = selfun__TkPlusMinus
    sideVar['AddTk_pm_mHad'] = 'massHadTks'
    binning['AddTk_pm_mHad'] = [30, 2.3, 3.75]

    def selfun__TkMinusMinus(ds):
        sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == -2, ds['N_goodAddTks'] == 2)
        sel = np.logical_and(ds['massVisTks'] < 5.3, sel)
        return sel
    sideSelecton['AddTk_mm_mHad'] = selfun__TkMinusMinus
    sideVar['AddTk_mm_mHad'] = 'massHadTks'
    binning['AddTk_mm_mHad'] = [15, 2.25, 3.6]

    def selfun__TkPlusPlus(ds):
        sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == +2, ds['N_goodAddTks'] == 2)
        sel = np.logical_and(ds['massVisTks'] < 5.3, sel)
        return sel
    sideSelecton['AddTk_pp_mHad'] = selfun__TkPlusPlus
    sideVar['AddTk_pp_mHad'] = 'massHadTks'
    binning['AddTk_pp_mHad'] = [15, 2.25, 3.75]

    # Fill control regions histrograms
    for k in sideSelecton.keys(): histo[k] = {}
    totalCounting = {}
    for n in processOrder:
        ds = dSetTkSide[n]
        if n == 'data': continue
        print '\n----------->', n, '<-------------'
        sMC = MCsample[n]
        wVar = {}
        weights = {}

        print 'Including pileup reweighting'
        weights['pileup'] = puReweighter.weightsPileupMC[ds['N_vtx'].astype(np.int)]
        print 'Including trigger corrections'
        weights['trg{}SF'.format(cat.trg)], wVar['trg{}SFUp'.format(cat.trg)], wVar['trg{}SFDown'.format(cat.trg)] = computeTrgSF(ds)
        print 'Including muon ID corrections'
        weights['muonIdSF'], _, _ = computeMuonIDSF(ds)
    #     print 'Including muon pT corrections'
    #     weights['MuPt'], auxVarDic = computePtWeights(ds, 'mu_pt', 'MuPt', cal_pT_mu)
    #     wVar.update(auxVarDic)

        if n in SamplesB0:
            print 'Including B0 pT corrections'
            if cal_pT_B0.kind == 'ratio':
                weights['B0pT'], wVar['B0pTUp'], wVar['B0pTDown'] = computePtWeights(ds, 'MC_B_pt', None, cal_pT_Bp)
            else:
                weights['B0pT'], auxVarDic = computePtWeights(ds, 'MC_B_pt', 'B0pT', cal_pT_B0)
                wVar.update(auxVarDic)
        if n in SamplesBp:
            print 'Including B +/- pT corrections'
            weights['BpPt'], wVar['BpPtUp'], wVar['BpPtDown'] = computePtWeights(ds, 'MC_B_pt', None, cal_pT_Bp)
        if n in ['mu', 'tau'] and schemeFF != 'NoFF':
            print 'Including FF corrections (Hammer)'
            weights['B2DstFF'] = ds['wh_'+schemeFF+'Central']*sMC.effCand['rate_den']/sMC.effCand['rate_'+schemeFF+'Central']
            for nPar in FreeParFF:
                for var in ['Up', 'Down']:
                    tag = schemeFF + nPar + var
                    wVar['B2Dst'+tag] = ds['wh_'+tag]/ds['wh_'+schemeFF+'Central']
                    wVar['B2Dst'+tag] *= sMC.effCand['rate_'+schemeFF+'Central']/sMC.effCand['rate_' + tag]
        #Dstst resonance mix
        if n == 'DstPip' or n == 'TauDstPip':
            print 'Including D**->D*Pi width variations'
            _, wVar['fDststWideUp'], wVar['fDststWideDown'] = computeBrVarWeights(ds, {'MC_munuSisterPdgId': -20423}, 0.6/2.7, keepNorm=True) #Gamma 14 pdg 2020

            _, wVar['D2420_10WidthUp'], wVar['D2420_10WidthDown'] = computeWidthVarWeights(ds, selItems=[[10423, 2.422, 0.020]], relScale=0.15)
            # _, wVar['D2430_10WidthUp'], wVar['D2430_10WidthDown'] = computeWidthVarWeights(ds, selItems=[[20423, 2.445, 0.250]], relScale=0.1)
            _, wVar['D2460_1StWidthUp'], wVar['D2460_1StWidthDown'] = computeWidthVarWeights(ds, selItems=[[425, 2.461, 0.043]], relScale=0.15)
        if n == 'DstPipPim' or n == 'DstPi0Pi0':
            print 'Including D**->D*PiPi width variations'
            widthMods = [[100413, 2.640, 0.200]]
            _, wVar['DstPiPiWidthUp'], wVar['DstPiPiWidthDown'] = computeWidthVarWeights(ds, selItems=widthMods, relScale=0.1)
        #Hc mix variations
        if n == 'DstmD0':
            _, wVar['BrB02DstD0KpUp'], wVar['BrB02DstD0KpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 421, 'MC_DstSisterPdgId_light': 321}, 0.21/2.47) #Gamma 169 pdg 2020
            _, wVar['BrB02DstD0KstpUp'], wVar['BrB02DstD0KstpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 421, 'MC_DstSisterPdgId_light': 323}, 0.5) # Guess
            _, wVar['BrB02DstDst0KpUp'], wVar['BrB02DstDst0KpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 423, 'MC_DstSisterPdgId_light': 321}, 0.09/1.06) #Gamma 170 pdg 2020
            _, wVar['BrB02DstDst0KstpUp'], wVar['BrB02DstDst0KstpDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 423, 'MC_DstSisterPdgId_light': 323}, 0.5) # Guess
            _, wVar['BrB02DstDstpK0Up'], wVar['BrB02DstDstpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 311}, 0.5/5.3) #Gamma 173 pdg 2020
            _, wVar['BrB02DstDstpKst0Up'], wVar['BrB02DstDstpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
        if n == 'DstmDp':
            _, wVar['BrB02DstDpK0Up'], wVar['BrB02DstDpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 411, 'MC_DstSisterPdgId_light': 311}, 0.5/3.2) #Gamma 172 pdg 2020
            _, wVar['BrB02DstDpKst0Up'], wVar['BrB02DstDpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 411, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
            _, wVar['BrB02DstDstpK0Up'], wVar['BrB02DstDstpK0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 311}, 0.2/2.7) #Gamma 173 pdg 2020
            _, wVar['BrB02DstDstpKst0Up'], wVar['BrB02DstDstpKst0Down'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 413, 'MC_DstSisterPdgId_light': 313}, 0.5) # Guess
        if n == 'DstmDsp':
            _, wVar['BrB02DstDsUp'], wVar['BrB02DstDsDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 431}, 1.1/8.0) #Gamma 83 pdg 2020
            _, wVar['BrB02DstDsstUp'], wVar['BrB02DstDsstDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 433}, .14/1.77) #Gamma 85 pdg 2020
            _, wVar['BrB02DstDs0stUp'], wVar['BrB02DstDs0stDown'] = computeBrVarWeights(ds, {'MC_DstSisterPdgId_heavy': 10431}, .6/1.5) #Gamma 95 pdg 2020

        # Correct the amount of random tracks from PV
        weights['tkPVfrac'], wVar['tkPVfracUp'], wVar['tkPVfracDown'] = computeTksPVweights(ds, relScale=0.1)

        print 'Computing total weights'
        weightsCentral = np.ones_like(ds['q2'])
        for w in weights.values():
            weightsCentral *= w
        wVar[''] = np.ones_like(weightsCentral)

        nGenExp = sMC.effMCgen['xsec'][0] * CMS_lumi.integrated_lumi * RDoMC_normRatio
        eff = [1, 0]
        for f, df in [sMC.effMCgen['effGEN'],
                      decayBR[n],
                      sMC.effCand['effCAND'],
                      sMC.getSkimEff(cat.name+'_trkCtrl_corr'),
                     ]:
            eff[0] *= f
            eff[1] += np.square(df/f)
        eff[1] = eff[0] * np.sqrt(eff[1])
        nTotExp = nGenExp*eff[0]

        sel = {}
        scale = {}

        latexTableString = {}
        for k, selFun in sideSelecton.iteritems():
            sel[k] = selFun(ds)
            nTotSel = float(np.sum(sel[k]))
            print 'N tot selected {}: {:.0f}'.format(k, nTotSel)
            nExp = nTotExp * nTotSel / sel[k].shape[0]
            print 'N tot expected {} (before weights): {:.0f}'.format(k, nExp)
            nAux = nTotExp * np.sum(weightsCentral[sel[k]]) / sel[k].shape[0]
            print 'N tot expected {} (after weights): {:.0f}'.format(k, nAux)
            latexTableString[k] = '{:.0f} ({:.0f})'.format(nAux, nTotSel)
            if not k in totalCounting.keys():
                totalCounting[k] = [0, 0]
            totalCounting[k][0] += nAux
            totalCounting[k][1] += nTotSel
            if nTotSel ==0:
                nTotSel+=1
            scale[k] = nExp/nTotSel
        s = ' & '.join([latexTableString['AddTk_'+s+'_mHad'] for s in ['p', 'm', 'pp', 'pm', 'mm']])
        print s
        eventCountingStr[n] += ' & ' + s + '\\\\'


        for name_wVar, v_wVar in wVar.iteritems():
            h_name = n
            if not name_wVar == '':
                h_name += '__' + name_wVar
            w = weightsCentral*v_wVar

            for k in sideVar.keys():
                histo[k][h_name] = create_TH1D(
                                               ds[sideVar[k]][sel[k]],
                                               name=h_name, title=h_name,
                                               binning=binning[k],
                                               opt='underflow',
                                               weights=w[sel[k]], scale_histo=scale[k]
                                              )

    s = ' & '.join(['{:.0f} ({:.0f})'.format(*totalCounting['AddTk_'+s+'_mHad']) for s in ['p', 'm', 'pp', 'pm', 'mm']]) + ' \\\\'
    eventCountingStr['tot'] += ' & ' + s
    with open(outdir + '/eventCounting.txt', 'w') as f:
        for p in processOrder + ['tot']:
            f.write(p + '   ' + eventCountingStr[p] + '\n')


    ######################################################
    ########## Create total MC histograms and Pseudo data if necessary
    ######################################################
    print 'Creating total MC'
    for k, hDic in histo.iteritems():
        h = hDic.values()[0].Clone('total')
        h.SetTitle('Total MC')
        h.Reset()
        for n, hMC in hDic.iteritems():
            if not '__' in n and not n == 'data':
                scale = SM_RDst if 'tau' in n else 1.
                h.Add(hMC, scale)
        histo[k]['total'] = h
        if args.asimov:
            hAsimov = h.Clone('data_obs')
            hAsimov.Sumw2(0)
            for ix in range(1, hAsimov.GetNbinsX()+1):
                if hAsimov.GetNbinsY() > 1:
                    for iy in range(1, hAsimov.GetNbinsY()+1):
                        hAsimov.SetBinContent(ix, iy, np.around(hAsimov.GetBinContent(ix, iy)))
                else:
                    hAsimov.SetBinContent(ix, np.around(hAsimov.GetBinContent(ix)))
            hAsimov.Sumw2()
            histo[k]['data'] = hAsimov

    ######################################################
    ########## Create Real data Histograms
    ######################################################
    if not args.asimov:
        print 'Creating data histos'
        ds = dSet['data']
        print 'N observed data signal region: {:.1f}k'.format(1e-3*ds['q2'].shape[0])
        if args.useMVA:
            print 'Evaluating MVA'
            sMVA = computeVarMVA(ds)
            histo['MVA']['data'] = create_TH1D(sMVA, name='data_obs', binning=binning['MVA'], opt='underflow,overflow')

        for i_q2 in range(len(binning['q2'])-1):
            q2_l = binning['q2'][i_q2]
            q2_h = binning['q2'][i_q2 + 1]
            sel_q2 = np.logical_and(ds['q2'] > q2_l, ds['q2'] < q2_h)
            for var in ['M2_miss', 'Est_mu', 'mu_pt', 'Dst_pt']:
                cat_name = var+'_q2bin'+str(i_q2)
                histo[cat_name]['data'] = create_TH1D(
                                                      ds[var][sel_q2],
                                                      name='data_obs', title='Data Obs',
                                                      binning=binning[var][i_q2],
                                                      opt='underflow,overflow'
                                                     )
                if not var == 'M2_miss':
                        continue
                auxS = np.column_stack((ds['M2_miss'][sel_q2], ds['Est_mu'][sel_q2]))
                cat2D = 'h2D_q2bin'+str(i_q2)
                histo[cat2D]['data'] = create_TH2D(auxS, name='data_obs', title='Data Obs',
                                                   binning=binning_2D[i_q2])
                catU = 'Unrolled_q2bin'+str(i_q2)
                validBins = nonEmptyIdxs[i_q2]
                hUnrolled = rt.TH1D('data_obs', 'Data Obs', len(validBins), 0.5, len(validBins)+0.5)
                for i, (ix, iy) in enumerate(validBins):
                    hUnrolled.SetBinContent(i+1, histo[cat2D]['data'].GetBinContent(ix, iy))
                    hUnrolled.SetBinError(i+1, histo[cat2D]['data'].GetBinError(ix, iy))
                histo[catU]['data'] = hUnrolled

        for var in ['B_pt', 'B_eta']:
            histo[var]['data'] = create_TH1D(ds[var], name='data_obs', binning=binning[var], opt='underflow,overflow')

        ds = dSetTkSide['data']
        for k in sideVar.keys():
            histo[k]['data'] = create_TH1D(
                                           ds[sideVar[k]][sideSelecton[k](ds)],
                                           name='data_obs', title='Data Obs',
                                           binning=binning[k],
                                           opt='underflow',
                                          )
        print 'N observed data control regions: ' + ' & '.join(['{:.0f}'.format(histo['AddTk_'+s+'_mHad']['data'].Integral()) for s in ['p', 'm', 'pp', 'pm', 'mm']]) + ' \\\\'

    ######################################################
    ########## Dump root file
    ######################################################

    if not os.path.isdir(histo_file_dir): os.makedirs(histo_file_dir)

    for cat_name, h_dic in histo.iteritems():
        tf = rt.TFile(histo_file_dir+'{}_{}.root'.format(card_name, cat_name), 'recreate')
        for v in h_dic.values():
            v.Write()
        tf.Close()


########################### -------- Create histrograms ------------------ #########################

def drawPlots(tag, hDic, scale_dic={}):
    outCanvas = []
    print 20*'-', 'Drawing plots', tag, 20*'-'
    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        nameU = 'Unrolled_q2bin'+str(i_q2)
        if not nameU in hDic.keys():
            print nameU, 'not found'
            continue
        print 'Creating', nameU
        hDic[nameU]['data'].GetXaxis().SetTitle('Unrolled 2D bins')
        hDic[nameU]['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic[nameU], scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=False,
                                   addText='Cat. '+category.capitalize()+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=True, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1,
                                   tag=tag+'Unrolled_q2bin'+str(i_q2),
                                   legLoc=[0.15, 0.5, 0.25, 0.8],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True),
                                   figsize = [900, 400]
                                   )
        cAux.SaveAs(outdir+'/fig/Unrolled_q2bin'+str(i_q2)+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/Unrolled_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)

    print 'Creating signal region grid'
    cAux = plot_gridVarQ2(CMS_lumi, binning, hDic, draw_pulls=True, scale_dic=scale_dic,
                          categoryText=cat.name.capitalize(), cNameTag=tag,
                          iq2_maskData=[] if args.unblinded else [2, 3])
    cAux.SaveAs(outdir+'/fig/signalRegion_'+tag+'.png')
    cAux.SaveAs(webFolder+'/signalRegion_'+tag+'.png')
    outCanvas.append(cAux)

    if 'MVA' in hDic.keys():
        print 'Creating MVA'
        cAux = plot_SingleCategory(CMS_lumi, hDic['MVA'], draw_pulls=True, scale_dic=scale_dic,
                                   addText=cat.name.capitalize(), logy=True, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1, tag=tag+'MVA', legLoc=[0.2, 0.1, 0.4, 0.4])
        cAux.SaveAs(outdir+'/fig/MVA_'+tag+'.png')
        cAux.SaveAs(webFolder+'/MVA_'+tag+'.png')
        outCanvas.append(cAux)

    if 'B_pt' in hDic.keys():
        print 'Creating B_pt'
        hDic['B_pt']['data'].GetXaxis().SetTitle('B p_{T} [GeV]')
        hDic['B_pt']['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['B_pt'], draw_pulls=True, scale_dic=scale_dic,
                                   addText=cat.name.capitalize(), logy=False, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1, tag=tag+'B_pt', legLoc=[0.65, 0.4, 0.9, 0.75])
        cAux.SaveAs(outdir+'/fig/B_pt_'+tag+'.png')
        cAux.SaveAs(webFolder+'/B_pt_'+tag+'.png')
        outCanvas.append(cAux)

    if 'B_eta' in hDic.keys():
        print 'Creating B_eta'
        hDic['B_eta']['data'].GetXaxis().SetTitle('B #eta')
        hDic['B_eta']['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['B_eta'], draw_pulls=True, scale_dic=scale_dic,
                                   addText=cat.name.capitalize(), logy=False, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1, tag=tag+'B_eta', legLoc=[0.44, 0.23, 0.63, 0.53])
        cAux.SaveAs(outdir+'/fig/B_eta_'+tag+'.png')
        cAux.SaveAs(webFolder+'/B_eta_'+tag+'.png')
        outCanvas.append(cAux)

    for k in np.sort([k for k in hDic.keys() if 'AddTk' in k]):
        print 'Creating', k
        legLoc = [0.67, 0.3, 0.93, 0.72]
        if 'Vis' in k:
            legLoc = [0.18, 0.4, 0.4, 0.75]
        cAux = plot_SingleCategory(CMS_lumi, hDic[k], scale_dic=scale_dic,
                                   xtitle=getControlXtitle(k),
                                   addText='Cat. '+cat.name.capitalize() + ', ' + getControlSideText(k),
                                   tag=k, legLoc=legLoc,
                                   draw_pulls=True
                                   )
        cAux.SaveAs(outdir+'/fig/'+k+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/'+k+'_'+tag+'.png')
        outCanvas.append(cAux)

    # 2D plot
    for i_q2 in range(len(binning['q2'])-1):
        name2D = 'h2D_q2bin'+str(i_q2)
        if not name2D in hDic.keys():
            continue
        print 'Creating', name2D
        hSum = hDic[name2D]['total']
        hSum.GetXaxis().SetTitle('m^{2}_{miss} [GeV^{2}]')
        hSum.GetYaxis().SetTitle('E_{#mu}* [GeV]')
        hSum.GetZaxis().SetTitle('Events')
        cAux = drawOnCMSCanvas(CMS_lumi, [hSum], ['colz'], tag=str(i_q2), mR=0.17)
        cAux.SetLogz()
        cAux.SaveAs(outdir+'/fig/M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        cAux.SaveAs(webFolder+'/M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        outCanvas.append(cAux)

    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        name = 'mu_pt_q2bin'+str(i_q2)
        if not name in hDic.keys(): continue
        print 'Creating', name
        hDic[name]['data'].GetXaxis().SetTitle('#mu p_{T} [GeV]')
        hDic[name]['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic[name], scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=True,
                                   addText='Cat. '+category.capitalize()+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=False, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1,
                                   tag=tag+'mu_pt_q2bin'+str(i_q2),
                                   legLoc=[0.7, 0.5, 0.9, 0.75],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                   )
        cAux.SaveAs(outdir+'/fig/muPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/muPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)

    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        name = 'Dst_pt_q2bin'+str(i_q2)
        if not name in hDic.keys(): continue
        print 'Creating', name
        hDic[name]['data'].GetXaxis().SetTitle('D* p_{T} [GeV]')
        hDic[name]['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic[name], scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=True,
                                   addText='Cat. '+category.capitalize()+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=False, legBkg=True,
                                   procOrder = ['tau', 'DstD', 'Dstst', 'mu'],
                                   min_y=1,
                                   tag=tag+'Dst_pt_q2bin'+str(i_q2),
                                   legLoc=[0.7, 0.4, 0.9, 0.75],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                   )
        cAux.SaveAs(outdir+'/fig/DstPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/DstPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)

    print 30*'-' + '\n\n'
    return outCanvas



########################### -------- Create the card ------------------ #########################

def createCard(histo, fitRegionsOnly=False):
    print '-----> Creating the datacard' + (' (fit regions only)' if fitRegionsOnly else '')
    processes = processOrder
    nProc = len(processes)

    categories = []
    for c in np.sort(histo.keys()):
        if c.startswith('h2'): continue
        if fitRegionsOnly:
            if c == 'AddTk_pm_mHad': continue
            aux = c.startswith('Unrolled')
            aux = aux or c.startswith('AddTk_')
            if not aux: continue
            if (not args.unblinded) and (c.endswith('_q2bin2') or c.endswith('_q2bin3')): continue
        categories.append(c)
    nCat = len(categories)

    ######################################################
    ########## Define categories (bin) and processed
    ######################################################

    card = 'imax *\njmax *\nkmax *\n'
    card += 60*'-'+'\n'

    for k in categories:
        fname = histo_file_dir+'{}_{}.root'.format(card_name, k)
        card += 'shapes * {} {} $PROCESS $PROCESS__$SYSTEMATIC\n'.format(k, fname)
    card += 60*'-'+'\n'

    # number of events observed
    card += 'bin ' + ' '.join(categories) + '\n'
    obs = map(lambda k: '{:.0f}'.format(histo[k]['data'].Integral()), categories)
    obs = ' '.join(obs)
    card += 'observation ' + obs + '\n'
    card += 60*'-'+'\n'

    # MC expected events
    aux_bin = ''
    aux_proc_name = ''
    aux_proc_id = ''
    aux_proc_rate = ''
    for c, p in itertools.product(categories, processes):
        aux_bin += ' '+c
        aux_proc_name += ' '+p
        aux_proc_id += ' '+str(np.argmax(np.array(processes) == p))
        aux_proc_rate += ' {:.2f}'.format(histo[c][p].Integral())

    card += 'bin' + aux_bin + '\n'
    card += 'process' + aux_proc_name + '\n'
    # Zero or negative for sig and positive for bkg
    card += 'process' + aux_proc_id + '\n'
    # Expected rate
    card += 'rate' + aux_proc_rate + '\n'
    card += 60*'-'+'\n'


    ######################################################
    ########## Scale systematics uncertainties
    ######################################################
    #pp -> bb cros-section * luminosity
    card += 'xsecpp2bbXlumi'+cat.trg+' lnN' + ' 1.2'*nProc*nCat + '\n'

    # Branching ratio uncertainty
    decayBR = pickle.load(open('/storage/user/ocerri/BPhysics/data/forcedDecayChannelsFactors.pickle', 'rb'))

    for n in processes:
        if n in ['tau',
                 'DstPip', 'DstPi0',
                 'DstPipPim', 'DstPipPi0', 'DstPi0Pi0',
                 'TauDstPip', 'TauDstPi0'
                 ]: continue
        val = ' {:.2f}'.format(1+decayBR[n][1]/decayBR[n][0])
        aux = ''
        for nn in processes:
            if nn == n or (n=='mu' and nn=='tau'):
                aux += val
            else:
                aux += ' -'
        card += n + 'Br lnN' + aux*nCat + '\n'


    # Branching ratio uncertainty with isospin symmetry constraint
    val = ' {:.2f}'.format(1+decayBR['DstPip'][1]/decayBR['DstPip'][0]) #DstPi0 has no info
    aux = ''
    for n in processes:
        if n in ['DstPip', 'DstPi0']: aux += val
        else: aux += ' -'
    card += 'DstPiBr lnN' + aux*nCat + '\n'

    val = ' {:.2f}'.format(1+decayBR['DstPipPim'][1]/decayBR['DstPipPim'][0]) #DstPipPim is the only one meaasured
    aux = ''
    for n in processes:
        if n in ['DstPipPim', 'DstPipPi0', 'DstPi0Pi0']: aux += val
        else: aux += ' -'
    card += 'DstPiPiBr lnN' + aux*nCat + '\n'

    val = ' {:.2f}'.format(1+decayBR['TauDstPip'][1]/decayBR['TauDstPip'][0]) #DstPi0 has no info
    aux = ''
    for n in processes:
        if n in ['TauDstPip', 'TauDstPi0']: aux += val
        else: aux += ' -'
    card += 'TauDstPiBr lnN' + aux*nCat + '\n'


    card += 60*'-'+'\n'


    ######################################################
    ########## Shape systematics uncertainties
    ######################################################

    # card += 'trg{}SF shape'.format(cat.trg) + ' 1.'*nProc*nCat + '\n'
    # card += 'muonIdSF shape' + ' 1.'*nProc*nCat + '\n'

    aux = ''
    for c in categories:
        if c.startswith('AddTk_'):
            aux += ' 1.'*nProc
        else: aux += ' -'*nProc
    card += 'tkPVfrac shape' + aux + '\n'

    # Mu pT spectrum
    # aux = ''
    # for p in processes:
    #     aux += ' 1.'
    # names = []
    # for k in histo.values()[0].keys():
    #     if k.startswith('mu__MuPt') and k.endswith('Up'):
    #         names.append(k[4:-2])
    # for n in sorted(names):
    #     card += n + ' shape' + aux*nCat + '\n'

    # B0 pT spectrum
    aux = ''
    for p in processes:
        if p in SamplesB0:
            aux += ' 0.7'
        else:
            aux += ' -'
    names = []
    for k in histo.values()[0].keys():
        if k.startswith(SamplesB0[0]+'__B0pT') and k.endswith('Up'):
            names.append(k[4:-2])
    if len(names) == 1 and names[0]=='B0pT':
        card += 'B0pT shape' + aux*nCat + '\n'
    else:
        for n in sorted(names):
            card += n + ' shape' + aux*nCat + '\n'

    # B +/- pT spectrum
    aux = ''
    for p in processes:
        if p in SamplesBp:
            aux += ' 1.'
        else:
            aux += ' -'
    card += 'BpPt shape' + aux*nCat + '\n'


    # Form Factors from Hammer
    for n_pFF in FreeParFF:
        aux = ''
        for p in processes:
            if p in ['tau', 'mu']:
                aux += ' 1.'
            else:
                aux += ' -'
        card += 'B2Dst'+schemeFF+'{} shape'.format(n_pFF) + aux*nCat + '\n'


    # Dstst mix composition
    aux = ''
    for p in processes:
        if p == 'DstPip' or p == 'TauDstPip':
            aux += ' 1.'
        else:
            aux += ' -'
    card += 'fDststWide shape' + aux*nCat + '\n'
    card += 'D2420_10Width shape' + aux*nCat + '\n'
    card += 'D2460_1StWidth shape' + aux*nCat + '\n'

    # Dstst->DstPiPi width
    aux = ''
    for p in processes:
        if p == 'DstPipPim' or p == 'DstPi0Pi0':
            aux += ' 1.'
        else: aux += ' -'
    card += 'DstPiPiWidth shape' + aux*nCat + '\n'


    # Hc mix composition
    aux = ''
    for p in processes:
        if p == 'DstmD0': aux += ' 1.'
        else: aux += ' -'
    card += 'BrB02DstD0Kp shape' + aux*nCat + '\n'
    card += 'BrB02DstD0Kstp shape' + aux*nCat + '\n'
    card += 'BrB02DstDst0Kp shape' + aux*nCat + '\n'
    card += 'BrB02DstDst0Kstp shape' + aux*nCat + '\n'

    aux = ''
    for p in processes:
        if p == 'DstmDp': aux += ' 1.'
        else: aux += ' -'
    card += 'BrB02DstDpK0 shape' + aux*nCat + '\n'
    card += 'BrB02DstDpKst0 shape' + aux*nCat + '\n'

    aux = ''
    for p in processes:
        if p == 'DstmDp' or p == 'DstmD0': aux += ' 1.'
        else: aux += ' -'
    card += 'BrB02DstDstpK0 shape' + aux*nCat + '\n'
    card += 'BrB02DstDstpKst0 shape' + aux*nCat + '\n'

    aux = ''
    for p in processes:
        if p == 'DstmDsp': aux += ' 1.'
        else: aux += ' -'
    card += 'BrB02DstDs shape' + aux*nCat + '\n'
    card += 'BrB02DstDsst shape' + aux*nCat + '\n'
    card += 'BrB02DstDs0st shape' + aux*nCat + '\n'

    card += 60*'-'+'\n'

    ######################################################
    ########## MC statistical uncertainties
    ######################################################

    if not args.noMCstats:
        card += 'AddTk_p_mHad autoMCStats 2 1 1\n'
        card += 'AddTk_m_mHad autoMCStats 2 1 1\n'
        # card += 'AddTk_pm_mHad autoMCStats 2 1 1\n'
        card += 'AddTk_pm_mVis autoMCStats 2 1 1\n'
        card += 'AddTk_pp_mHad autoMCStats 2 1 1\n'
        card += 'AddTk_mm_mHad autoMCStats 2 1 1\n'

        if args.useMVA:
            card += 'MVA autoMCStats 2 1 1\n'
        else:
        #     card += 'Est_mu* autoMCStats 0 1 1\n'
        #     card += 'M2_miss* autoMCStats 0 1 1\n'

            card += 'Unrolled_q2bin0 autoMCStats 0 1 1\n'
            card += 'Unrolled_q2bin1 autoMCStats 0 1 1\n'
            if args.unblinded:
                card += 'Unrolled_q2bin2 autoMCStats 0 1 1\n'
                card += 'Unrolled_q2bin3 autoMCStats 0 1 1\n'

        card += 60*'-'+'\n'

    ######################################################
    ########## Defining groups of systematics
    ######################################################

    # autoMCStats group = defined by default when using autoMCStats

    if len(FreeParFF):
        aux_FF = ' '.join(['B2Dst'+schemeFF+n for n in FreeParFF])
        card += 'B2DstFF group = ' + aux_FF + '\n'

    # cardParts = card.split(60*'-'+'\n')
    # scaleNuis = []
    # for ln in cardParts[4].split('\n'):
    #     scaleNuis.append(ln.split(' ')[0])
    # shapeNuis = []
    # for ln in cardParts[5].split('\n'):
    #     if 'autoMCStats' in ln:
    #         continue
    #     shapeNuis.append(ln.split(' ')[0])
    # card += 'allSys group = ' + ' '.join(scaleNuis+shapeNuis) + '\n'

    return card


########################### -------- Create the workspace ------------------ #########################

def createWorkspace(cardLoc):
    print '-----> Creating workspace'
    print cardLoc
    cmd = 'text2workspace.py ' + cardLoc
    cmd += ' -o ' + cardLoc.replace('.txt', '.root')
    cmd += ' --no-b-only --verbose 1 --channel-masks'
    # cmd += ' --no-wrappers'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    else:
        text_file = open(cardLoc.replace('.txt', '_text2workspace.out'), "w")
        text_file.write(output)
        text_file.close()


########################### -------- Bias studies ------------------ #########################

def biasToysScan(card, out, seed=1, nToys=10, rVal=SM_RDst, maskStr=''):
    # To be debugged
    print card, out, seed, nToys
    exit()
    cmd += 'cd ' + out + '; '
    cmd += 'combine -M GenerateOnly'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --seed ' + str(seed)
    cmd += ' --noMCbonly 1'
    cmd += ' --setParameters r={} --freezeParameters r'.format(rVal)
    cmd += ' --toysFrequentist -t {} --saveToys'.format(nToys)
    cmd += ' -n Toys -m {:.0f}'.format(100*rVal)
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise

    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=50'
    cmd += ' --robustFit 1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic'
    cmd += ' --seed ' + str(seed)
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --toysFile higgsCombineToys.GenerateOnly.mH{}.{}.root -t {}'.format(100*rVal, seed, nToys)
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' --setParameterRanges r=0.15,0.5'
    cmd += ' -n Scan -m {:.0f}'.format(100*rVal)
    cmd += ' --verbose 1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    for line in output.split('\n'):
        if 'ERROR' in line: print line.replace('ERROR', '\033[1m\x1b[31mError\x1b[0m')
    if status:
        print output
        raise

def collectBiasToysResults(scansLoc, rVal=SM_RDst):
    print '-----> Collectiong bias toys scans'
    if not scansLoc[-1] == '/': scansLoc += '/'
    fnames = glob(scansLoc + 'higgsCombineScan.MultiDimFit.mH{:.0f}.*.root'.format(100*rVal))
    res = np.zeros(0)
    for fname in fnames:
        res = np.concatenate((res, getUncertaintyFromLimitTree(fname, verbose=False)), axis=0)
    r = res[:,0]
    rLoErr = res[:,1]
    rHiErr = res[:,2]

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(8,6))
    plt.errorbar(np.arange(1, 1+r.shape[0]), r, yerr=np.column_stack((rLoErr, rHiErr)).T, fmt='o', color='#1f77b4', label='Toys fit results')
    m = np.mean(r)
    sm = np.std(r)/np.sqrt(r.shape[0])
    x = [0, r.shape[0]+1]
    plt.fill_between(x, 2*[m-sm], 2*[m+sm], color='#ff7f0e', alpha=0.3)
    plt.plot(x, 2*[m], color='#d62728', lw=1, label='Toys mean')
    plt.plot(x, [rVal, rVal], 'm--', lw=2, label='Injected value')
    plt.legend(loc='upper right', numpoints=1)
    plt.xlabel('Toy number')
    plt.ylabel(r'$R(D^*)$')
    plt.savefig(outdir + '/fig/biasStudy_toysResults.png')
    plt.savefig(webFolder + '/biasStudy_toysResults.png')

    z = (r - rVal)/(0.5*(rLoErr + rHiErr))
    h = create_TH1D(z, name='hZtest', binning=[int(2*np.sqrt(r.shape[0])), -4, 4], axis_title=['#hat{R(D*)} - R(D*) / #sigma', 'Number of toys'])
    h.Sumw2()
    h.Fit('gaus', 'ILQ')
    rt.gStyle.SetStatY(0.95)
    c = drawOnCMSCanvas(CMS_lumi, [h])
    c.SaveAs(outdir + '/fig/biasStudy_zTest.png')
    c.SaveAs(webFolder + '/biasStudy_zTest.png')


########################### -------- Likelihood scan ------------------ #########################

def runScan(tag, card, out, rVal=SM_RDst, rLimits=[0.1, 0.7], nPoints=50, maskStr='', strategy=1, draw=True):
    if not out[-1] == '/': out += '/'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit'
    cmd += ' --algo grid --points='+str(nPoints)
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy='+str(strategy)
    cmd += ' --X-rtd MINIMIZER_analytic'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --rMin={:.2f} --rMax={:.2f}'.format(*rLimits)
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' -n ' + tag
    cmd += ' --verbose -1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise

    if draw:
        json.dump({'r': 'R(D*)'}, open(out+'renameDicLikelihoodScan.json', 'w'))

        cmd = 'cd ' + out + '; '
        cmd += 'plot1DScan.py higgsCombine{t}.MultiDimFit.mH120.root -o scan{t}'.format(t=tag)
        cmd += ' --main-label "{} {}'.format('Obs.' if not args.asimov else 'Asimov', category.capitalize())
        if not args.unblinded: cmd += ' (blinded)'
        cmd += '"'
        cmd += ' --translate ' + out+'renameDicLikelihoodScan.json'
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise
        cmd = 'cp {}/scan{}.png {}/'.format(out, tag, webFolder)
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise

    print 'Extracting new POI boundaries'
    res = getUncertaintyFromLimitTree(out+'higgsCombine{}.MultiDimFit.mH120.root'.format(tag))
    rLimitsOut = [res[0] - 3*res[1], res[0] + 3*res[2]]
    rValOur = res[0]
    if args.showPlots:
        display(Image(filename=out+'/scan'+tag+'.png'))
    return rValOur, rLimitsOut



########################### -------- Fit Diagnostic ------------------ #########################

def runFitDiagnostic(tag, card, out, forceRDst=False, maskStr='', rVal=SM_RDst, rLimits=[0.1, 0.7], seed=6741):
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M FitDiagnostics'
    cmd += ' --robustFit 1 --robustHesse 1 --cminDefaultMinimizerStrategy 2 --X-rtd MINIMIZER_analytic'
    cmd += ' --seed ' + str(seed)
    cmd += ' -d ' + card.replace('.txt', '.root')
    if forceRDst:
        cmd += ' --setParameterRanges r={:.2f},{:.2f}'.format(0, 1)
        cmd += ' --customStartingPoint --setParameters r={:.3f}'.format(rVal)
    else:
        cmd += ' --skipBOnlyFit'
        cmd += ' --setParameterRanges r={:.2f},{:.2f}'.format(*rLimits)
        cmd += ' --setParameters r={:.3f}'.format(rVal)
    if maskStr:
        cmd += ',' + maskStr
    runName = tag + ('_RDstFixed' if forceRDst else '')
    cmd += ' -n ' + runName
    cmd += ' --saveShapes --saveWithUncertainties --saveNormalizations'
    cmd += ' --trackParameters rgx{.*}'
    cmd += ' --plots'
    cmd += ' --verbose 0'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status or 'There was a crash.' in output:
        print output
        raise
    for line in output.split('\n'):
            if 'ERROR' in line: print line.replace('ERROR', '\033[1m\x1b[31mError\x1b[0m')
            if line.startswith('customStartingPoint'): print line

    if not out[-1] == '/': out += '/'
    arr = rtnp.root2array(out + 'higgsCombine{}.FitDiagnostics.mH120.{}.root'.format(runName, seed), treename='limit')
    c, d, u, _ = arr['limit']
    print 'R(D*) = {:.3f} +{:.3f}/-{:.3f} [{:.1f} %]'.format(c, u-c, c-d, 100*(u-d)*0.5/c)

def getPostfitHistos(tag, out, forceRDst, histo_prefit):
    runName = tag + ('_RDstFixed' if forceRDst else '')

    # Get post-fit shapes
    if not out[-1] == '/': out += '/'
    n = out + 'fitDiagnostics{}.root'.format(runName)
    print 'Loading post-fit from:', n
    fFitDiagnostics = rt.TFile(n, 'READ')
    if forceRDst:
        fd = fFitDiagnostics.shapes_fit_b
    else:
        fd = fFitDiagnostics.shapes_fit_s

    histo_postfit = {}
    for regName in [k.GetTitle() for k in fd.GetListOfKeys()]:
        histo_postfit[regName] = {}

        for n, h in histo_prefit[regName].iteritems():
            if '__' in n or n =='total':
                continue
            h_post = h.Clone(h.GetName() + '_postfit')
            if 'data' in n:
                h_fit = fd.Get(regName+'/total')
                h_data = h.Clone(h.GetName() + '_data')
                for i in range(1, h_post.GetNbinsX()+1):
                    h_post.SetBinContent(i, h_fit.GetBinContent(i))
                    h_post.SetBinError(i, h_fit.GetBinError(i))

                histo_postfit[regName]['total'] = h_post
                histo_postfit[regName]['data'] = h_data
            else:
                h_fit = fd.Get(regName+'/'+n)
                if not h_fit:
                    print n+' missing from '+regName
                    continue
                for i in range(1, h_post.GetNbinsX()+1):
                    h_post.SetBinContent(i, h_fit.GetBinContent(i))
                    h_post.SetBinError(i, h_fit.GetBinError(i))

                histo_postfit[regName][n] = h_post

    h2 = fFitDiagnostics.Get('covariance_fit_' + ('b' if forceRDst else 's'))
    rt.gStyle.SetPaintTextFormat('.1f')

    N = h2.GetNbinsX()
    n=35
    h2.GetXaxis().SetRange(1, n)
    h2.GetYaxis().SetRangeUser(N-n, N)
    h2.SetMarkerSize(.8)
    h2.LabelsOption("v")
    h2.GetXaxis().SetLabelSize(0.03)
    h2.GetYaxis().SetLabelSize(0.03)
    CC = drawOnCMSCanvas(CMS_lumi, [h2, h2], ['colz', 'text same'], size=(900, 700), tag='tl', mL=0.22, mR=0.15, mB=0.25)
    CC.SaveAs(out+'fig/covariance_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')
    CC.SaveAs(webFolder+'/covariance_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')

    return histo_postfit, CC, fFitDiagnostics

def nuisancesDiff(tag, out, forceRDst):
    runName = tag + ('_RDstFixed' if forceRDst else '')
    cmd = 'python diffNuisances.py ' + out + '/fitDiagnostics{}.root'.format(runName)
    if not forceRDst:
        cmd += ' --skipFitB'
    cmd += ' --all'
    cmd += ' --abs'
    cmd += ' -g {}/nuisance_difference'.format(out) + runName + '.root'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    dumpDiffNuisances(output, out, tag='_RDstFixed' if forceRDst else '',
                      useBonlyResults=forceRDst, parsToPrint=100)
    cmd = 'cp {}nuisance_difference.txt {}/'.format(out, webFolder)
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise


########################### -------- Uncertainty breakdown ------------------ #########################

def runUncertaintyBreakDown(card, out, rVal=SM_RDst, rLimits=[0.1, 0.7], maskStr=''):
    if not out[-1] == '/': out += '/'
    print '--------> Running uncertainty breakdown <--------------'
    print '--------> Nominal scan'
    rValOut, rLimitsOut = runScan('Nominal', card, out, rVal, rLimits, nPoints=80, maskStr=maskStr, strategy=1, draw=False)
    sig = (rLimitsOut[1] - rValOut)/3
    rLimitsTight = [rValOut - 2*sig, rValOut + 2*sig]

    print '--------> Best fit snap'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --setParameters r={:.2f}'.format(rValOut)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' --setParameterRanges r={:.3f},{:.3f}'.format(*rLimitsTight)
    cmd += ' -n Bestfit'
    cmd += ' --saveWorkspace --verbose -1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise

    print '--------> Statistical uncertanty only'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=100'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' -d higgsCombineBestfit.MultiDimFit.mH120.root'
    cmd += ' --snapshotName MultiDimFit'
    cmd += ' --rMin={:.3f} --rMax={:.3f}'.format(*rLimitsTight)
    cmd += ' -n StatOnly'
    cmd += ' --freezeParameters allConstrainedNuisances'
    if maskStr: cmd += ' --setParameters ' + maskStr
    cmd += ' --fastScan' # To be added if there are no free parameters otherwise
    cmd += ' --verbose -1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    getUncertaintyFromLimitTree(out + 'higgsCombineStatOnly.MultiDimFit.mH120.root')

    print '--------> MC stats and Statistical uncertanty only'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=100'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' -d higgsCombineBestfit.MultiDimFit.mH120.root'
    cmd += ' --rMin={:.4f} --rMax={:.4f}'.format(*rLimitsTight)
    cmd += ' -n MCstat'
    cmd += ' --snapshotName MultiDimFit'
    if maskStr: cmd += ' --setParameters ' + maskStr
    cmd += ' --freezeNuisanceGroups=autoMCStats'
    cmd += ' --verbose -1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    getUncertaintyFromLimitTree(out + 'higgsCombineMCstat.MultiDimFit.mH120.root')


    json.dump({'r': 'R(D*)'}, open(out+'renameDicLikelihoodScan.json', 'w'))

    cmd = 'cd ' + out + '; '
    cmd += 'plot1DScan.py higgsCombineNominal.MultiDimFit.mH120.root'
    cmd += ' --main-label "{} {}'.format('Obs.' if not args.asimov else 'Asimov', category.capitalize())
    if not args.unblinded: cmd += ' (blinded)'
    cmd += '"'
    cmd += ' --others'
    cmd += ' "higgsCombineMCstat.MultiDimFit.mH120.root:Stat. + Syst.:4"'
    cmd += ' "higgsCombineStatOnly.MultiDimFit.mH120.root:Stat. only:2"'
    cmd += ' --breakdown "MC stat.,syst.,stat."'
    cmd += ' --translate ' + out+'renameDicLikelihoodScan.json'
    cmd += ' -o scanBreakdown'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    cmd = 'cp {}scanBreakdown.png {}/'.format(out, webFolder)
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    if args.showPlots:
        display(Image(filename=out+'scanBreakdown.png'))



########################### -------- Nuisances impact ------------------ #########################

def runNuisanceImpacts(card, out, maskStr='', rVal=SM_RDst, submit=True, collect=True):
    if not out[-1] == '/': out += '/'

    if submit:
        if os.path.isdir(out+'impactPlots'):
            os.system('rm -rf '+out+'impactPlots')
        os.mkdir(out+'impactPlots')

        cmd = 'cd {}impactPlots; '.format(out)
        cmd += ' combineTool.py -M Impacts --doInitialFit -m 0'
        cmd += ' --robustFit 1 --X-rtd MINIMIZER_analytic'
        cmd += ' -d ' + card.replace('.txt', '.root')
        cmd += ' --setParameters r={:.2f}'.format(rVal)
        if maskStr:
            cmd += ','+maskStr
        cmd += ' --setParameterRanges r=0.1,0.6'
        cmd += ' --verbose -1'
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise


        # If running on Tier2 condor remmeber to add this line to CombineToolBase.py ln 11
        # ``source /cvmfs/cms.cern.ch/cmsset_default.sh``
        cmd = 'cd {}impactPlots;'.format(out)
        cmd += ' combineTool.py -M Impacts --doFits -m 0'
        cmd += ' --robustFit 1 --X-rtd MINIMIZER_analytic'
        cmd += ' --parallel 100 --job-mode condor --task-name combineImpacts_'+category
        cmd += ' --sub-opts "{}"'.format(stringJubCustomizationCaltechT2.replace('"', '\\\"').replace('$', '\$'))
        cmd += ' -d ' + card.replace('.txt', '.root')
        cmd += ' --setParameters r={:.2f}'.format(rVal)
        if maskStr:
            cmd += ','+maskStr
        cmd += ' --verbose -1'
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise

    if collect:
        status, output = commands.getstatusoutput('condor_q')
        while 'combineImpacts_'+category in output:
            time.sleep(30)
            status, output = commands.getstatusoutput('condor_q')
            for l in output.split('\n'):
                if 'combineImpacts_'+category in l:
                    print l
                    sys.stdout.flush()
        cmd = 'cd {}impactPlots;'.format(out)
        cmd += ' combineTool.py -M Impacts -o impacts.json -m 0'
        cmd += ' -d ' + card.replace('.txt', '.root')
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise


        rename = {
        'r': 'R(D*)',
        'B0pT': 'B_{0} p_{T} spectrum',
        'B2DstCLNR0':'R_{0} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig1':'#lambda_{1} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig2':'#lambda_{2} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig3':'#lambda_{3} (CLN B#rightarrow D*l#nu)',
        'trgSF': 'Trigger scale factor',
        'xsecpp2bbXlumi': 'Luminosity*#sigma_{pp#rightarrowbb}',
        }
        for i in range(1,5):
            s = str(i)
            rename['B0pT_lam'+s] = 'B_{0} p_{T} #lambda_{'+s+'}'

        procName_dic = {
        'mu'        : 'B_{0}#rightarrow D*#mu#nu',
        'tau'       : 'B_{0}#rightarrow D*#tau#nu',
        'DstmD0'    : 'B^{+}#rightarrow D*D_{0}(#muY) + X',
        'DstmDp'    : 'B^{+}#rightarrow D*D^{+}(#muY) + X',
        'DstmDsp'   : 'B^{+}#rightarrow D*D_{s}^{+}(#muX)',
        'DstPip'    : 'B^{+}#rightarrow D*#pi^{+}#mu#nu',
        'DstPipPi0' : 'B^{+}#rightarrow D*#pi^{+}#pi^{0}#mu#nu',
        'DstPi0'    : 'B_{0}#rightarrow D*#pi^{0}#mu#nu',
        'DstPipPim' : 'B_{0}#rightarrow D*#pi^{+}#pi^{-}#mu#nu',
        'DstPi0Pi0' : 'B_{0}#rightarrow D*#pi^{0}#pi^{0}#mu#nu',
        'BpDstmHc'  : 'B^{+}#rightarrow D*D(#muX)',
        'BmDstmHc'  : 'B^{-}#rightarrow D*D(#muX)',
        'antiB0DstmHc'  : '#bar{B}_{0}#rightarrow D*D(#muX)',
        'DstPi'     : 'B #rightarrow D**(#rightarrow D*#pi)#mu#nu',
        'DstPiPi'   : 'B #rightarrow D**(#rightarrow D*#pi#pi)#mu#nu',
        }

        for n in procName_dic:
            rename[n+'Br'] = 'Branching fraction ' + procName_dic[n]

        d = json.load(open(out+'impactPlots/impacts.json', 'r'))
        for par in d['params']:
            name = str(par['name'])
            if not name.startswith('prop_bin'): continue
            label = name.replace('prop_bin', 'MC stat. ')
            label = label.replace('M2_miss_', 'M^{2}_{miss} ')
            label = label.replace('Est_mu_', 'E*_{#mu} ')
            label = label.replace('q2bin', '[b_{q^{2}}=')
            label = label.replace('_bin', '] ')
            rename[name] = label + 10*' '
        json.dump(rename, open(out+'impactPlots/rename.json', 'w'))

        cmd = 'cd {};'.format(out)
        cmd += 'plotImpacts.py -i impactPlots/impacts.json -o impacts -t impactPlots/rename.json --max-pages 1'
        print cmd
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise

        cmd = 'cp {}impacts.pdf {}/'.format(out, webFolder)
        status, output = commands.getstatusoutput(cmd)
        if status:
            print output
            raise


########################### -------- Goodness of Fit ------------------ #########################
def runCommand(cmd):
    status, output = commands.getstatusoutput(cmd)
    return [status, output]

def runGoodnessOfFit(tag, card, out, algo, maskEvalGoF='', fixRDst=False, rVal=SM_RDst):
    # Always to be tun with fitRegionsOnly cards
    tag += '_algo'+algo
    if fixRDst:
        tag += '_fixRDst'

    if not out[-1] == '/': out += '/'
    gofOutdir = out + 'goodnessOfFit'+tag
    if os.path.isdir(gofOutdir): os.system('rm -rf ' + gofOutdir)
    os.system('mkdir ' + gofOutdir)


    cmd = 'cd ' + gofOutdir + '; '
    cmd += 'combine -M GoodnessOfFit'
    cmd += ' --algo=saturated  --toysFrequentist' if algo=='Sat' else ' --algo='+algo
    cmd += ' --X-rtd MINIMIZER_analytic  --cminDefaultMinimizerStrategy=1'
    cmd += ' -d ' + card.replace('.txt', '.root')
    if fixRDst:
        cmd += ' --freezeParameters r --setParameters r={:.3f}'.format(rVal)
    if maskEvalGoF:
        cmd += ' --setParametersForEval ' + maskEvalGoF
    cmd += ' -n Obs'+tag
    cmd += ' -t 0 -s 100'
    cmd += ' --verbose -1'
    print cmd
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    arr = rtnp.root2array(gofOutdir+'/higgsCombineObs'+tag+'.GoodnessOfFit.mH120.100.root', treename='limit')
    s_obs = arr['limit'][0]
    print 'Observed test statistics: {:.1f}'.format(s_obs)


    # Run the test stat toy distribution
    cmdToys = cmd.replace('-n Obs', '-n Toys')
    cmdToys = cmdToys.replace('-t 0 -s 100', '-t 20 -s -1')
    print cmdToys

    Nrep = 20
    p = Pool(min(20,Nrep))
    outputs = p.map(runCommand, Nrep*[cmdToys])
    for s,o in outputs:
        if s: print o


    # Get the p-value
    s_toys = []
    for name_toys in glob(gofOutdir+'/higgsCombineToys'+tag+'.GoodnessOfFit.*.root'):
        arr = rtnp.root2array(name_toys, treename='limit')
        s_toys += list(arr['limit'])
    p_val = np.sum(s_toys > s_obs)/float(len(s_toys))

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(8,6))
    content, center, _ = plt.hist(s_toys, weights=np.ones_like(s_toys)/float(len(s_toys)),
                                  alpha=0.7, label='Toys ({:.0f})'.format(float(len(s_toys))))
    plt.plot([s_obs, s_obs], [0, np.max(content)], 'm--', label='Observed\np-val {:.1f}%'.format(100*p_val))
    plt.legend(loc='upper right')
    plt.xlabel('Test statistic (algo {})'.format(algo))
    plt.ylabel('Probability / {:.1f}'.format(0.5*(center[2]-center[1])))
    plt.savefig(out + 'fig/resultsGoF'+tag+'.png')
    plt.savefig(webFolder + '/resultsGoF'+tag+'.png')

    strRes = tag
    strRes += ' '*(55-len(strRes))
    strRes += '{:.2f}'.format(s_obs)
    strRes += ' '*(70-len(strRes))
    strRes += '{:.2f} {:.2f} {:.2f}'.format(np.percentile(s_toys, 50), np.percentile(s_toys, 95), np.percentile(s_toys, 99))
    print strRes
    os.system('echo "{}" >> {}GoF_results.txt'.format(strRes, out));




if __name__ == "__main__":

    if 'histos' in args.step:
        createHistograms()

    histo = loadHisto4CombineFromRoot(histo_file_dir, card_name)

    if 'preFitPlots' in args.step:
        cPre = drawPlots('prefit', histo, scale_dic={'tau': SM_RDst})

    if 'card' in args.step:
        for fitRegionsOnly in [False, True]:
            card = createCard(histo, fitRegionsOnly)
            if args.showCard:
                print 3*'\n'
                print '--------> Dumping Combine card (fitRegionsOnly {})'.format(fitRegionsOnly)
                print card
                print 3*'\n'
            cl = card_location.replace('.txt', '_fitRegionsOnly.txt') if fitRegionsOnly else card_location
            print 'Card location:', cl
            fc = open(cl, 'w')
            fc.write(card)
            fc.close()
        print '\n'

    if 'workspace' in args.step:
        createWorkspace(card_location)
        createWorkspace(card_location.replace('.txt', '_fitRegionsOnly.txt'))
        print '\n'

    if 'bias' in args.step:
        biasOut = outdir + '/biasStudyToys'
        biasToysScan(card_location.replace('.txt', '_fitRegionsOnly.txt'), biasOut, args.seed)
        collectBiasToysResults(biasOut)

    rDst_postFitRegion = args.RDstLims
    fit_RDst = SM_RDst

    globalChannelMasking = ['AddTk_pm_mHad', 'B_pt', 'B_eta']
    for n in ['Est_mu', 'M2_miss', 'mu_pt', 'Dst_pt']:
        for i in range(len(binning['q2'])-1):
            globalChannelMasking.append(n+'_q2bin'+str(i))

    if not args.unblinded:
        globalChannelMasking += ['Unrolled_q2bin2', 'Unrolled_q2bin3']
    globalChannelMaskingStr = ','.join(['mask_{}=1'.format(c) for c in globalChannelMasking])

    if 'scan' in args.step:
        if args.maskScan:
            maskList = []
            for n in args.maskScan:
                for kn in histo.keys():
                    if not re.match(n, kn) is None:
                        print 'Masking', kn
                        maskList.append('mask_'+kn+'=1')
            maskStr = ','.join(maskList)
            fit_RDst, rDst_postFitRegion = runScan(args.cardTag+args.tagScan, card_location, maskStr=maskStr,
                                                   out=outdir, rLimits=rDst_postFitRegion, strategy=1, draw=True)
        else:
            fit_RDst, rDst_postFitRegion = runScan(args.cardTag+'Base', card_location.replace('.txt', '_fitRegionsOnly.txt'),
                                                   out=outdir, rLimits=rDst_postFitRegion, strategy=1, draw=True)

    if 'fitDiag' in args.step:
        print '-----> Running fit diagnostic'
        runFitDiagnostic(args.cardTag, card_location, outdir,
                         forceRDst=args.forceRDst, maskStr=globalChannelMaskingStr,
                         rVal=fit_RDst, rLimits=rDst_postFitRegion)
    if 'postFitPlots' in args.step:
        print '-----> Getting postfit results'
        histo_post, _, _ = getPostfitHistos(args.cardTag, outdir, forceRDst=args.forceRDst, histo_prefit=histo)
        tag = 'postfit'+ ('_RDstFixed' if args.forceRDst else '')
        cPost = drawPlots(tag, histo_post)
        nuisancesDiff(args.cardTag, outdir, args.forceRDst)
        print '\n'

    if 'uncBreakdown' in args.step:
        runUncertaintyBreakDown(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir, rVal=fit_RDst, rLimits=rDst_postFitRegion)

    if 'impacts' in args.step:
        print '-----> Running impact plots'
        submit, collect = True, True
        if args.collectImpacts and not args.subOnlyImpacts:
            submit, collect = False, True
        elif not args.collectImpacts and args.subOnlyImpacts:
            submit, collect = True, False
        runNuisanceImpacts(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir, rVal=fit_RDst, submit=submit, collect=collect)

    if 'GoF' in args.step:
        print '-----> Goodnees of Fit'
        maskList = []
        if args.maskEvalGoF:
            if len(args.algoGoF) > 1 or args.algoGoF[0] != 'Sat':
                print 'Only saturated algorith accept masks. Running only algo=Sat'
                args.algoGoF = ['Sat']

            for n in args.maskEvalGoF:
                for kn in histo.keys():
                    if not re.match(n, kn) is None:
                        print 'Masking', kn
                        maskList.append('mask_'+kn+'=1')
        maskStr = ','.join(maskList)
        for algo in args.algoGoF:
            runGoodnessOfFit(args.tagGoF, card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                             algo, fixRDst=args.forceRDst, maskEvalGoF=maskStr)
            print '-'
