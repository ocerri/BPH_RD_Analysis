#!/usr/bin/env python

"""
######### Notes #########
To activate the environment run: cd ~/work/CMSSW_10_2_13/src/; cmsenv; cd ~/BPH_RD_Analysis/Combine/

To run more categories at once:
for cat in low mid high comb; do ./runCombine.py --cat $cat --submit; done
"""


import sys, os, pickle, time, json, yaml, itertools, commands, re
from datetime import datetime
from glob import glob
sys.path.append('../lib')
sys.path.append('../analysis')
from collections import defaultdict
from multiprocessing import Pool
from prettytable import PrettyTable
import numpy as np
import pandas as pd
from scipy.stats import chi2 as scipy_chi2
from scipy.stats import crystalball
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
from beamSpot_calibration import getBeamSpotCorrectionWeights
from pT_calibration_reader import pTCalReader as kinCalReader
from histo_utilities import create_TH1D, create_TH2D, std_color_list, make_ratio_plot
from gridVarQ2Plot import plot_gridVarQ2, plot_SingleCategory, getControlXtitle, getControlSideText
from lumi_utilities import getLumiByTrigger
from combine_utilities import getUncertaintyFromLimitTree, dumpDiffNuisances, stringJubCustomizationCaltechT2, loadHisto4CombineFromRoot, getResultsFromMultiDimFitSingles

import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "     Preliminary"
donotdelete = []

import argparse
parser = argparse.ArgumentParser(description='Script used to run combine on the R(D*) analysis.',
                                 epilog='Test example: ./runCombine.py -c low',
                                 add_help=True
                                 )
parser.add_argument ('--HELP', '-H', default=False, action='store_true', help='Print help message.')
parser.add_argument ('--cardTag', '-v', default='test_', help='Card name initial tag.')
parser.add_argument ('--category', '-c', type=str, default='high', choices=['single', 'low', 'mid', 'high', 'comb'], help='Category.')

# parser.add_argument ('--skimmedTag', default='_220117_lostInnerHits_fullM2miss', type=str, help='Tag to append to the skimmed directory.')
parser.add_argument ('--skimmedTag', default='', type=str, help='Tag to append to the skimmed directory.')
parser.add_argument ('--bareMC', default=True, type=bool, help='Use bare MC instead of the corrected one.')
parser.add_argument ('--maxEventsToLoad', default=None, type=int, help='Max number of MC events to load per sample.')
parser.add_argument ('--calBpT', default='poly', choices=['poly', 'none'], help='Form factor scheme to use.')
parser.add_argument ('--schemeFF', default='CLN', choices=['CLN', 'BLPR', 'NoFF'], help='Form factor scheme to use.')
parser.add_argument ('--lumiMult', default=1., type=float, help='Luminosity multiplier for asimov dataset. Only works when asimov=True')

parser.add_argument ('--useMVA', default=False, choices=[False, 'v0', 'v1'], help='Use MVA in the fit.')
parser.add_argument ('--signalRegProj1D', default='', choices=['M2_miss', 'Est_mu', 'U_miss'], help='Use 1D projections in signal region instead of the unrolled histograms')
parser.add_argument ('--unblinded', default=False, type=bool, help='Unblind the fit regions.')
parser.add_argument ('--noLowq2', default=False, action='store_true', help='Mask the low q2 signal regions.')
parser.add_argument ('--controlRegions', default=['p__mHad', 'm__mHad', 'pm_mHad', 'pp_mHad', 'mm_mHad'], help='Control regions to use', nargs='+')

parser.add_argument ('--correlate_tkPVfrac', default=False, action='store_true', help='Correlate tkPVfrac in all categories.')
parser.add_argument ('--freezeFF', default=False, action='store_true', help='Freeze form factors to central value.')
parser.add_argument ('--freeMuBr', default=True, help='Make muon branching fraction with a rate parameter (flat prior).')
parser.add_argument ('--asimov', default=False, action='store_true', help='Use Asimov dataset insted of real data.')
parser.add_argument ('--noMCstats', default=False, action='store_true', help='Do not include MC stat systematic.')

parser.add_argument ('--dumpWeightsTree', default=False, action='store_true', help='Dump tree with weights for skimmed events.')

availableSteps = ['clean', 'histos', 'preFitPlots', 'shapeVarPlots',
                  'card', 'workspace',
                  'bias', 'scan', 'catComp',
                  'fitDiag', 'postFitPlots',
                  'uncBreakdownScan', 'uncBreakdownTable',
                  'externalize',
                  'impacts', 'GoF']
defaultPipelineSingle = ['histos', 'preFitPlots', 'card', 'workspace', 'scan', 'GoF', 'fitDiag', 'postFitPlots', 'uncBreakdownScan']
defaultPipelineComb = ['preFitPlots', 'card', 'workspace', 'scan', 'catComp', 'uncBreakdownTable', 'GoF', 'fitDiag', 'postFitPlots', 'uncBreakdownScan']
# histos preFitPlots shapeVarPlots card workspace scan fitDiag postFitPlots uncBreakdownScan GoF
parser.add_argument ('--step', '-s', type=str, default=[], choices=availableSteps, help='Analysis steps to run.', nargs='+')
parser.add_argument ('--submit', default=False, action='store_true', help='Submit a job instead of running the call interactively.')
parser.add_argument ('--runInJob', default=False, action='store_true', help='Not for user, reserved to jobs calls.')

parser.add_argument ('--validateCard', default=False, action='store_true', help='Run combine card validation.')
parser.add_argument ('--decorrelateFFpars', default=False, action='store_true', help='Decorrelte form factors parameters')

parser.add_argument ('--forceRDst', default=False, action='store_true', help='Perform fit fixing R(D*) to 0.295')
parser.add_argument ('--seed', default=6741, type=int, help='Seed used by Combine')
parser.add_argument ('--RDstLims', default=[], type=float, help='Initial boundaries for R(D*).', nargs='+')

# Bias options
parser.add_argument ('--runBiasToys', default=False, action='store_true', help='Only generate toys and run scans for bias, do not collect results.')
parser.add_argument ('--nToys', default=10, type=int, help='Number of toys to run')
parser.add_argument ('--toysRDst', default=0.295, type=float, help='R(D*) value used to generate the toys.')
parser.add_argument ('--runBiasAnalysis', default=False, action='store_true', help='Only analyze bias scans which have been previously produced.')

# Scan options
parser.add_argument ('--scanStrategy', default=1, type=int, help='Minimizer strategy for the scan.')
parser.add_argument ('--maskScan', type=str, default=[], nargs='+', help='Channels to mask during likelyhood scan. If this list is non empty, the full card is used (default is fitregionsOnly).')
parser.add_argument ('--scanTag', type=str, default='')
parser.add_argument ('--freezeParsScan', type=str, default=[], nargs='+')

# Externalization options
parser.add_argument ('--externPars', default=['B2DstCLNeig1', 'B2DstCLNeig2', 'B2DstCLNeig3'], type=str, help='Parameters to externalize.', nargs='+')
parser.add_argument ('--externSigma', default=1., type=float, help='Externalization sigmas.')
parser.add_argument ('--externTag', default='FF', type=str, help='Externalization tag.')
parser.add_argument ('--externCenter', default='postFit', type=str, choices=['preFit', 'postFit'], help='Externalization tag.')

# Impacts options
parser.add_argument ('--collectImpacts', default=False, action='store_true', help='Only collect impact fits which have been previously run')
parser.add_argument ('--subOnlyImpacts', default=False, action='store_true', help='Only submit impact fits, do not collect results')

# Goodness Of Fit options
parser.add_argument ('--algoGoF', type=str, default=['Sat', 'KS'], choices=['Sat', 'AD', 'KS'], help='Algorithm to be used for the goodness of fit test', nargs='+')
parser.add_argument ('--maskEvalGoF', type=str, default=[], nargs='+', help='Additional channels to mask during GoF evaluation')
parser.add_argument ('--tagGoF', type=str, default='all')

parser.add_argument ('--showPlots', default=False, action='store_true', help='Show plots by setting ROOT batch mode OFF (default ON)')
parser.add_argument ('--showCard', default=False, action='store_true', help='Dump card on std outoput')
parser.add_argument ('--verbose', default=0, type=int, help='Run verbosity.')

args = parser.parse_args()
if args.HELP:
    parser.print_help()
    exit()

if len(args.step) == 0:
    if args.category == 'comb':
        args.step = defaultPipelineComb
    else: args.step = defaultPipelineSingle

    if args.cardTag == 'test_' and (not args.submit) and (not args.runInJob):
        args.step = ['clean'] + args.step

    if not args.unblinded:
        for s in args.step:
            if 'uncBreakdown' in s:
                args.step.remove(s)

    # if args.runInJob and not args.category == 'comb':
    #     args.step.append('shapeVarPlots')
    print 'Running default steps: ' + ', '.join(args.step)

schemeFF = args.schemeFF
if not args.showPlots:
    rt.gROOT.SetBatch(True)
    plt.ioff()
    plt.switch_backend('Agg')

if len(args.RDstLims) > 2:
    raise
elif len(args.RDstLims) == 2:
    if args.RDstLims[1] <= args.RDstLims[0]:
        raise

########################### ---- Define standards ----------- ########################
categoriesToCombine = ['low', 'mid', 'high']

binning = {'q2': array('d', [0, 3.5, 6, 9.4, 12])}

SM_RDst = 0.295
expectedLumi = {'Low':6.4, 'Mid':20.7, 'High':26.4, 'Single':20.7} #fb^-1
if args.lumiMult != 1.:
    print 'Multipling the expected luminosity by {:.1f}'.format(args.lumiMult)
    for n in expectedLumi.keys():
        expectedLumi[n] *= args.lumiMult
    print expectedLumi

FreeParFF = {
   'CLN': ['R0', 'eig1', 'eig2', 'eig3'],
   'BLPR': ['eig1', 'eig2', 'eig3', 'eig4', 'eig5', 'eig6'],
   'NoFF': []
}[schemeFF]

processOrder = [
    'tau', 'mu',
    'Bu_MuDstPi', 'Bd_MuDstPi',
    'Bd_MuDstPiPi', 'Bu_MuDstPiPi',
    'Bu_TauDstPi', 'Bd_TauDstPi',
    'Bd_TauDstPiPi', 'Bu_TauDstPiPi',
    'Bs_MuDstK', 'Bs_TauDstK',
    'Bd_DstDu', 'Bu_DstDu',
    'Bd_DstDd', 'Bu_DstDd',
    'Bd_DstDs', 'Bs_DstDs',
    # 'dataSS_DstMu'
]


samples_Bd = [p  for p in processOrder if (p[:2] == 'Bd' or p in ['tau', 'mu'])]
samples_Bu = [p  for p in processOrder if p[:2] == 'Bu']
samples_Bs = [p  for p in processOrder if p[:2] == 'Bs']


########################### --------------------------------- #########################
if args.asimov:
    CMS_lumi.extraText = "     Simulation Preliminary"

def createCardName(a):
    c = a.cardTag + a.category + '_' + a.schemeFF
    if a.decorrelateFFpars:
        c += 'decorr'
    if a.freezeFF:
        c += 'frozen'
    if a.useMVA:
        c += '_MVA'+useMVA
    if a.asimov:
        c += '_Asimov'
        if args.lumiMult != 1.:
            c += '{:.0f}'.format(args.lumiMult)
    if not a.unblinded:
        c += '_blinded'
    if a.noMCstats:
        c += '_NoMCstats'
    if not a.freeMuBr:
        c += '_muBrPDG'
    # else:
    #     c += '_freeMuBr'
    return c

card_name = createCardName(args)
print 'Card name:', card_name
print 'Control regions:', args.controlRegions

basedir = os.path.dirname(os.path.abspath(__file__)).replace('Combine', '')

outdir = basedir + 'Combine/results/' + card_name
if not os.path.isdir(outdir):
    os.system('mkdir -p ' + outdir + '/fig')
card_location = basedir + 'Combine/cards/{}.txt'.format(card_name)
histo_file_dir = basedir + 'data/_root/histos4combine/'

userName = basedir[basedir.find('/user/')+6:].split('/')[0]
webFolder = '/storage/af/user/'+userName+'/public_html/BPH_RDst/Combine/' + card_name
if not os.path.isdir(webFolder):
    os.makedirs(webFolder)
    os.system('cp '+webFolder+'/../index.php '+webFolder)
# Log command line call and arguments
if not args.submit:
    with open(webFolder + '/callsLog.txt', 'a') as f:
        f.write(50*'#'+'\n')
        f.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write( ' '.join(sys.argv) + '\n')
        f.write(5*'>'+' Arguments\n')
        f.write( yaml.dump(args.__dict__) )
        f.write(5*'<'+'\n')
        f.write(50*'-'+ '\n')

def runCommandSafe(command, printCommand=True):
    if printCommand:
        print command
    status, output = commands.getstatusoutput(command)

    def raiseFlag(inputText):
        flag = 'Warning: Did not find a parameter' in inputText
        flag = flag or ('WARNING: cannot freeze nuisance' in inputText)
        flag = flag or ('WARNING: MultiDimFit failed' in inputText)
        flag = flag or ('ERROR' in inputText and not 'Messages of type ERROR : 0' in inputText)
        flag = flag or ('There was a crash.' in inputText)
        # flag = flag or ('Error' in inputText)
        return flag
    flagged = raiseFlag(output)
    if status or flagged:
        print output, '\n\n'
        print '==================================================================='
        print '====================== Breaking the execution ====================='
        print '==================================================================='
        if status:
            print '\033[1m\x1b[31mStatus:\x1b[0m', status

        if flagged:
            for line in output.split('\n'):
                if raiseFlag(line):
                    print '\033[1m\x1b[31mFlagged line:\x1b[0m', line
        raise
    return output
########################### -------- Clean previous results ------------------ #########################

def cleanPreviousResults():
    os.system('rm -v '+card_location.replace('.txt', '*'))

    os.system('rm '+histo_file_dir+os.path.basename(card_location).replace('.txt', '_*'))

    os.system('rm -rf '+outdir)
    os.system('mkdir -p ' + outdir + '/fig')

    os.system('rm -rf '+webFolder)
    os.makedirs(webFolder)
    os.system('cp '+webFolder+'/../index.php '+webFolder)


########################### -------- Create histrograms ------------------ #########################
controlRegSel = {}
def selfun__TkPlus(ds):
    sel = np.logical_and(ds['N_goodAddTks'] == 1, ds['tkCharge_0'] > 0)
    return sel
controlRegSel['p_'] = selfun__TkPlus

def selfun__TkMinus(ds):
    sel = np.logical_and(ds['N_goodAddTks'] == 1, ds['tkCharge_0'] < 0)
    return sel
controlRegSel['m_'] = selfun__TkMinus

def selfun__TkPlusMinus(ds):
    sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == 0, ds['N_goodAddTks'] == 2)
    sel = np.logical_and(ds['massVisTks'] < 5.55, sel)
    return sel
controlRegSel['pm'] = selfun__TkPlusMinus

def selfun__TkMinusMinus(ds):
    sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == -2, ds['N_goodAddTks'] == 2)
    sel = np.logical_and(ds['massVisTks'] < 5.3, sel)
    return sel
controlRegSel['mm'] = selfun__TkMinusMinus

def selfun__TkPlusPlus(ds):
    sel = np.logical_and(ds['tkCharge_0']+ds['tkCharge_1'] == +2, ds['N_goodAddTks'] == 2)
    sel = np.logical_and(ds['massVisTks'] < 5.3, sel)
    return sel
controlRegSel['pp'] = selfun__TkPlusPlus

corrScaleFactors = {}
def loadDatasets(category, loadRD):
    print 'Loading MC datasets'
    #They all have to be produced with the same pileup
    # candDir='ntuples_B2DstMu_mediumId_lostInnerHits'
    candDir='ntuples_B2DstMu_220225'
    print 'Using candDir =', candDir
    print 'Using skim = skimmed'+args.skimmedTag
    MCsample = {
    ######## Signals
    'tau': DSetLoader('Bd_TauNuDst', candDir=candDir, skimmedTag=args.skimmedTag),
    'mu': DSetLoader('Bd_MuNuDst', candDir=candDir, skimmedTag=args.skimmedTag),
    ######## D** background
    'Bu_MuDstPi': DSetLoader('Bu_MuNuDstPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_MuDstPi': DSetLoader('Bd_MuNuDstPi', candDir=candDir, skimmedTag=args.skimmedTag),
    # 'Bd_MuDstPiPi': DSetLoader('Bd_MuNuDstPiPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_MuDstPiPi': DSetLoader('Bd_MuNuDstPiPi_v3', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bu_MuDstPiPi': DSetLoader('Bu_MuNuDstPiPi_v3', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bu_TauDstPi': DSetLoader('Bu_TauNuDstPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_TauDstPi': DSetLoader('Bd_TauNuDstPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_TauDstPiPi': DSetLoader('Bd_TauNuDstPiPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bu_TauDstPiPi': DSetLoader('Bu_TauNuDstPiPi', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bs_MuDstK': DSetLoader('Bs_MuNuDstK', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bs_TauDstK': DSetLoader('Bs_TauNuDstK', candDir=candDir, skimmedTag=args.skimmedTag),
    ######## D*Hc background
    'Bd_DstDu': DSetLoader('Bd_DstDu', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bu_DstDu': DSetLoader('Bu_DstDu', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_DstDd': DSetLoader('Bd_DstDd', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bu_DstDd': DSetLoader('Bu_DstDd', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bd_DstDs': DSetLoader('Bd_DstDs', candDir=candDir, skimmedTag=args.skimmedTag),
    'Bs_DstDs': DSetLoader('Bs_DstDs', candDir=candDir, skimmedTag=args.skimmedTag),
    }

    dSet = {}
    dSetTkSide = {}
    mcType = 'bare' if args.bareMC else 'corr'
    print 'mcType:', mcType
    if not args.maxEventsToLoad is None:
        print 'Limiting events per MC sample to', args.maxEventsToLoad

    for n, s in MCsample.iteritems():
        if not n in processOrder:
            print n, 'not declarted in processOrder'
            raise
        dSet[n] = pd.DataFrame(rtnp.root2array(s.skimmed_dir + '/{}_{}.root'.format(category.name, mcType),
                                               stop=args.maxEventsToLoad
                                              )
                              )
        dSetTkSide[n] = pd.DataFrame(rtnp.root2array(s.skimmed_dir + '/{}_trkCtrl_{}.root'.format(category.name, mcType),
                                                     stop=args.maxEventsToLoad
                                                    )
                                    )

    dataDir = '/storage/af/group/rdst_analysis/BPhysics/data/cmsRD'
    # locRD = dataDir+'/skimmed'+args.skimmedTag+'/B2DstMu_SS_211014_{}'.format(category.name)
    # locRD = dataDir+'/skimmed/B2DstMu_SS_211014_{}'.format(category.name)
    # dSet['dataSS_DstMu'] = pd.DataFrame(rtnp.root2array(locRD + '_corr.root'))
    # dSetTkSide['dataSS_DstMu'] = pd.DataFrame(rtnp.root2array(locRD + '_trkCtrl_corr.root'))


    if loadRD:
        print 'Loading real data datasets'

        creation_date = '220220'
        locRD = dataDir+'/skimmed'+args.skimmedTag+'/B2DstMu_{}_{}'.format(creation_date, category.name)
        dSet['data'] = pd.DataFrame(rtnp.root2array(locRD + '_corr.root'))
        dSetTkSide['data'] = pd.DataFrame(rtnp.root2array(locRD + '_trkCtrl_corr.root'))
        # Not used
        # datasets_loc = glob(dataDir + '/ParkingBPH*/*RDntuplizer_B2DstMu_{}_CAND.root'.format(creation_date))
        # lumi_tot = getLumiByTrigger(datasets_loc, category.trg, verbose=True)

    if args.dumpWeightsTree:
        print 'Skipping on the flight cuts (if any).'
    else:
        addCuts = [
        ['M2_miss', 0.4, 1e3],
        ['mu_eta', -0.8, 0.8],
        # ['mu_eta', -1.5, 1.5],
        # ['B_eta', -1., 1.],
        # ['K_pt', 1., 1e3],
        # ['pi_pt', 1., 1e3],
        # ['pis_pt', 1., 1e3],
        ['mu_lostInnerHits', -2, 1],
        ['K_lostInnerHits', -2, 1],
        ['pi_lostInnerHits', -2, 1],
        ['pis_lostInnerHits', -2, 1],
        # ['ctrl_tk_pval_0', 0.2, 1.0],
        # ['ctrl_tk_pval_1', 0.2, 1.0],
        # ['ctrl_pm_massVisTks', 0, 3.8],
        # ['ctrl_pm_massHadTks', 2.6, 10],
        # ['ctrl_pm_index', 3, 0],
        ]
        if len(addCuts) > 0:
            with open(webFolder + '/callsLog.txt', 'a') as f:
                f.write(5*'>'+' Cuts at loading time\n')
                for v, m, M in addCuts:
                    if 'index' in v:
                        f.write('{} % {} <= {}\n'.format(v, m, M))
                    else:
                        f.write('{} < {} < {}\n'.format(m, v, M))
                f.write(5*'<'+'\n')
                f.write(50*'-'+ '\n')
            for k in dSet.keys():
                sel = np.ones_like(dSet[k]['q2']).astype(np.bool)
                for var, low, high in addCuts:
                    if var.startswith('ctrl_'): continue
                    if var not in dSet[k].columns:
                        print var, 'not in', k, 'main dataset'
                        raise
                    sel = np.logical_and(sel, np.logical_and(dSet[k][var] > low, dSet[k][var] < high))

                dSet[k] = dSet[k][sel]
                corrScaleFactors[k] = np.sum(sel)/float(sel.shape[0])

                sel = np.ones_like(dSetTkSide[k]['q2']).astype(np.bool)
                for var, low, high in addCuts:
                    region = 'all'
                    if var.startswith('ctrl_'):
                        var = var.replace('ctrl_', '')
                        if re.match('[pm][pm_]_', var[:3]):
                            region = var[:2]
                            var = var[3:]
                    if var not in dSetTkSide[k].columns:
                        print var, 'not in', k, 'control regions dataset'
                        raise
                    if var == 'index':
                        thisSel = np.mod(dSetTkSide[k][var], low) <= high
                    else:
                        thisSel = np.logical_and(dSetTkSide[k][var] > low, dSetTkSide[k][var] < high)
                    if region == 'all':
                        sel = np.logical_and(sel, thisSel)
                    else:
                        auxSel = np.logical_or( np.logical_not(controlRegSel[region](dSetTkSide[k])), thisSel)
                        sel = np.logical_and(sel, auxSel)

                dSetTkSide[k] = dSetTkSide[k][sel]
                corrScaleFactors[k+'_tk'] = np.sum(sel)/float(sel.shape[0])


    return MCsample, dSet, dSetTkSide

def computeBrVarWeights(ds, selItems={}, relScale=0.2, centralVal=1., keepNorm=False, absVal=True):
    sel = np.ones_like(ds['mu_pt']).astype(np.bool)
    for var, val in selItems.iteritems():
        if absVal:
            sel = np.logical_and(np.abs(ds[var].astype(np.int)) == val, sel)
        else:
            sel = np.logical_and(ds[var].astype(np.int) == val, sel)
    w = np.where(sel, centralVal, 1.)
    up = np.where(sel, centralVal*(1.+relScale), 1.)
    down = np.where(sel, centralVal*max(0, 1.-relScale), 1.)
    if keepNorm:
        up = (float(up.shape[0])/np.sum(up)) * up
        down = (float(down.shape[0])/np.sum(down)) * down
    return w, up/w, down/w

def computeWidthVarWeights(ds, selItems=[], newGamma=None, relScale=0.1, keepNorm=False): #Gamma modification factor
    # See EvtGen BW implementation: https://github.com/alberto-sanchez/evtgen-cms/blob/master/src/EvtGenBase/EvtBreitWignerPdf.cpp#L41
    # selItems=[ [pdgId, mass, Gamma] ]
    w = np.ones_like(ds['mu_pt'])
    up = np.ones_like(ds['mu_pt'])
    down = np.ones_like(ds['mu_pt'])
    for i, (pdgId, mass, gamma) in enumerate(selItems):
        # print pdgId, mass, gamma
        dx2 = np.clip(np.square(ds['MC_MassCharmedBDaughter'] - mass), 0, 9*(gamma**2)) # Non relativistic
        # dx2 = np.clip(np.square(np.square(ds['MC_MassCharmedBDaughter']) - mass**2), 0, 9*((mass*gamma)**2))

        if not (newGamma is None) and not (newGamma[i] is None):
                gNew = newGamma[i]
                wCentral = ((dx2 + (gamma**2)/4.0)*gNew)/(gamma*(dx2 + (gNew**2)/4.0)) # Non relativistic
                # wCentral = (gNew/gamma) * ( (dx2 + (mass*gamma)**2) / (dx2 + (mass*gNew)**2) )
                gUp = gNew*(1+relScale)
                gDown = gNew*(1-relScale)
        else:
            wCentral = np.ones_like(dx2)
            gUp = gamma*(1+relScale)
            gDown = gamma*(1-relScale)

        # Non relativistic
        wUp = ((dx2 + (gamma**2)/4.0)*gUp)/(gamma*(dx2 + (gUp**2)/4.0))
        wDown = ((dx2 + (gamma**2)/4.0)*gDown)/(gamma*(dx2 + (gDown**2)/4.0))
        # wUp = (gUp/gamma) * ( (dx2 + (mass*gamma)**2) / (dx2 + (mass*gUp)**2) )
        # wDown = (gDown/gamma) * ( (dx2 + (mass*gamma)**2) / (dx2 + (mass*gDown)**2) )

        sel = np.abs(ds['MC_DstMotherPdgId'].astype(np.int)) == np.abs(pdgId)
        w = np.where(sel, wCentral, w)
        up = np.where(sel, wUp, up)
        down = np.where(sel, wDown, down)

    if keepNorm:
        w *= float(w.shape[0])/np.sum(w)
        up *= float(w.shape[0])/np.sum(up)
        down *=  float(w.shape[0])/np.sum(down)
    return w, up/w, down/w

def computeTksPVweights(ds, relScale=0.05, centralVal=0.39/0.10):
    # selPdgID0 = np.logical_and(np.abs(ds['MC_tkMotherPdgId_0']) < 6, ds['MC_tkMotherPdgId_0'] != 0)
    # selPdgID0 = np.logical_or(selPdgID0, ds['MC_tkMotherPdgId_0']==2212)
    # selPdgID0 = np.logical_and(selPdgID0, ds['MC_tkFlag_0'] == 1)
    # selPdgID1 = np.logical_and(np.abs(ds['MC_tkMotherPdgId_1']) < 6, ds['MC_tkMotherPdgId_1'] != 0)
    # selPdgID1 = np.logical_or(selPdgID1, ds['MC_tkMotherPdgId_1']==2212)
    # selPdgID1 = np.logical_and( selPdgID1, ds['MC_tkFlag_1'] == 1)

    selPdgID0 = np.logical_and( ds['MC_tkFlag_0'] == 1, ds['MC_tkFromMainB_0'] == 0)
    selPdgID1 = np.logical_and( ds['MC_tkFlag_1'] == 1, ds['MC_tkFromMainB_1'] == 0)
    exponent = selPdgID0.astype(np.int) + selPdgID1.astype(np.int)
    w = np.power(centralVal, exponent)
    up = np.power(centralVal*(1+relScale), exponent)/w
    down = np.power(centralVal*max(0, 1-relScale), exponent)/w
    return w, up, down

def createHistograms(category):
    MCsample, dSet, dSetTkSide = loadDatasets(category, not args.asimov)
    mcType = 'bare' if args.bareMC else 'corr'
    ######################################################
    ########## Load calibrations
    ######################################################
    from pileup_utilities import pileupReweighter
    skimmedFile_loc = MCsample['tau'].skimmed_dir + '/{}_{}.root'.format(category.name, mcType)
    puReweighter = pileupReweighter(skimmedFile_loc, 'hAllNTrueIntMC', trg=category.trg)

    fname = '/storage/af/group/rdst_analysis/BPhysics/data/calibration/beamSpot/crystalball_calibration_v1_'+category.name.capitalize()+'.yaml'
    beamSpotParam = yaml.load(open(fname, 'r'))

    dataDir = '/storage/af/group/rdst_analysis/BPhysics/data'
    decayBR = pickle.load(open(dataDir+'/forcedDecayChannelsFactors_v2.pickle', 'rb'))

    loc = dataDir+'/calibration/triggerScaleFactors/'
    fTriggerSF = rt.TFile.Open(loc + 'HLT_' + category.trg + '_SF_v22_count.root', 'READ')
    hTriggerSF = fTriggerSF.Get('hSF_HLT_' + category.trg)
    def computeTrgSF(ds, hSF, selection=None):
        trgSF = np.ones_like(ds['q2'])
        trgSFUnc = np.zeros_like(ds['q2'])
        ptmax = hSF.GetXaxis().GetXmax() - 0.01
        ipmax = hSF.GetYaxis().GetXmax() - 0.01
        etamax = hSF.GetZaxis().GetXmax() - 0.01
        x = np.column_stack((ds['mu_pt'], ds['mu_eta'], ds['mu_sigdxy']))
        if not selection is None:
            x = x[selection]
        for i, (pt, eta, ip) in enumerate(x):
            ix = hSF.GetXaxis().FindBin(min(ptmax, pt))
            iy = hSF.GetYaxis().FindBin(min(ipmax, ip))
            iz = hSF.GetZaxis().FindBin(min(etamax, np.abs(eta)))
            trgSF[i] = hSF.GetBinContent(ix, iy, iz)
            ib = hSF.GetBin(ix, iy, iz)
            trgSFUnc[i] = hSF.GetBinError(ib)
            if trgSF[i] == 0:
                print pt, ip, np.abs(eta)
                raise

        # Divide them for the weight so later you can simply multiply back to get the value
        up = 1 + trgSFUnc/trgSF
        down = 1 - trgSFUnc/trgSF
        return trgSF, up, down

    fMuonIDSF = rt.TFile.Open(dataDir+'/calibration/muonIDscaleFactors/Run2018ABCD_SF_MuonID_Jpsi.root', 'READ')
    hMuonIDSF = fMuonIDSF.Get('NUM_MediumID_DEN_genTracks_pt_abseta')
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

    if args.calBpT == 'none':
        print 'Not using any B pT calibration'
    elif args.calBpT == 'poly':
        cal_pT_Bd = kinCalReader(calibration_file=dataDir+'/calibration/kinematicCalibration_Bd/pt_polyCoeff_'+category.name+'_v2.pkl')

    cal_eta_B = kinCalReader(calibration_file=dataDir+'/calibration/kinematicCalibration_Bd/eta_polyCoeff_'+category.name+'_v2.pkl')


    def computeKinCalWeights(ds, var, tag, kinCal):
        if kinCal.kind == 'poly':
            # The denominator (sum of weights) for this weights is not known but it cancel out in the ratio
            w = kinCal.getWeights(ds[var], shape=0)
            if np.sum(w==0):
                print np.sum(w==0)
                raise

            varDic = {}
            for iShape in range(1, kinCal.nVar+1):
                varDic[tag+'_lam{}Down'.format(iShape)] = kinCal.getWeights(ds[var], shape= -iShape, scale=1.)/w
                varDic[tag+'_lam{}Up'.format(iShape)] = kinCal.getWeights(ds[var], shape= iShape, scale=1.)/w

            return w, varDic
        elif kinCal.kind == 'ratio':
            w = kinCal.f['C'](ds[var])
            if np.sum(w==0):
                print np.sum(w==0)
                raise
            up = kinCal.f['Up'](ds[var])/w
            down = kinCal.f['Down'](ds[var])/w
            return w, up, down
        else:
            print 'Unknown calibration'
            raise


    parsSoftTracks = {'s':[0.2, 0.15], 'w':[0.9, 0.05]}
    def fSoftTrackEff(x, w, s):
        # Model beta
        x = [0, 3.5, 1.5]
        yLin = x[2]*(1 - w)/x[1] + w
        y = [w, 1., yLin + s*(1-yLin) ]
        beta = np.polyfit(x, y, 2)

        sel = np.logical_or(x < 0.2, x > 3.5)
        return np.where(sel, np.ones_like(x), np.polyval(beta, x))

    softPtUnc= [[0.5, 0.6, 0.10],
                [0.6, 0.7, 0.07],
                [0.7, 0.8, 0.05],
                [0.8, 0.9, 0.04],
                [0.9, 1.0, 0.03],
                [1.0, 1.2, 0.02],
                [1.2, 2.0, 0.01],
                ]
    def binnnedSoftTrackEff(x, bin, size):
        sel = np.logical_or(x < softPtUnc[bin][0], x > softPtUnc[bin][1])
        return np.where(sel, np.ones_like(x), 1+size*softPtUnc[bin][2])

    def weightsSoftTrackEff(ds, ptList, w=None, s=None, bin=None, size=None):
        if w is None or s is None:
            weight = binnnedSoftTrackEff(ds[ptList[0]], bin, size)
            for v in ptList[1:]:
                weight *= binnnedSoftTrackEff(ds[v], bin, size)
        else:
            weight = fSoftTrackEff(ds[ptList[0]], w, s)
            for v in ptList[1:]:
                weight *= fSoftTrackEff(ds[v], w, s)
        return weight


    if args.useMVA:
        fname = dataDir+'../plot_scripts/kinObsMVA/clfGBC_tauVall_{}{}.p'.format(args.useMVA, category.name)
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
    data_to_MC_prescale = 0.8

    ######################################################
    ########## Signal region
    ######################################################
    n_q2bins = len(binning['q2'])-1
    negSide = [-2.5, -1.5, -1.0, -0.6, -0.4, -0.2]
    binning['M2_miss'] = [
            array('d', negSide + [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 4] ),
            array('d', negSide + [0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.5, 0.2)) + [8] ),
            array('d', negSide + list(np.arange(0, 6, 0.2)) + [8] ),
            array('d', negSide + list(np.arange(0, 7.8, 0.2)) + [8] ),
        ]
    # binning['M2_miss'] = [
    #         array('d', [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 4] ),
    #         array('d', [0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.5, 0.2)) + [8] ),
    #         array('d', list(np.arange(0, 6, 0.2)) + [8] ),
    #         array('d', list(np.arange(0, 7.8, 0.2)) + [8] ),
    #     ]
    binning['Est_mu'] = [
            array('d', [0.3] + list(np.arange(0.5, 2.2, 0.05)) + [2.3] ),
            array('d', [0.3] + list(np.arange(0.5, 2.2, 0.05)) + [2.3] ),
            array('d', [0.3] + list(np.arange(0.5, 2.2, 0.05)) + [2.2] ),
            [24, 0.3, 2.0],
        ]

    binning['U_miss'] = 4*[[30, -0.1, 0.18]]

    binning['mu_pt'] = n_q2bins*[{'Low': array('d', list(np.arange(7.2, 9.2, 0.05))+[9.2] ),
                                  'Mid': array('d', list(np.arange(9.2, 12.2, 0.05)) +[12.2] ),
                                  'High': array('d', list(8+np.logspace(np.log10(12.2-8), np.log10(60), 30)) )
                                 }[category.name]]

    binning['Dst_pt'] = n_q2bins*[{'Low':  array('d', list(np.arange(3, 35, 1)) ),
                                   'Mid':  array('d', list(np.arange(4, 40, 1)) ),
                                   'High': array('d', list(np.arange(5, 50, 1)) )
                                  }[category.name]]

    binning['K_pt'] = n_q2bins*[{'Low':  array('d', list(np.arange(0.8, 15, 0.4)) ),
                                   'Mid':  array('d', list(np.arange(0.8, 20, 0.4)) ),
                                   'High': array('d', list(np.arange(0.8, 30, 0.4)) )
                                  }[category.name]]

    binning['pi_pt'] = n_q2bins*[{'Low':  array('d', list(np.arange(0.8, 15, 0.4)) ),
                                   'Mid':  array('d', list(np.arange(0.8, 20, 0.4)) ),
                                   'High': array('d', list(np.arange(0.8, 30, 0.4)) )
                                  }[category.name]]

    binning['pis_pt'] = n_q2bins*[{'Low':  array('d', list(np.arange(0.5, 3.5, 0.1)) ),
                                   'Mid':  array('d', list(np.arange(0.5, 4, 0.1)) ),
                                   'High': array('d', list(np.arange(0.5, 5, 0.1)) )
                                  }[category.name]]

    binning['mass_D0pismu'] = n_q2bins*[[50, 2.1, 5.28]]

    negSide = [-2.5, -1.0, -0.4, -0.2]
    binning_2D = [
        [
            array('d', negSide + [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5, 4] ),
            array('d', [0.3] + list(np.arange(0.7, 2.3, 0.3)) + [2.3] )
        ],
        [
            array('d', negSide + [0.0, 0.1, 0.2, 0.3] + list(np.arange(0.4, 3.0, 0.2)) + [8] ),
            array('d', [0.3] + list(np.arange(0.7, 2.2, 0.3)) )
        ],
        [
            array('d', negSide + list(np.arange(0, 5.6, 0.4)) + [8] ),
            array('d', [0.3] + list(np.arange(0.5, 2.1, 0.3)) + [2.1] )
        ],
        [
            array('d', negSide + list(np.arange(0, 7.6, 0.4)) + [8] ),
            array('d', list(np.linspace(0.3, 2.0, 10)) )
        ]

    ]
    binning['B_pt'] = {'Low': array('d', list(np.arange(10, 75, 2)) ), 'Mid': array('d', list(np.arange(14, 90, 2)) ), 'High': array('d', list(np.arange(18, 110, 2)))}[category.name]

    binning['mu_eta'] = array('d', list(np.arange(-1.6, 1.61, 0.05)) )
    binning['B_eta'] = array('d', list(np.arange(-1.9, 1.91, 0.05)) )
    binning['K_eta'] = array('d', list(np.arange(-2.1, 2.11, 0.05)) )
    binning['pi_eta'] = array('d', list(np.arange(-2.1, 2.11, 0.05)) )
    binning['pis_eta'] = array('d', list(np.arange(-2.1, 2.11, 0.05)) )

    binning['specQ2'] = array('d', list(np.arange(-2, 11.4, 0.2)))

    binning['mass_piK'] = [75, 1.81, 1.92]
    binning['deltaM_DstD'] = [75, 0.1434, 0.1475]

    if args.unblinded:
        binning['MVA'] = array('d', list(np.arange(0., 0.83, 0.02)))
    else:
        binning['MVA'] = array('d', list(np.arange(0., 0.51, 0.03)))

    additionalVariablesToPlot = ['B_pt', 'B_eta', 'specQ2']
    additionalVariablesToPlot += ['mu_eta', 'K_eta', 'pi_eta', 'pis_eta']
    additionalVariablesToPlot += ['mass_piK', 'deltaM_DstD']


    totalCounting = [0,0]
    print '---------> Fill signal region histograms'
    for n in processOrder:
        ds = dSet[n]
        if n == 'data': continue
        print '\n----------->', n, '<-------------'
        wVar = {}
        weights = {}
        if n == 'dataSS_DstMu':
            nTotSelected = ds['q2'].shape[0]
            nTotExp = ds['q2'].shape[0]
        else:
            sMC = MCsample[n]

            nTotSelected = ds['q2'].shape[0]
            print 'N tot selected: {:.1f}k'.format(1e-3*nTotSelected)
            totalCounting[1] += 1e-3*nTotSelected
            nGenExp = sMC.effMCgen['xsec'][0] * expectedLumi[category.name] * data_to_MC_prescale
            eff = [1, 0]
            for f, df in [sMC.effMCgen['effGEN'], decayBR[n], sMC.effCand['effCAND'], sMC.getSkimEff(category.name+'_'+mcType)]:
                eff[0] *= f
                eff[1] += np.square(df/f)
            eff[1] = eff[0] * np.sqrt(eff[1])
            if n in corrScaleFactors.keys():
                eff[0] *= corrScaleFactors[n]
                print 'Using scale factor from a posteriori selection: {:.3f}'.format(corrScaleFactors[n])
            nTotExp = nGenExp*eff[0]
        print 'N tot expected (before weights): {:.2f}k'.format(1e-3*nTotExp)

        if not 'data' in n:
            print 'Including pileup reweighting'
            weights['pileup'] = puReweighter.getPileupWeights(ds['MC_nInteractions'])

            print 'Including beam spot correction'
            weights['beamSpot'] = getBeamSpotCorrectionWeights(ds, beamSpotParam)

            print 'Including trigger corrections'
            nameSF = 'trg{}SF'.format(category.trg)
            weights[nameSF], wSfUp, wSfDw = computeTrgSF(ds, hTriggerSF)
            auxOnes = np.ones_like(wSfUp)
            # for i_eta in range(1, hTriggerSF.GetNbinsZ()+1):
            #     c_eta = hTriggerSF.GetZaxis().GetBinCenter(i_eta)
            #     w_eta = hTriggerSF.GetZaxis().GetBinWidth(i_eta)
            # for i_ip in range(1, hTriggerSF.GetNbinsY()+1):
            #     c_ip = hTriggerSF.GetYaxis().GetBinCenter(i_ip)
            #     w_ip = hTriggerSF.GetYaxis().GetBinWidth(i_ip)
            #     if c_ip + 0.5*w_ip <= category.minIP:
            #         continue
            for i_pt in range(1, hTriggerSF.GetNbinsX()+2):
                if i_pt > hTriggerSF.GetNbinsX() and category.name == 'High':
                    sel = ds['mu_pt'] > hTriggerSF.GetXaxis().GetXmax()
                else:
                    c_pt = hTriggerSF.GetXaxis().GetBinCenter(i_pt)
                    w_pt = hTriggerSF.GetXaxis().GetBinWidth(i_pt)
                    if (c_pt + 0.5*w_pt <= category.min_pt) or (c_pt - 0.5*w_pt >= category.max_pt):
                        continue

                    sel = np.abs(ds['mu_pt'] - c_pt) < w_pt
                # sel = np.logical_and(sel, np.abs(ds['mu_sigdxy'] - c_ip) < w_ip)
                # sel = np.logical_and(sel, np.abs(ds['mu_eta'] - c_eta) < w_eta)
                # binName = '_pt{}ip{}eta{}'.format(i_pt, i_ip, i_eta)
                # print 'Trg SF', i_pt, i_ip, i_eta, '-> selected {}'.format(np.sum(sel))
                # binName = '_pt{}ip{}'.format(i_pt, i_ip)
                binName = '_pt{}'.format(i_pt)
                wVar[nameSF+binName+'Up'] = np.where(sel, wSfUp, auxOnes)
                wVar[nameSF+binName+'Down'] = np.where(sel, wSfDw, auxOnes)

            print 'Including muon ID corrections'
            weights['muonIdSF'], _, _ = computeMuonIDSF(ds)

            print 'Including soft track pT corrections'
            # l = ['K_pt', 'pi_pt', 'pis_pt']
            # weights['softTrkEff'] = weightsSoftTrackEff(ds, l, parsSoftTracks['w'][0], parsSoftTracks['s'][0])
            # for mod in [+1, -1]:
            #     varName = 'Up' if mod > 0 else 'Down'
            #     w = parsSoftTracks['w'][0] + mod*parsSoftTracks['w'][1]
            #     wVar['softTrkEff_w'+varName] = weightsSoftTrackEff(ds, l, w, parsSoftTracks['s'][0])/weights['softTrkEff']
            #     s = parsSoftTracks['s'][0] + mod*parsSoftTracks['s'][1]
            #     wVar['softTrkEff_s'+varName] = weightsSoftTrackEff(ds, l, parsSoftTracks['w'][0], s)/weights['softTrkEff']
            partList = ['K_pt', 'pi_pt', 'pis_pt']
            for nBin in range(len(softPtUnc)):
                refPt = '{:.0f}'.format(np.round(np.mean(softPtUnc[nBin][:-1])*1e3))
                wVar['softTrkEff_'+refPt+'Up'] = weightsSoftTrackEff(ds, partList, bin=nBin, size=+1)
                wVar['softTrkEff_'+refPt+'Down'] = weightsSoftTrackEff(ds, partList, bin=nBin, size=-1)


            ############################
            # B phase space corrections
            ############################
            print 'Including B eta corrections'
            cname = 'B_eta'+category.name
            weights[cname], auxVarDic = computeKinCalWeights(ds, 'MC_B_eta', cname, cal_eta_B)
            wVar.update(auxVarDic)

            if (not args.calBpT == 'none') and (n in samples_Bd):
                print 'Including Bd pT corrections'
                cname = 'BdpT'+category.name
                if cal_pT_Bd.kind == 'ratio':
                    weights[cname], wVar[cname+'Up'], wVar[cname+'Down'] = computeKinCalWeights(ds, 'MC_B_pt', None, cal_pT_Bd)
                else:
                    weights[cname], auxVarDic = computeKinCalWeights(ds, 'MC_B_pt', cname, cal_pT_Bd)
                    wVar.update(auxVarDic)

        ############################
        # Form factor correction
        ############################
        if n in ['mu', 'tau'] and schemeFF != 'NoFF':
            print 'Including B-> D*EllNu FF corrections (Hammer)'
            weights['B2DstFF'] = ds['wh_'+schemeFF+'Central']*sMC.effCand['rate_den']/sMC.effCand['rate_'+schemeFF+'Central']
            for nPar in FreeParFF:
                for var in ['Up', 'Down']:
                    tag = schemeFF + nPar + var
                    wVar['B2Dst'+tag] = ds['wh_'+tag]/ds['wh_'+schemeFF+'Central']
                    wVar['B2Dst'+tag] *= sMC.effCand['rate_'+schemeFF+'Central']/sMC.effCand['rate_' + tag]
                    wVar['B2Dst'+tag][np.isnan(wVar['B2Dst'+tag])] = 0.

        if '_MuDstPi' in n:
            print 'Including B -> D**MuNu FF corrections (Hammer)'
            weights['purgeUwantedEvts'] = np.array(ds['DststProc_id'] >= 0).astype(np.int)
            DststProc_id = np.array(ds['DststProc_id']).astype(np.int)

            ratesRatio = np.ones(50)
            ratesRatio[1] = sMC.effCand['rate_den_D1']/sMC.effCand['rate_D1BLR_central']
            ratesRatio[2] = sMC.effCand['rate_den_D1st']/sMC.effCand['rate_D1stBLR_central']
            ratesRatio[3] = sMC.effCand['rate_den_D2st']/sMC.effCand['rate_D2stBLR_central']
            selDstst = np.logical_and(DststProc_id >= 1, DststProc_id <= 3)
            weights['B2DststFF'] = np.where(selDstst, ds['wh_Dstst_BLRCentral']*ratesRatio[DststProc_id], 1.)

            selDstst = np.logical_or(DststProc_id == 1, DststProc_id == 3)
            for nEig in [1,2,3,4]:
                for var in ['Up', 'Down']:
                    tag = 'DststN_BLReig'+str(nEig)+var
                    ratesRatio = np.ones(50)
                    ratesRatio[1] = sMC.effCand['rate_D1BLR_central']/sMC.effCand['rate_D1BLR_eig'+str(nEig)+var]
                    ratesRatio[3] = sMC.effCand['rate_D2stBLR_central']/sMC.effCand['rate_D2stBLR_eig'+str(nEig)+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_Dstst_BLRCentral'])*ratesRatio[DststProc_id], 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.

            selDstst = DststProc_id == 2
            for nEig in [1,2,3]:
                for var in ['Up', 'Down']:
                    tag = 'DststW_BLReig'+str(nEig)+var
                    ratesRatio = sMC.effCand['rate_D1stBLR_central']/sMC.effCand['rate_D1stBLR_eig'+str(nEig)+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_Dstst_BLRCentral'])*ratesRatio, 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.

        if n.endswith('_MuDstPiPi'):
            print 'Including B -> D(2S)MuNu FF corrections (Hammer)'
            DststProc_id = np.array(ds['DststProc_id']).astype(np.int)

            ratesRatio = np.ones(50)
            ratesRatio[4] = sMC.effCand['rate_den_D2sst']/sMC.effCand['rate_D2sstBLOP_central']
            ratesRatio[11] = ratesRatio[4]
            ratesRatio[12] = ratesRatio[4]
            ratesRatio[13] = ratesRatio[4]
            ratesRatio[23] = sMC.effCand['rate_den_D2s']/sMC.effCand['rate_D2sBLOP_central']

            selDstst = np.logical_and(DststProc_id >= 4, DststProc_id <= 23)
            weights['B2D2S_FF'] = np.where(selDstst, ds['wh_D2S_BLOPCentral']*ratesRatio[DststProc_id], 1.)

            for parName in ['RhoSq', 'chi11', 'chi21', 'chi31', 'eta1']:
                for var in ['Up', 'Down']:
                    tag = 'D2S_BLOP'+parName+var
                    ratesRatio = np.ones(50)
                    ratesRatio[4] = sMC.effCand['rate_D2sstBLOP_central']/sMC.effCand['rate_D2sstBLOP_'+parName+var]
                    ratesRatio[11] = ratesRatio[4]
                    ratesRatio[12] = ratesRatio[4]
                    ratesRatio[13] = ratesRatio[4]
                    ratesRatio[23] = sMC.effCand['rate_D2sBLOP_central']/sMC.effCand['rate_D2sBLOP_'+parName+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_D2S_BLOPCentral'])*ratesRatio[DststProc_id], 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.


        ############################
        # Dstst resonance mix
        ############################
        if not re.search('DstPi\Z', n) is None:
            print 'Including D**->D*Pi branching ratio and width variations'

            infate = 2
            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10423}, infate*0.2/3.0)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10413}, infate*0.14/1.40)
            wVar['brB_D2420MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2420MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20423}, infate*0.6/2.7)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20413}, infate*0.45/1.55)
            wVar['brB_D2430MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2430MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 425}, infate*0.24/1.01)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 415}, infate*0.06/0.34)
            wVar['brB_D2460MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2460MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_DstMotherPdgId': 521}, infate*0.8)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_DstMotherPdgId': 511}, infate*0.8)
            wVar['brB_DstPiMuNuUp'] = wNeuUp * wChUp
            wVar['brB_DstPiMuNuDown'] = wNeuDw * wChDw


            w, wVar['D2420_widthUp'], wVar['D2420_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[10423, 2.421, 0.0274],
                                                                                                 [10413, 2.423, 0.020]],
                                                                                       newGamma=[0.030, 0.025],
                                                                                       relScale=0.15)
            weights['D2420_width'] = w

            w, wVar['D2430_widthUp'], wVar['D2430_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[20423, 2.445, 0.250],
                                                                                                 [20413, 2.445, 0.250]],
                                                                                       newGamma=[0.3, 0.3],
                                                                                       relScale=0.1)
            weights['D2430_width'] = w

            w, wVar['D2460_widthUp'], wVar['D2460_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[425, 2.461, 0.049],
                                                                                                 [415, 2.460, 0.037]],
                                                                                       newGamma=[0.045, 0.045],
                                                                                       relScale=0.15)
            weights['D2460_width'] = w

        if not re.search('MuDstPiPi\Z', n) is None:
            print 'Including B -> MuNuD*PiPi width and Br variations'
            widthMods = [[x, 2.640, 0.300] for x in [100413, 100423, 100411, 100421]]
            uncN = 'Dst2S_width'
            weights[uncN], wVar[uncN+'Up'], wVar[uncN+'Down'] = computeWidthVarWeights(ds, selItems=widthMods, relScale=0.3)

            keepNorm=False
            nnn = 'brDstPiPi_NR'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 0}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 1}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1st'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 2}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2st'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 3}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 4}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D1Pi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 11}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D1stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 12}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 13}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2S_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 23}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1Pi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 31}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 32}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 33}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            # wNeu, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_1': 111}, relScale=1.0, centralVal=2., keepNorm=keepNorm)
            # wCh, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_1': 211}, relScale=1.0, centralVal=2., keepNorm=keepNorm)
            # wVar['brDstst_DststPiUp'] = wNeuUp * wChUp
            # wVar['brDstst_DststPiDown'] = wNeuDw * wChDw
            # weights['DstPiPi_test'] = wNeu*wCh
            #
            # w, wD1Up, wD1Dw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.1, keepNorm=keepNorm)
            # wVar['brD2420_DstPiPiUp'], wVar['brD2420_DstPiPiDown'] = wD1Up, wD1Dw
            # weights['DstPiPi_test'] *= w
            # w, wD2Up, wD2Dw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 415, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=3.0, keepNorm=keepNorm)
            # wVar['brD2460_DstPiPiUp'], wVar['brD2460_DstPiPiDown'] = wD2Up, wD2Dw
            # weights['DstPiPi_test'] *= w
            #
            # w, wUp, wDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.5, keepNorm=keepNorm)
            # wVar['brD2430_DstPiPiUp'], wVar['brD2430_DstPiPiDown'] = wUp, wDw
            # weights['DstPiPi_test'] *= w
            # w, wUp, wDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 100413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.5, keepNorm=keepNorm)
            # wVar['brDst2S_DstPiPiUp'], wVar['brDst2S_DstPiPiDown'] = wUp, wDw
            # weights['DstPiPi_test'] *= w


        ############################
        # Hc mix variations
        ############################
        # From https://www.overleaf.com/read/ykppfynnfxdt
        inflateRate = 1.5
        if n == 'Bd_DstDu': #1
            print 'Including Bd->D*Du br variations'
            # 1.1
            _, wVar['brBd_DstDuKUp'], wVar['brBd_DstDuKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421,
                                                                                         'MC_StrangeDstSisPdgId': 321},
                                                                                    inflateRate*0.21/2.47)
            # 1.3, 1.5-8
            _, x1_3u, x1_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 321}, inflateRate*0.9/10.6)
            _, x1_5u, x1_5d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 311}, inflateRate*0.47/5.4)
            wVar['brBd_DstDustKUp'] = x1_3u * x1_5u
            wVar['brBd_DstDustKDown'] = x1_3d * x1_5d

            # 1.7-10
            _, wVar['brBd_DstDustUp'], wVar['brBd_DstDustDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                           'MC_StrangeDstSisPdgId': 0},
                                                                                      inflateRate*0.4/5.4)

            # Correlate all K* decays (1.2, 1.4, 1.6-9)
            x1_2, x1_2u, x1_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            x1_4, x1_4u, x1_4d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            x1_6, x1_6u, x1_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            weights['brBd_DstDuKst'] = x1_2 * x1_4 * x1_6
            wVar['brBd_DstDuKstUp'] = x1_2u * x1_4u * x1_6u
            wVar['brBd_DstDuKstDown'] = x1_2d * x1_4d * x1_6d
        if n == 'Bd_DstDd': #2
            # 2.1-2
            _, wVar['brBd_DstDdKUp'], wVar['brBd_DstDdKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                         'MC_StrangeDstSisPdgId': 311},
                                                                                    inflateRate*0.25/3.2)
            # 2.5-8
            _, wVar['brBd_DstDdstKUp'], wVar['brBd_DstDdstKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                             'MC_StrangeDstSisPdgId': 311},
                                                                                        inflateRate*0.23/2.7)
            # 2.7-10
            _, wVar['brBd_DstDdstUp'], wVar['brBd_DstDdstDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                           'MC_StrangeDstSisPdgId': 0},
                                                                                      inflateRate*0.2/2.6)
            # 2.11, 2.12
            _, wVar['brBd_DstDdUp'], wVar['brBd_DstDdDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                       'MC_StrangeDstSisPdgId': 0},
                                                                                  inflateRate*1.5/6.1)
            # Correlate all K* decays (2.3-4, 2.6-9)
            x2_3, x2_3u, x2_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x2_6, x2_6u, x2_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            weights['brBd_DstDdKst'] = x2_3 * x2_6
            wVar['brBd_DstDdKstUp'] = x2_3u * x2_6u
            wVar['brBd_DstDdKstDown'] = x2_3u * x2_6u
        if n == 'Bd_DstDs': #3
            print 'Including Bd->D*Ds br variations'
            # 3.1
            _, wVar['brBd_DstDsUp'], wVar['brBd_DstDsDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 431}, inflateRate*1.1/8.0)
            # 3.2
            _, wVar['brBd_DstDsstUp'], wVar['brBd_DstDsstDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 433}, inflateRate*0.14/17.7)
            # 3.3
            _, wVar['brBd_DstDsst0Up'], wVar['brBd_DstDsst0Down'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 10431}, inflateRate*0.7/1.5)
        if n == 'Bu_DstDu': #4
            # 4.1
            _, wVar['brBu_DstDuKUp'], wVar['brBu_DstDuKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421,
                                                                                         'MC_StrangeDstSisPdgId': 311},
                                                                                    inflateRate*0.4/3.8)
            # 4.3, 4.5-7
            _, x4_3u, x4_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 311}, inflateRate*1.2/9.2)
            _, x4_5u, x4_5d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 321}, inflateRate*1.2/8.8)
            wVar['brBu_DstDustKUp'] = x4_3u * x4_5u
            wVar['brBu_DstDustKDown'] = x4_3d * x4_5d

            # Correlate all K* decays (4.2, 4.4, 4.6-8)
            x4_2, x4_2u, x4_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x4_4, x4_4u, x4_4d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x4_6, x4_6u, x4_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            weights['brBu_DstDuKst'] = x4_2 * x4_4 * x4_6
            wVar['brBu_DstDuKstUp'] = x4_2u * x4_4u * x4_6u
            wVar['brBu_DstDuKstDown'] = x4_2d * x4_4d * x4_6d

            # 4.9, 4.10
            _, x4_9u, x4_9d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 0}, inflateRate*0.5/3.9)
            _, x4_10u, x4_10d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 0}, inflateRate*1.7/8.1)
            wVar['brBu_DstDuUp'] = x4_9u * x4_10u
            wVar['brBu_DstDuDown'] = x4_9d * x4_10d
        if n == 'Bu_DstDd': #5
            print 'Including Bu->D*Dd br variations'
            # 5.1-3
            _, wVar['brBu_DstDdKUp'], wVar['brBu_DstDdKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                         'MC_StrangeDstSisPdgId': 321},
                                                                                    inflateRate*1.2/6.0)
            # 5.5-7
            _, wVar['brBu_DstDdstKUp'], wVar['brBu_DstDdstKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                             'MC_StrangeDstSisPdgId': 321},
                                                                                        inflateRate*0.6/4.4)
            # Correlate all K* decays (5.2-4, 5.6-8)
            x5_2, x5_2u, x5_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.0)
            x5_6, x5_6u, x5_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.0)
            weights['brBu_DstDdKst'] = x5_2 * x5_6
            wVar['brBu_DstDdKstUp'] = x5_2u * x5_6u
            wVar['brBu_DstDdKstDown'] = x5_2d * x5_6d
        if n == 'Bs_DstDs': #6
            x6_2_4, x6_2_4u, x6_2_4d = computeBrVarWeights(ds, {'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.0)
            weights['brBs_DstDsKst'] = x6_2_4
            wVar['brBs_DstDsKstUp'] = x6_2_4u
            wVar['brBs_DstDsKstDown'] = x6_2_4d

        print 'Computing total weights'
        weightsCentral = np.ones_like(ds['q2'])
        for wName, w in weights.iteritems():
            if np.max(np.abs(w)) > 3:
                iMax = np.argmax(np.abs(w))
                print '[WARNING] Max weights', wName, ': {:.3}'.format(w[iMax])
            weightsCentral *= w
        print 'N tot expected (after weights): {:.3f}k'.format(1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected)
        totalCounting[0] += 1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected
        evCountStr = '{:.2f} ({:.2f})'.format(1e-3*nTotExp*np.sum(weightsCentral)/nTotSelected, 1e-3*nTotSelected)
        eventCountingStr[n] = evCountStr
        wVar[''] = np.ones_like(weightsCentral)

        if args.dumpWeightsTree and not 'data' in n:
            weightsDir = os.path.join(MCsample[n].skimmed_dir, 'weights')
            if not os.path.isdir(weightsDir):
                os.makedirs(weightsDir)

            mcType = 'bare' if args.bareMC else 'corr'
            auxName = category.name + '_' + mcType + '_' + card_name + '.root'

            wDf = pd.DataFrame.from_dict({'central': weightsCentral*nTotExp/nTotSelected})
            for k, w in wVar.iteritems():
                if k:
                    wDf[k] = w*wDf['central']
            print 'Dumping weights tree in', os.path.join(weightsDir, auxName)
            rtnp.array2root(wDf.to_records(), os.path.join(weightsDir, auxName), treename='weights', mode='RECREATE')


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
            for var in ['M2_miss', 'Est_mu', 'U_miss', 'mu_pt', 'Dst_pt', 'K_pt', 'pi_pt', 'pis_pt', 'mass_D0pismu']:
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
        for var in additionalVariablesToPlot:
            if not var in histo.keys():
                histo[var] = {}
            for name_wVar, v_wVar in wVar.iteritems():
                h_name = n
                if not name_wVar == '':
                    h_name += '__' + name_wVar
                w = weightsCentral*v_wVar
                scale = nTotExp/nTotSelected
                varName = var
                if var == 'specQ2':
                    varName = 'q2'
                histo[var][h_name] = create_TH1D(ds[varName], name=h_name, weights=w, scale_histo=scale,
                                                    binning=binning[var], opt='underflow,overflow')

    evCountStr = '{:.1f} ({:.1f})'.format(*totalCounting)
    eventCountingStr['tot'] = evCountStr

    # Do the unrolling
    print '\n\n########### Unrolling 2D histograms ###########'
    unrolledBins = []
    unrollingCutOff = 3
    for i_q2 in range(len(binning['q2'])-1):
        unrolledBins.append([])
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
                    unrolledBins[i_q2].append([ix, iy])
                else:
                    nDroppedBins += 1
                    nExpectedDroppedEvents += hSum.GetBinContent(ix, iy)
        print 'Dropped bins:', nDroppedBins
        print 'Expected dropped candidates:', nExpectedDroppedEvents

        nameU = 'Unrolled_q2bin'+str(i_q2)
        histo[nameU] = {}
        validBins = unrolledBins[i_q2]
        for n, h in histo[name2D].iteritems():
            hUnrolled = rt.TH1D(h.GetName(), h.GetTitle(), len(validBins), 0.5, len(validBins)+0.5)
            for i, (ix, iy) in enumerate(validBins):
                hUnrolled.SetBinContent(i+1, h.GetBinContent(ix, iy))
                hUnrolled.SetBinError(i+1, h.GetBinError(ix, iy))
            histo[nameU][n] = hUnrolled
    if not args.dumpWeightsTree:
        pickle.dump(unrolledBins, open(outdir+'/unrolledBinsMap.pkl', 'wb'))


    ######################################################
    ########## Control region
    ######################################################
    ctrlVar = {}
    ctrlVar['ctrl_p__mHad'] = 'massHadTks'
    binning['ctrl_p__mHad'] = [35, 2.13, 2.83]

    ctrlVar['ctrl_p_mVis'] = 'massVisTks'
    binning['ctrl_p_mVis'] = array('d', [2.8] + list(np.arange(3., 5.5, 0.12)) + [5.5] )

    ctrlVar['ctrl_m__mHad'] = 'massHadTks'
    binning['ctrl_m__mHad'] = [40, 2.1, 3.3]

    ctrlVar['ctrl_pm_mVis'] = 'massVisTks'
    # binning['ctrl_pm_mVis'] = array('d', [2.8] + list(np.arange(3., 5.5, 0.07)) + [5.55] )
    binning['ctrl_pm_mVis'] = array('d', [2.8] + list(np.arange(3., 5.5, 0.12)) + [5.5] )

    ctrlVar['ctrl_pm_mHad'] = 'massHadTks'
    # binning['ctrl_pm_mHad'] = [75, 2.3, 3.8]
    binning['ctrl_pm_mHad'] = [50, 2.3, 3.8]

    ctrlVar['ctrl_mm_mHad'] = 'massHadTks'
    binning['ctrl_mm_mHad'] = [20, 2.25, 3.6]

    ctrlVar['ctrl_pp_mHad'] = 'massHadTks'
    binning['ctrl_pp_mHad'] = [20, 2.25, 3.7]

    ctrlVar['ctrl_m__tk_pt_0'] = 'tkPt_0'
    binning['ctrl_m__tk_pt_0'] = array('d', list(np.arange(.5, 2., 0.1)) + list(np.arange(2., 7.51, 0.25)) )
    ctrlVar['ctrl_mm_tk_pt_0'] = 'tkPt_0'
    binning['ctrl_mm_tk_pt_0'] = [12, 0.5, 10]
    ctrlVar['ctrl_mm_tk_pt_1'] = 'tkPt_1'
    binning['ctrl_mm_tk_pt_1'] = [20, 0.5, 4]
    ctrlVar['ctrl_p__tk_pt_0'] = 'tkPt_0'
    binning['ctrl_p__tk_pt_0'] = [50, 0.5, 15]
    ctrlVar['ctrl_pm_tk_pt_0'] = 'tkPt_0'
    binning['ctrl_pm_tk_pt_0'] = [40, 0.5, 18]
    ctrlVar['ctrl_pm_tk_pt_1'] = 'tkPt_1'
    binning['ctrl_pm_tk_pt_1'] = [30, 0.5, 8]
    ctrlVar['ctrl_pp_tk_pt_0'] = 'tkPt_0'
    binning['ctrl_pp_tk_pt_0'] = [30, 0.5, 15]
    ctrlVar['ctrl_pp_tk_pt_1'] = 'tkPt_1'
    binning['ctrl_pp_tk_pt_1'] = [30, 0.5, 4]


    for s in ['m_', 'p_', 'mm', 'pm', 'pp']:
        ctrlVar['ctrl_'+s+'_mu_pt'] = 'mu_pt'
        binning['ctrl_'+s+'_mu_pt'] = binning['mu_pt'][0]

        ctrlVar['ctrl_'+s+'_umiss'] = 'UmissTks'
        binning['ctrl_'+s+'_umiss'] = [40, -0.08, 0.15]

        ctrlVar['ctrl_'+s+'_M2miss'] = 'M2_miss'
        binning['ctrl_'+s+'_M2miss'] = [25, -2, 8]

        ctrlVar['ctrl_'+s+'_q2'] = 'M2_miss'
        binning['ctrl_'+s+'_q2'] = [50, -2, 15]

        if not s[1] == '_':
            ctrlVar['ctrl_'+s+'_mPiPi'] = 'massTks_pipi'
            binning['ctrl_'+s+'_mPiPi'] = [30, 0.25, 1.]

            ctrlVar['ctrl_'+s+'_mDstPi_0'] = 'tkMassHad_0'
            binning['ctrl_'+s+'_mDstPi_0'] = [30, 2.1, 3.]

    # Figuring out the mod to avoid double counting
    ctrlVar_mod = defaultdict(lambda : None)
    ctrlVar_counter = defaultdict(lambda : 0)
    for r in args.controlRegions:
        if not 'ctrl_'+r in ctrlVar.keys():
            print '[ERROR] Region "{}" not defined'.format(r)
            raise

        ctrlVar_mod['ctrl_'+r] = [ -1, ctrlVar_counter[r[:2]] ]
        ctrlVar_counter[r[:2]]+= 1

    for k in ctrlVar_mod.keys():
        ctrlVar_mod[k][0] = ctrlVar_counter[k[5:7]]

    print '[DEBUG] Control regions modulos:', ','.join([k+':'+str(v) for k, v in ctrlVar_mod.iteritems()])


    print '---------> Fill control regions histograms'
    for k in ctrlVar.keys():
        histo[k] = {}

    totalCounting = {}
    for n in processOrder:
        ds = dSetTkSide[n]
        if n == 'data': continue
        print '\n----------->', n, '<-------------'
        wVar = {}
        weights = {}
        if n == 'dataSS_DstMu':
            nTotExp = ds['q2'].shape[0]
        else:
            sMC = MCsample[n]

            nGenExp = sMC.effMCgen['xsec'][0] * expectedLumi[category.name] * data_to_MC_prescale
            eff = [1, 0]
            for f, df in [sMC.effMCgen['effGEN'],
                          decayBR[n],
                          sMC.effCand['effCAND'],
                          sMC.getSkimEff(category.name+'_trkCtrl_'+mcType),
                         ]:
                eff[0] *= f
                eff[1] += np.square(df/f)
            eff[1] = eff[0] * np.sqrt(eff[1])
            if n+'_tk' in corrScaleFactors.keys():
                eff[0] *= corrScaleFactors[n+'_tk']
                print 'Using scale factor from a posteriori selection: {:.3f}'.format(corrScaleFactors[n+'_tk'])
            nTotExp = nGenExp*eff[0]

            print 'Including pileup reweighting'
            weights['pileup'] = puReweighter.getPileupWeights(ds['MC_nInteractions'])

            print 'Including beam spot correction'
            weights['beamSpot'] = getBeamSpotCorrectionWeights(ds, beamSpotParam)

            print 'Including trigger corrections'
            nameSF = 'trg{}SF'.format(category.trg)
            weights[nameSF], wSfUp, wSfDw = computeTrgSF(ds, hTriggerSF)
            auxOnes = np.ones_like(wSfUp)
            # for i_eta in range(1, hTriggerSF.GetNbinsZ()+1):
            #     c_eta = hTriggerSF.GetZaxis().GetBinCenter(i_eta)
            #     w_eta = hTriggerSF.GetZaxis().GetBinWidth(i_eta)
            # for i_ip in range(1, hTriggerSF.GetNbinsY()+1):
            #     c_ip = hTriggerSF.GetYaxis().GetBinCenter(i_ip)
            #     w_ip = hTriggerSF.GetYaxis().GetBinWidth(i_ip)
            #     if c_ip + 0.5*w_ip <= category.minIP:
            #         continue
            for i_pt in range(1, hTriggerSF.GetNbinsX()+2):
                if i_pt > hTriggerSF.GetNbinsX() and category.name == 'High':
                    sel = ds['mu_pt'] > hTriggerSF.GetXaxis().GetXmax()
                else:
                    c_pt = hTriggerSF.GetXaxis().GetBinCenter(i_pt)
                    w_pt = hTriggerSF.GetXaxis().GetBinWidth(i_pt)
                    if (c_pt + 0.5*w_pt <= category.min_pt) or (c_pt - 0.5*w_pt >= category.max_pt):
                        continue

                    sel = np.abs(ds['mu_pt'] - c_pt) < w_pt
                # sel = np.logical_and(sel, np.abs(ds['mu_sigdxy'] - c_ip) < w_ip)
                # sel = np.logical_and(sel, np.abs(ds['mu_eta'] - c_eta) < w_eta)
                # binName = '_pt{}ip{}eta{}'.format(i_pt, i_ip, i_eta)
                # print 'Trg SF', i_pt, i_ip, i_eta, '-> selected {}'.format(np.sum(sel))
                # binName = '_pt{}ip{}'.format(i_pt, i_ip)
                binName = '_pt{}'.format(i_pt)
                wVar[nameSF+binName+'Up'] = np.where(sel, wSfUp, auxOnes)
                wVar[nameSF+binName+'Down'] = np.where(sel, wSfDw, auxOnes)

            print 'Including muon ID corrections'
            weights['muonIdSF'], _, _ = computeMuonIDSF(ds)

            print 'Including soft track pT corrections'
            # l = ['K_pt', 'pi_pt', 'pis_pt', 'tkPt_0', 'tkPt_1']
            # weights['softTrkEff'] = weightsSoftTrackEff(ds, l, parsSoftTracks['w'][0], parsSoftTracks['s'][0])
            # for mod in [+1, -1]:
            #     varName = 'Up' if mod > 0 else 'Down'
            #     w = parsSoftTracks['w'][0] + mod*parsSoftTracks['w'][1]
            #     wVar['softTrkEff_w'+varName] = weightsSoftTrackEff(ds, l, w, parsSoftTracks['s'][0])/weights['softTrkEff']
            #     s = parsSoftTracks['s'][0] + mod*parsSoftTracks['s'][1]
            #     wVar['softTrkEff_s'+varName] = weightsSoftTrackEff(ds, l, parsSoftTracks['w'][0], s)/weights['softTrkEff']
            partList = ['K_pt', 'pi_pt', 'pis_pt', 'tkPt_0', 'tkPt_1']
            for nBin in range(len(softPtUnc)):
                refPt = '{:.0f}'.format(np.round(np.mean(softPtUnc[nBin][:-1])*1e3))
                wVar['softTrkEff_'+refPt+'Up'] = weightsSoftTrackEff(ds, partList, bin=nBin, size=+1)
                wVar['softTrkEff_'+refPt+'Down'] = weightsSoftTrackEff(ds, partList, bin=nBin, size=-1)

            ############################
            # B phase space corrections
            ############################
            print 'Including B eta corrections'
            cname = 'B_eta'+category.name
            weights[cname], auxVarDic = computeKinCalWeights(ds, 'MC_B_eta', cname, cal_eta_B)
            wVar.update(auxVarDic)

            if (not args.calBpT == 'none') and (n in samples_Bd):
                print 'Including Bd pT corrections'
                cname = 'BdpT'+category.name
                if cal_pT_Bd.kind == 'ratio':
                    weights[cname], wVar[cname+'Up'], wVar[cname+'Down'] = computeKinCalWeights(ds, 'MC_B_pt', None, cal_pT_Bd)
                else:
                    weights[cname], auxVarDic = computeKinCalWeights(ds, 'MC_B_pt', cname, cal_pT_Bd)
                    wVar.update(auxVarDic)

            # Correct the amount of random tracks from PV
            aux = '' if args.correlate_tkPVfrac else category.name
            weights['tkPVfrac'], wVar['tkPVfrac'+aux+'Up'], wVar['tkPVfrac'+aux+'Down'] = computeTksPVweights(ds, relScale=0.3, centralVal=1.4)
            # weights['tkPVfrac'], wVar['tkPVfrac'+category.name+'Up'], wVar['tkPVfrac'+category.name+'Down'] = computeTksPVweights(ds, relScale=0.05, centralVal=2.3)
            print 'Average tkPVfrac weight: {:.2f}'.format(np.mean(weights['tkPVfrac']))

        ############################
        # Form factor correction
        ############################
        if n in ['mu', 'tau'] and schemeFF != 'NoFF':
            print 'Including FF corrections (Hammer)'
            weights['B2DstFF'] = ds['wh_'+schemeFF+'Central']*sMC.effCand['rate_den']/sMC.effCand['rate_'+schemeFF+'Central']
            for nPar in FreeParFF:
                for var in ['Up', 'Down']:
                    tag = schemeFF + nPar + var
                    wVar['B2Dst'+tag] = ds['wh_'+tag]/ds['wh_'+schemeFF+'Central']
                    wVar['B2Dst'+tag] *= sMC.effCand['rate_'+schemeFF+'Central']/sMC.effCand['rate_' + tag]
                    wVar['B2Dst'+tag][np.isnan(wVar['B2Dst'+tag])] = 0.


        if '_MuDstPi' in n:
            print 'Including B -> D**MuNu FF corrections (Hammer)'
            weights['purgeUwantedEvts'] = np.array(ds['DststProc_id'] >= 0).astype(np.int)
            DststProc_id = np.array(ds['DststProc_id']).astype(np.int)

            ratesRatio = np.ones(50)
            ratesRatio[1] = sMC.effCand['rate_den_D1']/sMC.effCand['rate_D1BLR_central']
            ratesRatio[2] = sMC.effCand['rate_den_D1st']/sMC.effCand['rate_D1stBLR_central']
            ratesRatio[3] = sMC.effCand['rate_den_D2st']/sMC.effCand['rate_D2stBLR_central']
            selDstst = np.logical_and(DststProc_id >= 1, DststProc_id <= 3)
            weights['B2DststFF'] = np.where(selDstst, ds['wh_Dstst_BLRCentral']*ratesRatio[DststProc_id], 1.)

            selDstst = np.logical_or(DststProc_id == 1, DststProc_id == 3)
            for nEig in [1,2,3,4]:
                for var in ['Up', 'Down']:
                    tag = 'DststN_BLReig'+str(nEig)+var
                    ratesRatio = np.ones(50)
                    ratesRatio[1] = sMC.effCand['rate_D1BLR_central']/sMC.effCand['rate_D1BLR_eig'+str(nEig)+var]
                    ratesRatio[3] = sMC.effCand['rate_D2stBLR_central']/sMC.effCand['rate_D2stBLR_eig'+str(nEig)+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_Dstst_BLRCentral'])*ratesRatio[DststProc_id], 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.

            selDstst = DststProc_id == 2
            for nEig in [1,2,3]:
                for var in ['Up', 'Down']:
                    tag = 'DststW_BLReig'+str(nEig)+var
                    ratesRatio = sMC.effCand['rate_D1stBLR_central']/sMC.effCand['rate_D1stBLR_eig'+str(nEig)+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_Dstst_BLRCentral'])*ratesRatio, 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.

        if n.endswith('_MuDstPiPi'):
            print 'Including B -> D(2S)MuNu FF corrections (Hammer)'
            DststProc_id = np.array(ds['DststProc_id']).astype(np.int)

            ratesRatio = np.ones(50)
            ratesRatio[4] = sMC.effCand['rate_den_D2sst']/sMC.effCand['rate_D2sstBLOP_central']
            ratesRatio[11] = ratesRatio[4]
            ratesRatio[12] = ratesRatio[4]
            ratesRatio[13] = ratesRatio[4]
            ratesRatio[23] = sMC.effCand['rate_den_D2s']/sMC.effCand['rate_D2sBLOP_central']

            selDstst = np.logical_and(DststProc_id >= 4, DststProc_id <= 23)
            weights['B2D2S_FF'] = np.where(selDstst, ds['wh_D2S_BLOPCentral']*ratesRatio[DststProc_id], 1.)

            for parName in ['RhoSq', 'chi11', 'chi21', 'chi31', 'eta1']:
                for var in ['Up', 'Down']:
                    tag = 'D2S_BLOP'+parName+var
                    ratesRatio = np.ones(50)
                    ratesRatio[4] = sMC.effCand['rate_D2sstBLOP_central']/sMC.effCand['rate_D2sstBLOP_'+parName+var]
                    ratesRatio[11] = ratesRatio[4]
                    ratesRatio[12] = ratesRatio[4]
                    ratesRatio[13] = ratesRatio[4]
                    ratesRatio[23] = sMC.effCand['rate_D2sBLOP_central']/sMC.effCand['rate_D2sBLOP_'+parName+var]
                    wVar['FF_B2'+tag] = np.where(selDstst, (ds['wh_'+tag]/ds['wh_D2S_BLOPCentral'])*ratesRatio[DststProc_id], 1.)
                    wVar['FF_B2'+tag][np.isnan(wVar['FF_B2'+tag])] = 0.

        ############################
        # Dstst resonance mix
        ############################
        if not re.search('DstPi\Z', n) is None:
            print 'Including D**->D*Pi branching ratio and width variations'

            infate = 2
            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10423}, infate*0.2/3.0)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10413}, infate*0.14/1.40)
            wVar['brB_D2420MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2420MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20423}, infate*0.6/2.7)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20413}, infate*0.45/1.55)
            wVar['brB_D2430MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2430MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 425}, infate*0.24/1.01)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 415}, infate*0.06/0.34)
            wVar['brB_D2460MuNuUp'] = wNeuUp * wChUp
            wVar['brB_D2460MuNuDown'] = wNeuDw * wChDw

            _, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_DstMotherPdgId': 521}, infate*0.5)
            _, wChUp, wChDw = computeBrVarWeights(ds, {'MC_DstMotherPdgId': 511}, infate*0.5)
            wVar['brB_DstPiMuNuUp'] = wNeuUp * wChUp
            wVar['brB_DstPiMuNuDown'] = wNeuDw * wChDw


            w, wVar['D2420_widthUp'], wVar['D2420_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[10423, 2.421, 0.0274],
                                                                                                 [10413, 2.423, 0.020]],
                                                                                       newGamma=[0.030, 0.025],
                                                                                       relScale=0.15)
            weights['D2420_width'] = w

            w, wVar['D2430_widthUp'], wVar['D2430_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[20423, 2.445, 0.250],
                                                                                                 [20413, 2.445, 0.250]],
                                                                                       newGamma=[0.3, 0.3],
                                                                                       relScale=0.1)
            weights['D2430_width'] = w

            w, wVar['D2460_widthUp'], wVar['D2460_widthDown'] = computeWidthVarWeights(ds,
                                                                                       selItems=[[425, 2.461, 0.049],
                                                                                                 [415, 2.460, 0.037]],
                                                                                       newGamma=[0.045, 0.045],
                                                                                       relScale=0.15)
            weights['D2460_width'] = w

        if not re.search('MuDstPiPi\Z', n) is None:
            print 'Including B -> MuNuD*PiPi width and Br variations'
            widthMods = [[x, 2.640, 0.300] for x in [100413, 100423, 100411, 100421]]
            uncN = 'Dst2S_width'
            weights[uncN], wVar[uncN+'Up'], wVar[uncN+'Down'] = computeWidthVarWeights(ds, selItems=widthMods, relScale=0.3)

            keepNorm=False
            nnn = 'brDstPiPi_NR'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 0}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 1}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1st'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 2}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2st'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 3}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 4}, relScale=0.5, centralVal=1., keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D1Pi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 11}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D1stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 12}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2Sst_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 13}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2S_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 23}, relScale=0.5, centralVal=0.5, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1Pi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 31}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D1stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 32}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            nnn = 'brDstPiPi_D2stPi'
            weights[nnn], wVar[nnn+'Up'], wVar[nnn+'Down'] = computeBrVarWeights(ds, {'DststProc_id': 33}, relScale=0.5, centralVal=1.0, keepNorm=keepNorm)

            # print 'Including D**->D*PiPi width variations'
            # widthMods = [[x, 2.640, 0.400] for x in [100413, 100423]]
            # uncN = 'Dst2S_width'
            # weights[uncN], wVar[uncN+'Up'], wVar[uncN+'Down'] = computeWidthVarWeights(ds, selItems=widthMods, relScale=0.5)
            #
            # keepNorm=False
            # wNeu, wNeuUp, wNeuDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_1': 111}, relScale=1.0, centralVal=2., keepNorm=keepNorm)
            # wCh, wChUp, wChDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_1': 211}, relScale=1.0, centralVal=2., keepNorm=keepNorm)
            # wVar['brDstst_DststPiUp'] = wNeuUp * wChUp
            # wVar['brDstst_DststPiDown'] = wNeuDw * wChDw
            # weights['DstPiPi_test'] = wNeu*wCh
            #
            # w, wD1Up, wD1Dw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 10413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.1, keepNorm=keepNorm)
            # wVar['brD2420_DstPiPiUp'], wVar['brD2420_DstPiPiDown'] = wD1Up, wD1Dw
            # weights['DstPiPi_test'] *= w
            # w, wD2Up, wD2Dw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 415, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=3.0, keepNorm=keepNorm)
            # wVar['brD2460_DstPiPiUp'], wVar['brD2460_DstPiPiDown'] = wD2Up, wD2Dw
            # weights['DstPiPi_test'] *= w
            #
            # w, wUp, wDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 20413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.5, keepNorm=keepNorm)
            # wVar['brD2430_DstPiPiUp'], wVar['brD2430_DstPiPiDown'] = wUp, wDw
            # weights['DstPiPi_test'] *= w
            # w, wUp, wDw = computeBrVarWeights(ds, {'MC_munuSisterPdgId_0': 100413, 'MC_munuSisterPdgId_1': 0}, relScale=1.0, centralVal=0.5, keepNorm=keepNorm)
            # wVar['brDst2S_DstPiPiUp'], wVar['brDst2S_DstPiPiDown'] = wUp, wDw
            # weights['DstPiPi_test'] *= w


        ############################
        # Hc mix variations
        ############################
        # From https://www.overleaf.com/read/ykppfynnfxdt
        inflateRate = 1.5
        if n == 'Bd_DstDu': #1
            print 'Including Bd->D*Du br variations'
            # 1.1
            _, wVar['brBd_DstDuKUp'], wVar['brBd_DstDuKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421,
                                                                                         'MC_StrangeDstSisPdgId': 321},
                                                                                    inflateRate*0.21/2.47)
            # 1.3, 1.5-8
            _, x1_3u, x1_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 321}, inflateRate*0.9/10.6)
            _, x1_5u, x1_5d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 311}, inflateRate*0.47/5.4)
            wVar['brBd_DstDustKUp'] = x1_3u * x1_5u
            wVar['brBd_DstDustKDown'] = x1_3d * x1_5d

            # 1.7-10
            _, wVar['brBd_DstDustUp'], wVar['brBd_DstDustDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                           'MC_StrangeDstSisPdgId': 0},
                                                                                      inflateRate*0.4/5.4)

            # Correlate all K* decays (1.2, 1.4, 1.6-9)
            x1_2, x1_2u, x1_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            x1_4, x1_4u, x1_4d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            x1_6, x1_6u, x1_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            weights['brBd_DstDuKst'] = x1_2 * x1_4 * x1_6
            wVar['brBd_DstDuKstUp'] = x1_2u * x1_4u * x1_6u
            wVar['brBd_DstDuKstDown'] = x1_2d * x1_4d * x1_6d
        if n == 'Bd_DstDd': #2
            # 2.1-2
            _, wVar['brBd_DstDdKUp'], wVar['brBd_DstDdKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                         'MC_StrangeDstSisPdgId': 311},
                                                                                    inflateRate*0.25/3.2)
            # 2.5-8
            _, wVar['brBd_DstDdstKUp'], wVar['brBd_DstDdstKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                             'MC_StrangeDstSisPdgId': 311},
                                                                                        inflateRate*0.23/2.7)
            # 2.7-10
            _, wVar['brBd_DstDdstUp'], wVar['brBd_DstDdstDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                           'MC_StrangeDstSisPdgId': 0},
                                                                                      inflateRate*0.2/2.6)
            # 2.11, 2.12
            _, wVar['brBd_DstDdUp'], wVar['brBd_DstDdDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                       'MC_StrangeDstSisPdgId': 0},
                                                                                  inflateRate*1.5/6.1)
            # Correlate all K* decays (2.3-4, 2.6-9)
            x2_3, x2_3u, x2_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x2_6, x2_6u, x2_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            weights['brBd_DstDdKst'] = x2_3 * x2_6
            wVar['brBd_DstDdKstUp'] = x2_3u * x2_6u
            wVar['brBd_DstDdKstDown'] = x2_3u * x2_6u
        if n == 'Bd_DstDs': #3
            print 'Including Bd->D*Ds br variations'
            # 3.1
            _, wVar['brBd_DstDsUp'], wVar['brBd_DstDsDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 431}, inflateRate*1.1/8.0)
            # 3.2
            _, wVar['brBd_DstDsstUp'], wVar['brBd_DstDsstDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 433}, inflateRate*0.14/17.7)
            # 3.3
            _, wVar['brBd_DstDsst0Up'], wVar['brBd_DstDsst0Down'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 10431}, inflateRate*0.7/1.5)
        if n == 'Bu_DstDu': #4
            # 4.1
            _, wVar['brBu_DstDuKUp'], wVar['brBu_DstDuKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421,
                                                                                         'MC_StrangeDstSisPdgId': 311},
                                                                                    inflateRate*0.4/3.8)
            # 4.3, 4.5-7
            _, x4_3u, x4_3d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 311}, inflateRate*1.2/9.2)
            _, x4_5u, x4_5d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 321}, inflateRate*1.2/8.8)
            wVar['brBu_DstDustKUp'] = x4_3u * x4_5u
            wVar['brBu_DstDustKDown'] = x4_3d * x4_5d

            # Correlate all K* decays (4.2, 4.4, 4.6-8)
            x4_2, x4_2u, x4_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x4_4, x4_4u, x4_4d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.)
            x4_6, x4_6u, x4_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.)
            weights['brBu_DstDuKst'] = x4_2 * x4_4 * x4_6
            wVar['brBu_DstDuKstUp'] = x4_2u * x4_4u * x4_6u
            wVar['brBu_DstDuKstDown'] = x4_2d * x4_4d * x4_6d

            # 4.9, 4.10
            _, x4_9u, x4_9d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 421, 'MC_StrangeDstSisPdgId': 0}, inflateRate*0.5/3.9)
            _, x4_10u, x4_10d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 423, 'MC_StrangeDstSisPdgId': 0}, inflateRate*1.7/8.1)
            wVar['brBu_DstDuUp'] = x4_9u * x4_10u
            wVar['brBu_DstDuDown'] = x4_9d * x4_10d
        if n == 'Bu_DstDd': #5
            print 'Including Bu->D*Dd br variations'
            # 5.1-3
            _, wVar['brBu_DstDdKUp'], wVar['brBu_DstDdKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411,
                                                                                         'MC_StrangeDstSisPdgId': 321},
                                                                                    inflateRate*1.2/6.0)
            # 5.5-7
            _, wVar['brBu_DstDdstKUp'], wVar['brBu_DstDdstKDown'] = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413,
                                                                                             'MC_StrangeDstSisPdgId': 321},
                                                                                        inflateRate*0.6/4.4)
            # Correlate all K* decays (5.2-4, 5.6-8)
            x5_2, x5_2u, x5_2d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 411, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.0)
            x5_6, x5_6u, x5_6d = computeBrVarWeights(ds, {'MC_CharmedDstSisPdgId': 413, 'MC_StrangeDstSisPdgId': 323}, 0.5, centralVal=4.0)
            weights['brBu_DstDdKst'] = x5_2 * x5_6
            wVar['brBu_DstDdKstUp'] = x5_2u * x5_6u
            wVar['brBu_DstDdKstDown'] = x5_2d * x5_6d
        if n == 'Bs_DstDs': #6
            x6_2_4, x6_2_4u, x6_2_4d = computeBrVarWeights(ds, {'MC_StrangeDstSisPdgId': 313}, 0.5, centralVal=4.0)
            weights['brBs_DstDsKst'] = x6_2_4
            wVar['brBs_DstDsKstUp'] = x6_2_4u
            wVar['brBs_DstDsKstDown'] = x6_2_4d

        print 'Computing total weights'
        weightsCentral = np.ones_like(ds['q2'])
        for w in weights.values():
            weightsCentral *= w
        wVar[''] = np.ones_like(weightsCentral)

        if args.dumpWeightsTree and not 'data' in n:
            weightsDir = os.path.join(MCsample[n].skimmed_dir, 'weights')
            if not os.path.isdir(weightsDir):
                os.makedirs(weightsDir)

            mcType = 'bare' if args.bareMC else 'corr'
            auxName = category.name + 'trkCtrl_' + mcType + '_' + card_name + '.root'

            wDf = pd.DataFrame.from_dict({'central': weightsCentral*nTotExp/float(weightsCentral.shape[0]) })
            for k, w in wVar.iteritems():
                if k:
                    wDf[k] = w*wDf['central']

            print 'Dumping weights tree in', os.path.join(weightsDir, auxName)
            rtnp.array2root(wDf.to_records(), os.path.join(weightsDir, auxName), treename='weights', mode='RECREATE')

        sel = {}
        scale = {}
        latexTableString = {}
        table = PrettyTable()
        table.field_names = ['Region', 'Selected', 'Exp. bare', 'Exp. weights']
        for k, selFun in controlRegSel.iteritems():
            sel[k] = selFun(ds)
            nTotSel = float(np.sum(sel[k]))
            nExp = nTotExp * nTotSel / sel[k].shape[0]
            if n == 'dataSS_DstMu' and k == 'mm':
                nExp *= 1e-3
            nAux = nTotExp * np.sum(weightsCentral[sel[k]]) / sel[k].shape[0]
            table.add_row([k] + '{:.0f} {:.0f} {:.0f}'.format(nTotSel, nExp, nAux).split(' '))
            latexTableString[k] = '{:.0f} ({:.0f})'.format(nAux, nTotSel)
            if not k in totalCounting.keys():
                totalCounting[k] = [0, 0]
            totalCounting[k][0] += nAux
            totalCounting[k][1] += nTotSel
            if nTotSel ==0:
                nTotSel+=1
            scale[k] = nExp/nTotSel
        print table.get_string(sortby="Region")
        s = ' & '.join([latexTableString[s] for s in ['m_', 'p_', 'mm', 'pm', 'pp']])
        print s
        eventCountingStr[n] += ' & ' + s + '\\\\'


        for name_wVar, v_wVar in wVar.iteritems():
            h_name = n
            if not name_wVar == '':
                h_name += '__' + name_wVar
            w = weightsCentral*v_wVar

            for k, var in ctrlVar.iteritems():
                region = k[5:7]

                auxSel = sel[region]
                if not ctrlVar_mod[k] is None:
                    m, j = ctrlVar_mod[k]
                    if m > 1:
                        auxSel = np.logical_and( np.mod(ds['index'], m) == j, auxSel )

                histo[k][h_name] = create_TH1D(ds[var][auxSel],
                                               name=h_name, title=h_name,
                                               binning=binning[k],
                                               opt='',
                                               weights=w[auxSel], scale_histo=scale[region]
                                              )


    if args.dumpWeightsTree:
        print 'Exiting runCombine.py. To run further remove dumpWeightsTree option.'
        exit()

    s = ' & '.join(['{:.0f} ({:.0f})'.format(*totalCounting[s]) for s in ['p_', 'm_', 'mm', 'pm', 'pp']]) + ' \\\\'
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
            for var in ['M2_miss', 'Est_mu', 'U_miss', 'mu_pt', 'Dst_pt', 'K_pt', 'pi_pt', 'pis_pt', 'mass_D0pismu']:
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
                validBins = unrolledBins[i_q2]
                hUnrolled = rt.TH1D('data_obs', 'Data Obs', len(validBins), 0.5, len(validBins)+0.5)
                for i, (ix, iy) in enumerate(validBins):
                    hUnrolled.SetBinContent(i+1, histo[cat2D]['data'].GetBinContent(ix, iy))
                    hUnrolled.SetBinError(i+1, histo[cat2D]['data'].GetBinError(ix, iy))
                histo[catU]['data'] = hUnrolled

        for var in additionalVariablesToPlot:
            varName = var
            if var == 'specQ2':
                varName = 'q2'
            histo[var]['data'] = create_TH1D(ds[varName], name='data_obs', binning=binning[var], opt='underflow,overflow')

        ds = dSetTkSide['data']
        for k, var in ctrlVar.iteritems():
            region = k[5:7]
            auxSel = controlRegSel[region](ds)
            if not ctrlVar_mod[k] is None:
                m, j = ctrlVar_mod[k]
                if m > 1:
                    auxSel = np.logical_and( np.mod(ds['index'], m) == j, auxSel )

            histo[k]['data'] = create_TH1D(ds[var][auxSel],
                                           name='data_obs', title='Data Obs',
                                           binning=binning[k],
                                           opt='',
                                          )

        print 'N observed data control regions: ' + ' & '.join(['{:.0f}'.format(histo['ctrl_'+s+'_mHad']['data'].Integral()) for s in ['p_', 'm_', 'mm', 'pm', 'pp']]) + ' \\\\'

    ######################################################
    ########## Dump root file
    ######################################################
    print 'Dumping histos in root files'
    if not os.path.isdir(histo_file_dir): os.makedirs(histo_file_dir)

    for cat_name, h_dic in histo.iteritems():
        tf = rt.TFile(histo_file_dir+'{}_{}.root'.format(card_name, cat_name), 'recreate')
        for v in h_dic.values():
            v.Write()
        tf.Close()


########################### -------- Create histrograms ------------------ #########################

def drawPlots(tag, hDic, catName, scale_dic={}):
    print 20*'-', 'Drawing plots', tag, 20*'-'
    outCanvas = []
    CMS_lumi.integrated_lumi = expectedLumi[catName]

    # Draw unrolled histograms
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
                                   addText='Cat. '+catName+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=True, legBkg=True,
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
    cAux = plot_gridVarQ2(CMS_lumi, binning, hDic, draw_pulls=True,
                          pullsRatio=True, pulls_ylim=[0.85, 1.1],
                          scale_dic=scale_dic,
                          categoryText=catName, cNameTag=tag,
                          iq2_maskData=[] if args.unblinded else [2, 3])
    cAux.SaveAs(outdir+'/fig/signalRegion_'+tag+'.png')
    cAux.SaveAs(webFolder+'/signalRegion_'+tag+'.png')
    outCanvas.append(cAux)

    if 'MVA' in hDic.keys():
        print 'Creating MVA'
        cAux = plot_SingleCategory(CMS_lumi, hDic['MVA'], draw_pulls=True, scale_dic=scale_dic,
                                   addText='Cat. '+catName, logy=True, legBkg=True,
                                   min_y=1, tag=tag+'MVA', legLoc=[0.2, 0.1, 0.4, 0.4])
        cAux.SaveAs(outdir+'/fig/MVA_'+tag+'.png')
        cAux.SaveAs(webFolder+'/MVA_'+tag+'.png')
        outCanvas.append(cAux)

    if 'B_pt' in hDic.keys():
        print 'Creating B_pt'
        hDic['B_pt']['data'].GetXaxis().SetTitle('B p_{T} [GeV]')
        hDic['B_pt']['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['B_pt'], draw_pulls=True, pullsRatio=True, scale_dic=scale_dic,
                                   addText='Cat. '+catName, logy=False, legBkg=True,
                                   min_y=1, tag=tag+'B_pt', legLoc=[0.65, 0.4, 0.9, 0.75])
        cAux.SaveAs(outdir+'/fig/B_pt_'+tag+'.png')
        cAux.SaveAs(webFolder+'/B_pt_'+tag+'.png')
        outCanvas.append(cAux)

        hDic['B_pt']['data'].GetYaxis().SetTitle('Normalized events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['B_pt'], draw_pulls=True, pullsRatio=True, scale_dic=scale_dic,
                                   density=True,
                                   addText='Cat. '+catName, logy=False, legBkg=True,
                                   min_y=0, tag=tag+'B_pt', legLoc=[0.65, 0.4, 0.9, 0.75])
        cAux.SaveAs(outdir+'/fig/B_pt_norm_'+tag+'.png')
        cAux.SaveAs(webFolder+'/B_pt_norm_'+tag+'.png')
        outCanvas.append(cAux)

    if 'specQ2' in hDic.keys():
        print 'Creating q2 spectrum'
        hDic['specQ2']['data'].GetXaxis().SetTitle('q^{2} [GeV]')
        hDic['specQ2']['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['specQ2'], draw_pulls=True, pullsRatio=True, scale_dic=scale_dic,
                                   addText='Cat. '+catName, logy=False, legBkg=True,
                                   min_y=1, tag=tag+'specQ2', legLoc=[0.75, 0.4, 0.93, 0.75])
        cAux.SaveAs(outdir+'/fig/q2_'+tag+'.png')
        cAux.SaveAs(webFolder+'/q2_'+tag+'.png')
        outCanvas.append(cAux)

        hDic['specQ2']['data'].GetYaxis().SetTitle('Normalized events')
        cAux = plot_SingleCategory(CMS_lumi, hDic['specQ2'], draw_pulls=True, pullsRatio=True, scale_dic=scale_dic,
                                   density=True,
                                   addText='Cat. '+catName, logy=False, legBkg=True,
                                   min_y=0, tag=tag+'specQ2', legLoc=[0.75, 0.4, 0.93, 0.75])
        cAux.SaveAs(outdir+'/fig/q2_norm_'+tag+'.png')
        cAux.SaveAs(webFolder+'/q2_norm_'+tag+'.png')
        outCanvas.append(cAux)

    varsToPlot = ['B_eta', 'mu_eta', 'K_eta', 'pi_eta', 'pis_eta','mass_piK','deltaM_DstD']
    varsToPlot += [v for v in hDic.keys() if v.startswith('ctrl_')]
    for var in varsToPlot:
        if var in hDic.keys():
            print 'Creating '+var
            title = var[8:] if var.startswith('ctrl_') else var
            title = title.replace('_eta', ' #eta').replace('st', '*').replace('_', ' ')
            hDic[var]['data'].GetXaxis().SetTitle(title)
            hDic[var]['data'].GetYaxis().SetTitle('Events')
            ctrl_text = (', '+getControlSideText(var)) if var.startswith('ctrl_') else ''
            cAux = plot_SingleCategory(CMS_lumi, hDic[var],
                                       draw_pulls=True, pullsRatio=True, pulls_ylim='auto',
                                       scale_dic=scale_dic,
                                       addText='Cat. '+catName+ctrl_text,
                                       logy=False, legBkg=True,
                                       min_y=1, tag=tag+var, legLoc=[0.77, 0.55, 0.94, 0.75])
            cAux.SaveAs(outdir+'/fig/'+var+'_'+tag+'.png')
            cAux.SaveAs(webFolder+'/'+var+'_'+tag+'.png')
            outCanvas.append(cAux)

            if var.endswith('eta'):
                h = hDic[var]['data']
                nbins = h.GetNbinsX()
                xmin = h.GetBinCenter(1) - 0.5*h.GetBinWidth(1)
                xmax = h.GetBinCenter(nbins) + 0.5*h.GetBinWidth(nbins)
                if nbins%2 != 0 or np.abs(xmax+xmin) > 1e-3:
                    continue
                title = title.replace('#eta', '|#eta|')
                hDicAux = {}
                for k, h in hDic[var].iteritems():
                    hDicAux[k] = rt.TH1D(h.GetName()+'abs', h.GetTitle(), int(nbins/2), 0, xmax)
                    for ib in range(1, int(nbins/2)+1):
                        c = h.GetBinContent(int(nbins/2) + ib) + h.GetBinContent(int(nbins/2) + 1 - ib)
                        dc = np.hypot(h.GetBinError(int(nbins/2) + ib), h.GetBinError(int(nbins/2) + 1 - ib))
                        hDicAux[k].SetBinContent(ib, c)
                        hDicAux[k].SetBinError(ib, dc)
                hDicAux['data'].GetXaxis().SetTitle(title)
                hDicAux['data'].GetYaxis().SetTitle('Events')
                cAux = plot_SingleCategory(CMS_lumi, hDicAux,
                                           draw_pulls=True, pullsRatio=True, pulls_ylim='auto',
                                           scale_dic=scale_dic,
                                           addText='Cat. '+catName, logy=False, legBkg=True,
                                           min_y=1, tag=tag+var, legLoc=[0.77, 0.55, 0.94, 0.75])
                cAux.SaveAs(outdir+'/fig/'+var+'abs_'+tag+'.png')
                cAux.SaveAs(webFolder+'/'+var+'abs_'+tag+'.png')
                outCanvas.append(cAux)
        else:
            print '[WARNING] ' + var + ' missing from histos dict'




    # Draw control regions
    # for k in np.sort([k for k in hDic.keys() if 'AddTk' in k]):
    #     print 'Creating', k
    #     legLoc = [0.67, 0.3, 0.93, 0.72]
    #     if 'Vis' in k:
    #         legLoc = [0.18, 0.4, 0.4, 0.75]
    #     cAux = plot_SingleCategory(CMS_lumi, hDic[k], scale_dic=scale_dic,
    #                                xtitle=getControlXtitle(k),
    #                                addText='Cat. '+catName + ', ' + getControlSideText(k),
    #                                procOrder = ['tau', 'DstHc', 'dataSS_DstMu', 'mu', 'Dstst'],
    #                                tag=k, legLoc=legLoc,
    #                                draw_pulls=True,
    #                                pullsRatio=True if 'prefit' in tag else False,
    #                                pulls_ylim='auto'
    #                                )
    #     cAux.SaveAs(outdir+'/fig/'+k+'_'+tag+'.png')
    #     cAux.SaveAs(webFolder+'/'+k+'_'+tag+'.png')
    #     outCanvas.append(cAux)

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
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.DrawLatexNDC(0.17, 0.8, 'Cat. '+catName)
        cAux.SetLogz()
        cAux.SaveAs(outdir+'/fig/M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        cAux.SaveAs(webFolder+'/M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        outCanvas.append(cAux)

    #Re-rolling histos
    unrolledBins = None
    if args.category == 'comb':
        unrolledBins = pickle.load(open(outdir.replace('comb', catName.lower())+'/unrolledBinsMap.pkl', 'rb'))
    else:
        unrolledBins = pickle.load(open(outdir+'/unrolledBinsMap.pkl', 'rb'))
    hDic_reRollProj = {}
    for i_q2 in range(len(binning['q2'])-1):
        name2D = 'h2D_q2bin'+str(i_q2)
        nameU = 'Unrolled_q2bin'+str(i_q2)
        name2Dr = 'h2Dr_q2bin'+str(i_q2)
        hDic[name2Dr] = {}
        hDic_reRollProj['Est_mu_q2bin'+str(i_q2)] = {}
        hDic_reRollProj['M2_miss_q2bin'+str(i_q2)] = {}
        for n in hDic[nameU]:
            hDic[name2Dr][n] = hDic[name2D][n].Clone()
            hDic[name2Dr][n].Reset()
            for idx, (ix, iy) in enumerate(unrolledBins[i_q2]):
                hDic[name2Dr][n].SetBinContent(ix, iy, hDic[nameU][n].GetBinContent(idx+1))
                hDic[name2Dr][n].SetBinError(ix, iy, hDic[nameU][n].GetBinError(idx+1))
            # Re-rolled projections
            auxN = 'Est_mu_q2bin'+str(i_q2)
            hDic_reRollProj[auxN][n] = hDic[name2Dr][n].ProjectionY('hre_'+n+auxN, 1, -1, 'e')
            auxN = 'M2_miss_q2bin'+str(i_q2)
            hDic_reRollProj[auxN][n] = hDic[name2Dr][n].ProjectionX('hre_'+n+auxN, 1, -1, 'e')

    # Re-rolled 2D plot
    for i_q2 in range(len(binning['q2'])-1):
        name2D = 'h2Dr_q2bin'+str(i_q2)
        print 'Creating re-rolled histo', name2D
        hSum = hDic[name2D]['total']
        hSum.GetXaxis().SetTitle('m^{2}_{miss} [GeV^{2}]')
        hSum.GetYaxis().SetTitle('E_{#mu}* [GeV]')
        hSum.GetZaxis().SetTitle('Events')
        cAux = drawOnCMSCanvas(CMS_lumi, [hSum], ['colz'], tag=str(i_q2), mR=0.17)
        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.DrawLatexNDC(0.17, 0.8, 'Cat. '+catName)
        cAux.SetLogz()
        cAux.SaveAs(outdir+'/fig/reRolled_M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        cAux.SaveAs(webFolder+'/reRolled_M2Miss_vs_EstMu_q2bin{}_{}_TotMC.png'.format(i_q2, tag))
        outCanvas.append(cAux)
    print 'Creating re-rolled signal region grid'
    cAux = plot_gridVarQ2(CMS_lumi, binning, hDic_reRollProj, draw_pulls=True,
                          pullsRatio=False, pulls_ylim=[0.9, 1.1],
                          scale_dic=scale_dic,
                          categoryText=catName, cNameTag=tag+'_reRolled',
                          iq2_maskData=[] if args.unblinded else [2, 3])
    cAux.SaveAs(outdir+'/fig/reRolled_signalRegion_'+tag+'.png')
    cAux.SaveAs(webFolder+'/reRolled_signalRegion_'+tag+'.png')
    outCanvas.append(cAux)

    hMuPt_all = None
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
                                   addText='Cat. '+catName+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=False, legBkg=True,
                                   min_y=1,
                                   tag=tag+'mu_pt_q2bin'+str(i_q2),
                                   legLoc=[0.7, 0.5, 0.9, 0.75],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                   )
        cAux.SaveAs(outdir+'/fig/muPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/muPt_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)
        if i_q2 == 0:
            hMuPt_all = {}
            for n in hDic[name].keys():
                hMuPt_all[n] = hDic[name][n].Clone()
        else:
            for n in hMuPt_all.keys():
                hMuPt_all[n].Add(hDic[name][n])

    print 'Creating mu_pt_all'
    hMuPt_all['data'].GetYaxis().SetTitle('Normalized events')
    cAux = plot_SingleCategory(CMS_lumi, hMuPt_all, scale_dic=scale_dic,
                               draw_pulls=True, pullsRatio=True,
                               addText='Cat. '+catName,
                               logy=False, legBkg=True,
                               min_y=0,
                               max_y='data',
                               pulls_ylim=[0.85, 1.15],
                               density=True,
                               tag=tag+'mu_pt_all',
                               legLoc=[0.7, 0.3, 0.9, 0.55],
                               maskData = False
                               )
    cAux.SaveAs(outdir+'/fig/muPt_allNorm_'+tag+'.png')
    cAux.SaveAs(webFolder+'/muPt_allNorm_'+tag+'.png')
    outCanvas.append(cAux)

    axName = {'Dst': 'D*', 'K':'K', 'pi':'#pi', 'pis':'#pi_{soft}'}
    for n in ['Dst', 'K', 'pi', 'pis']:
        h_all = None
        for i_q2 in range(len(binning['q2'])-1):
            q2_l = binning['q2'][i_q2]
            q2_h = binning['q2'][i_q2 + 1]
            name = n+'_pt_q2bin'+str(i_q2)
            if not name in hDic.keys(): continue
            print 'Creating', name
            hDic[name]['data'].GetXaxis().SetTitle(axName[n] + ' p_{T} [GeV]')
            hDic[name]['data'].GetYaxis().SetTitle('Events')
            cAux = plot_SingleCategory(CMS_lumi, hDic[name], scale_dic=scale_dic,
                                       draw_pulls=True, pullsRatio=False,
                                       addText='Cat. '+catName+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                       logy=False, legBkg=True,
                                       min_y=1,
                                       tag=tag+n+'_pt_q2bin'+str(i_q2),
                                       legLoc=[0.7, 0.4, 0.9, 0.75],
                                       maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                       )
            cAux.SaveAs(outdir+'/fig/'+n+'Pt_q2bin'+str(i_q2)+'_'+tag+'.png')
            cAux.SaveAs(webFolder+'/'+n+'Pt_q2bin'+str(i_q2)+'_'+tag+'.png')
            outCanvas.append(cAux)
            if i_q2 == 0:
                h_all = {}
                for nn in hDic[name].keys():
                    h_all[nn] = hDic[name][nn].Clone()
            else:
                for nn in h_all.keys():
                    h_all[nn].Add(hDic[name][nn])
        print 'Creating '+n+'Pt_all'
        h_all['data'].GetYaxis().SetTitle('Normalized events')
        cAux = plot_SingleCategory(CMS_lumi, h_all, scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=True,
                                   addText='Cat. '+catName,
                                   logy=False, legBkg=True,
                                   min_y=0,
                                   max_y='data',
                                   pulls_ylim=[0.85, 1.15],
                                   density=True,
                                   tag=tag+n+'_pt_all',
                                   legLoc=[0.7, 0.3, 0.9, 0.55],
                                   maskData = False
                                   )
        cAux.SaveAs(outdir+'/fig/'+n+'Pt_allNorm_'+tag+'.png')
        cAux.SaveAs(webFolder+'/'+n+'Pt_allNorm_'+tag+'.png')
        outCanvas.append(cAux)


    h_all = None
    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        name = 'mass_D0pismu_q2bin'+str(i_q2)
        if not name in hDic.keys(): continue
        print 'Creating', name
        hDic[name]['data'].GetXaxis().SetTitle('mass(D*#mu) [GeV]')
        hDic[name]['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic[name], scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=False,
                                   addText='Cat. '+catName+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=False, legBkg=True,
                                   min_y=1,
                                   tag='mass_D0pismu_q2bin'+str(i_q2),
                                   legLoc=[0.16, 0.45, 0.33, 0.8],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                   )
        cAux.SaveAs(outdir+'/fig/mass_D0pismu_q2bin'+str(i_q2)+'_'+tag+'.png')
        cAux.SaveAs(webFolder+'/mass_D0pismu_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)
        if i_q2 == 0:
            h_all = {}
            for nn in hDic[name].keys():
                h_all[nn] = hDic[name][nn].Clone()
        else:
            for nn in h_all.keys():
                h_all[nn].Add(hDic[name][nn])
    print 'Creating mass_D0pismu_all'
    h_all['data'].GetYaxis().SetTitle('Normalized events')
    cAux = plot_SingleCategory(CMS_lumi, h_all, scale_dic=scale_dic,
                               draw_pulls=True, pullsRatio=True,
                               addText='Cat. '+catName,
                               logy=False, legBkg=True,
                               min_y=0,
                               max_y='data',
                               pulls_ylim=[0.85, 1.15],
                               density=True,
                               tag=tag+'mass_D0pismu_all',
                               legLoc=[0.2, 0.45, 0.4, 0.8],
                               maskData = False
                               )
    cAux.SaveAs(outdir+'/fig/mass_D0pismu_allNorm_'+tag+'.png')
    cAux.SaveAs(webFolder+'/mass_D0pismu_allNorm_'+tag+'.png')
    outCanvas.append(cAux)

    for i_q2 in range(len(binning['q2'])-1):
        q2_l = binning['q2'][i_q2]
        q2_h = binning['q2'][i_q2 + 1]
        name = 'U_miss_q2bin'+str(i_q2)
        if not name in hDic.keys(): continue
        print 'Creating', name
        hDic[name]['data'].GetXaxis().SetTitle('U_{miss} [GeV]')
        hDic[name]['data'].GetYaxis().SetTitle('Events')
        cAux = plot_SingleCategory(CMS_lumi, hDic[name], scale_dic=scale_dic,
                                   draw_pulls=True, pullsRatio=False,
                                   addText='Cat. '+catName+', {:.1f} <  q^{{2}}  < {:.1f} GeV^{{2}}'.format(q2_l, q2_h),
                                   logy=False, legBkg=True,
                                   min_y=1,
                                   tag='Umiss_q2bin'+str(i_q2),
                                   legLoc=[0.16, 0.45, 0.33, 0.8],
                                   maskData = (not args.unblinded) and (False if i_q2 < 2 else True)
                                   )
        cAux.SaveAs(webFolder+'/U_miss_q2bin'+str(i_q2)+'_'+tag+'.png')
        outCanvas.append(cAux)


    print 30*'-' + '\n\n'
    return outCanvas

def drawShapeVarPlots(card, tag=''):
    print '-----> Creating shape variations plots for', card

    plotsDir = webFolder + '/shapeUncertainties' + tag
    if os.path.isdir(plotsDir):
        os.system('rm -rf '+plotsDir)
    os.system('mkdir -p '+plotsDir)
    os.system('cp {d}/../index.php {d}/'.format(d=plotsDir))

    print '--> REGION:'
    for fn in glob(histo_file_dir + card + '_*.root'):
        region = os.path.basename(fn).replace(card+'_', '').replace('.root', '')
        if region.startswith('h2D'):
            continue
        # if not region.startswith('ctrl_p') and not region.stratswith('M2_miss'):
        if not region.startswith('ctrl_m'):
            continue
        print ' ', region

        auxOut = os.path.join(plotsDir, region)
        os.system('mkdir -p '+auxOut)
        os.system('cp {d}/../index.php {d}/'.format(d=auxOut))

        tfile = rt.TFile.Open(fn, 'READ')
        keys = []
        for k in tfile.GetListOfKeys():
            keys.append(k.GetName())

        processes = [k for k in keys if not '__' in k and not k == 'data_obs']
        for p in processes:
            hCentral = tfile.Get(p).Clone('hCentral')
            hCentral.SetTitle('Central')
            hCentral.GetYaxis().SetTitle('Events')
            hCentral.GetXaxis().SetTitle(region)
            hCentral.GetYaxis().SetRangeUser(0, 1.3*hCentral.GetMaximum())
            if p == 'tau':
                hCentral.Scale(SM_RDst)

            shapeVarNames = []
            for k in keys:
                if (k.startswith(p + '__') or p=='total') and k.endswith('Up'):
                    shapeVarNames.append(k.split('__')[1][:-2])
            shapeVarNames = list(np.sort(np.unique(shapeVarNames)))

            for sVar in shapeVarNames:
                if p == 'total':
                    hUp = hCentral.Clone('hUp')
                    hUp.Reset()
                    hDown = hCentral.Clone('hDown')
                    hDown.Reset()
                    for pp in processes:
                        if pp=='total': continue
                        scale = SM_RDst if pp == 'tau' else 1.
                        if pp+'__'+sVar+'Up' in keys:
                            hUp.Add(tfile.Get(pp+'__'+sVar+'Up'), scale)
                            hDown.Add(tfile.Get(pp+'__'+sVar+'Down'), scale)
                        else:
                            hUp.Add(tfile.Get(pp), scale)
                            hDown.Add(tfile.Get(pp), scale)
                else:
                    hUp = tfile.Get(p+'__'+sVar+'Up').Clone('hUp')
                    hDown = tfile.Get(p+'__'+sVar+'Down').Clone('hDown')
                    if p == 'tau':
                        hUp.Scale(SM_RDst)
                        hDown.Scale(SM_RDst)

                hUp.SetLineColor(rt.kRed-4)
                hDown.SetLineColor(rt.kAzure+1)
                if hCentral.Integral() > 0:
                    hUp.SetTitle('+1#sigma (n:{:.0f}%)'.format(100*hUp.Integral()/hCentral.Integral()))
                    hDown.SetTitle('-1#sigma (n:{:.0f}%)'.format(100*hDown.Integral()/hCentral.Integral()))
                else:
                    hUp.SetTitle('+1#sigma (n:{:.0f})'.format(hUp.Integral()))
                    hDown.SetTitle('-1#sigma (n:{:.0f})'.format(hDown.Integral()))

                for i in range(1, hCentral.GetNbinsX()+1):
                    if hCentral.GetBinContent(i) == 0:
                        hCentral.SetBinContent(i, 1e-6)
                        hCentral.SetBinError(i, 1e-6)
                    if hUp.GetBinContent(i) == 0:
                        hUp.SetBinContent(i, 1e-6)
                    if hDown.GetBinContent(i) == 0:
                        hDown.SetBinContent(i, 1e-6)
                hl = [hCentral, hUp, hDown]

                for ih, h in enumerate(hl):
                    h.Sumw2(0)
                    h.SetMarkerStyle(2)
                    h.SetMarkerSize(2)
                    h.SetMarkerColor(h.GetLineColor())
                c = make_ratio_plot(hl,
                                    draw_opt='',
                                    leg_pos=[0.75,0.73,0.95,0.92],
                                    marginTop=0.062,
                                    ratio_bounds='auto')

                # # Drowing MC uncertianty
                h_tot = tfile.Get(p).Clone('h_tot')
                c.pad1.cd()
                h_tot.SetFillColor(rt.kGray)
                h_tot.SetLineWidth(0)
                h_tot.SetMarkerStyle(0)
                h_tot.Draw('sameE2')
                for h in hl:
                    h.Draw('same')

                g_up = rt.TGraph()
                g_up.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), 1)
                g_down = rt.TGraph()
                g_down.SetPoint(0, h_tot.GetBinCenter(1)-0.5*h_tot.GetBinWidth(1), 1)
                for i in range(1, h_tot.GetNbinsX()+1):
                    x_low = h_tot.GetBinCenter(i) - 0.5*h_tot.GetBinWidth(i)
                    x_up = h_tot.GetBinCenter(i) + 0.5*h_tot.GetBinWidth(i)
                    c_MC = h_tot.GetBinContent(i)
                    e_MC = h_tot.GetBinError(i)
                    if c_MC > 0:
                        y_up = (c_MC+e_MC)/c_MC
                        y_down = (c_MC-e_MC)/c_MC
                    else:
                        y_up = 1.
                        y_down = 1.
                    g_up.SetPoint(2*i-1, x_low, y_up)
                    g_up.SetPoint(2*i, x_up, y_up)
                    g_down.SetPoint(2*i-1, x_low, y_down)
                    g_down.SetPoint(2*i, x_up, y_down)
                g_up.SetPoint(2*i+1, x_up, 1)
                g_down.SetPoint(2*i+1, x_up, 1)
                g_up.SetFillColorAlpha(rt.kGray, 0.5)
                g_up.SetFillStyle(1)
                g_down.SetFillColorAlpha(rt.kGray, 0.5)
                g_down.SetFillStyle(1)
                g_up.SetLineWidth(0)
                c.leg.AddEntry(g_up, 'MC stat.', 'f')
                c.statMCratio = [g_up, g_down]

                c.pad2.cd()
                # c.cd(2)
                g_up.Draw('F')
                g_down.Draw('F')
                c.ln.DrawLine(h.GetXaxis().GetXmin(), 1, h.GetXaxis().GetXmax(), 1)
                for h in c.hratio_list:
                    h.Draw('same')

                c.pad1.cd()
                txt = rt.TLatex()
                txt.SetTextSize(0.05)
                txt.SetTextAlign(13)
                txt.DrawLatexNDC(0.18, 0.91, 'Region: '+region)
                txt.DrawLatexNDC(0.18, 0.85, 'Process: '+p)
                txt.DrawLatexNDC(0.18, 0.79, 'Shape: '+sVar)

                CMS_lumi.extraText = "       Simulation Internal"
                CMS_lumi.CMS_lumi(c, -1, 0)
                c.SaveAs(auxOut + '/'+p+'__'+sVar+'.png')
        tfile.Close()
    return

def drawPrePostFitComparison(histoPre, histoPost, tag=''):
    print 20*'-', 'Pre/Post fit comaprison', tag, 20*'-'
    CMS_lumi.integrated_lumi = None
    outCanvas = []

    for p in processOrder + ['total']:
        print p
        auxOutDir = outdir+'/fig/prePostFitComparison/'+p
        if not os.path.isdir(auxOutDir):
            os.system('mkdir -p '+auxOutDir)
        auxWebDir = webFolder+'/prePostFitComparison/'+p
        if not os.path.isdir(auxWebDir):
            os.system('mkdir -p '+auxWebDir)
            os.system('cp {t}/../../index.php {t}/../ '.format(t=auxWebDir))
            os.system('cp {t}/../../index.php {t}/ '.format(t=auxWebDir))

        for c in histoPost.keys():
            if '2D' in c:
                continue
            if not p in histoPre[c].keys():
                continue
            if not p in histoPost[c].keys():
                continue
            hPre = histoPre[c][p].Clone()
            hPre.SetTitle('Prefit')
            hPre.SetLineColor(rt.kRed-4)
            normPre = hPre.Integral()
            hPre.Scale(1./hPre.Integral(), 'width')
            hPre.Sumw2(0)
            hPre.GetXaxis().SetTitle(c)
            hPre.GetYaxis().SetTitle('Normalized entries')
            hPost = histoPost[c][p].Clone()
            hPost.SetLineColor(rt.kAzure+1)
            normPost = hPost.Integral()
            if hPost.Integral() != 0:
                hPost.Scale(1./hPost.Integral(), 'width')
            hPost.SetTitle('Postfit ({:.1f}%)'.format(100*normPost/float(normPre)))
            hPost.Sumw2(0)

            for i in range(1, hPre.GetNbinsX()+1):
                if hPre.GetBinContent(i) == 0:
                    hPre.SetBinContent(i, 1e-9)
                    hPost.SetBinContent(i, 1e-9)

            can = make_ratio_plot([hPre, hPost],
                                 draw_opt='',
                                 leg_pos=[0.6,0.8,0.8,0.92],
                                 marginTop=0.062,
                                 label = c+p+tag,
                                 ratio_bounds='auto')
            CMS_lumi.CMS_lumi(can, -1, 33)
            can.SaveAs(auxOutDir+'/'+c+'_'+tag+'.png')
            can.SaveAs(auxWebDir+'/'+c+'_'+tag+'.png')

            outCanvas.append(can)

    return outCanvas

########################### -------- Create the card ------------------ #########################

def createSingleCard(histo, category, fitRegionsOnly=False):

    processes = processOrder
    nProc = len(processes)

    categories = []
    for c in np.sort(histo.keys()):
        if c.startswith('h2'): continue
        if fitRegionsOnly:
            if c.startswith('ctrl_'):
                if not c[5:] in args.controlRegions:
                    continue
            if args.signalRegProj1D:
                aux = c.startswith(args.signalRegProj1D)
            else:
                aux = c.startswith('Unrolled')
            aux = aux or c.startswith('ctrl_')
            if not aux:
                continue
            if (not args.unblinded) and (c.endswith('_q2bin2') or c.endswith('_q2bin3')):
                continue
            if args.noLowq2 and (c.endswith('_q2bin0') or c.endswith('_q2bin1')):
                continue
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
    #### pp -> bb cros-section * luminosity
    card += 'overallMcNorm'+category.trg+' rateParam * mu 1.\n'
    card += 'overallMcNorm'+category.trg+' rateParam * tau 1.\n'
    card += 'overallMcNorm'+category.trg+' rateParam * B[usd]_* 1.\n'

    # card += 'overallMcNormBd'+category.trg+' rateParam * mu 1.\n'
    # card += 'overallMcNormBd'+category.trg+' rateParam * tau 1.\n'
    # card += 'overallMcNormBd'+category.trg+' rateParam * Bd_* 1.\n'
    # card += 'overallMcNormBu'+category.trg+' rateParam * Bu_* 1.\n'
    # card += 'overallMcNormBs'+category.trg+' rateParam * Bs_* 1.\n'

    # Relax green unc
    # card += 'overallBToDstHc rateParam * B[dsu]_DstD[dsu] 1.\n'


    #### Combinatorial background norm
    val = ' 1.50'
    aux = ''
    for n in processes:
        if n == 'dataSS_DstMu':
            aux += val
        else:
            aux += ' -'
    # card += 'normDataSS'+category.trg+' lnN' + aux*nCat + '\n'

    #### Tracking efficiency uncertainty
    uncVal = 1.03
    card += 'trkEff lnN'
    for c in categories:
        val = ''
        if re.match('ctrl_[pm]__', c):
            val = ' {:.3}'.format(uncVal)
        elif re.match('ctrl_[pm]{2}_', c):
            val = ' {:.3}'.format(uncVal**2)
        else:
            val = ' -'
        card += val*len(processes)
    card += '\n'

    #### Branching ratio uncertainty
    brPklLoc = '/storage/af/group/rdst_analysis/BPhysics/data/forcedDecayChannelsFactors_v2.pickle'
    decayBR = pickle.load(open(brPklLoc, 'rb'))

    def brScaleSys(name, relevantProcesses=[], relUnc=0.01):
        val = ' {:.2f}'.format(1+relUnc)
        aux = ''
        for nn in processes:
            if nn in relevantProcesses:
                aux += val
            else:
                aux += ' -'
        return name + ' lnN' + aux*nCat + '\n'

    if args.freeMuBr:
        pass
        # card += 'mutauNorm rateParam * tau 1. 0.9,1.1\n'
        # card += 'mutauNorm rateParam * mu 1. 0.9,1.1\n'
    else:
        card += brScaleSys('brBd_DstMuNu', ['mu', 'tau'], relUnc=1.2/50.6)

    # card += brScaleSys('DstPiBr', ['Bu_MuDstPi', 'Bd_MuDstPi', 'Bu_TauDstPi', 'Bd_TauDstPi'], relUnc=2*0.4/6.0) # Moved to shape

    ## Temporary
    # card += brScaleSys('DstPiPiBr', ['Bu_MuDstPiPi', 'Bd_MuDstPiPi', 'Bu_TauDstPiPi', 'Bd_TauDstPiPi'], relUnc=0.3/0.96)
    card += brScaleSys('DstPiPiBr', ['Bu_MuDstPiPi', 'Bu_TauDstPiPi', 'Bd_TauDstPiPi'], relUnc=0.3/0.96)
    ##
    card += brScaleSys('DstKBr', ['Bs_MuDstK', 'Bs_TauDstK'], relUnc=1.5/5.9)

    card += brScaleSys('RDs_stst', ['Bu_TauDstPi', 'Bd_TauDstPi', 'Bd_TauDstPiPi', 'Bu_TauDstPiPi', 'Bs_TauDstK'], relUnc=0.3)

    card += brScaleSys('DuMuBr', ['Bd_DstDu', 'Bu_DstDu'], relUnc=2.5/60.8)
    card += brScaleSys('DdMuBr', ['Bd_DstDd', 'Bu_DstDd'], relUnc=2.7/158.8)
    card += brScaleSys('DsMuBr', ['Bd_DstDs', 'Bs_DstDs'], relUnc=2.1/75.4)

    card += brScaleSys('Bs_DstDsBr', ['Bs_DstDs'], relUnc=0.5) #They have not been observed so we variate them alltogether like this

    card += 60*'-'+'\n'

    ######################################################
    ########## Shape systematics uncertainties
    ######################################################
    mcProcStr = ''
    for p in processes:
        if 'data' in p: mcProcStr += ' -'
        else: mcProcStr += ' 1.'

    nameSF = 'trg{}SF'.format(category.trg)
    counter = 0
    for k in histo.values()[0].keys():
        if k.startswith(processOrder[0]+'__'+nameSF + '_pt') and k.endswith('Up'):
            n = k[k.find('__')+2:-2]
            card += n+' shape' + mcProcStr*nCat + '\n'
            counter += 1
    print 'Trigger SF unc', counter

    # card += 'muonIdSF shape' + ' 1.'*nProc*nCat + '\n'

    aux = ''
    for c in categories:
        if c.startswith('ctrl_'):
            aux += mcProcStr
        else: aux += ' -'*nProc
    if '1.' in aux:
        auxTag = '' if args.correlate_tkPVfrac else category.name
        card += 'tkPVfrac'+auxTag+' shape' + aux + '\n'

    # Soft track efficiency
    # card += 'softTrkEff_w shape' + mcProcStr*nCat + '\n'
    # card += 'softTrkEff_s shape' + mcProcStr*nCat + '\n'
    for k in histo.values()[0].keys():
        if k.startswith(processOrder[0]+'__softTrkEff') and k.endswith('Up'):
            n = k[k.find('__')+2:-2]
            card += n+' shape' + mcProcStr*nCat + '\n'

    # B eta uncertainty
    names = []
    for k in histo.values()[0].keys():
        if k.startswith(samples_Bd[0]+'__B_eta'+category.name) and k.endswith('Up'):
            names.append(k[len(samples_Bd[0]) + 2:-2])
    for ii, n in enumerate(sorted(names)):
        card += n + ' shape' + mcProcStr*nCat + '\n'
        if ii == 2: break

    # B pT uncertainty
    if not args.calBpT == 'none':
        # Bd pT spectrum
        aux = ''
        for p in processes:
            if p in samples_Bd:
                aux += ' 1.'
            else:
                aux += ' -'
        names = []
        for k in histo.values()[0].keys():
            if k.startswith(samples_Bd[0]+'__BdpT'+category.name) and k.endswith('Up'):
                names.append(k[len(samples_Bd[0]) + 2:-2])
        for n in sorted(names):
            card += n + ' shape' + aux*nCat + '\n'


    # Form Factors from Hammer
    if not args.freezeFF:
        for n_pFF in FreeParFF:
            aux = ''
            for p in processes:
                if p in ['tau', 'mu']:
                    aux += ' 1.'
                else:
                    aux += ' -'
            card += 'B2Dst'+schemeFF+'{} shape'.format(n_pFF) + aux*nCat + '\n'


    aux = ''
    for p in processes:
        if '_MuDstPi' in p:
            aux += ' 1.'
        else:
            aux += ' -'
    for nEig in [1,2,3,4]:
        card += 'FF_B2DststN_BLReig{} shape'.format(nEig) + aux*nCat + '\n'
        if nEig <= 3:
            card += 'FF_B2DststW_BLReig{} shape'.format(nEig) + aux*nCat + '\n'


    aux = ''
    for p in processes:
        if '_MuDstPiPi' in p:
            aux += ' 1.'
        else:
            aux += ' -'
    for parName in ['RhoSq', 'chi11', 'chi21', 'chi31', 'eta1']:
        card += 'FF_B2D2S_BLOP{} shape'.format(parName) + aux*nCat + '\n'


    # Dstst mix composition
    aux = ''
    for p in processes:
        if not re.search('DstPi\Z', p) is None:
            aux += ' 1.'
        else: aux += ' -'
    # card += 'fDststWide shape' + aux*nCat + '\n'
    auxList = ['brB_D2420MuNu','brB_D2430MuNu','brB_D2460MuNu','brB_DstPiMuNu']
    auxList += ['D2420_width', 'D2430_width', 'D2460_width']
    for shapeName in auxList:
        card += shapeName+' shape' + aux*nCat + '\n'


    aux = ''
    for p in processes:
        if not re.search('MuDstPiPi\Z', p) is None:
            aux += ' 1.'
        else: aux += ' -'
    for nnn in ['Dst2S_width',
                'brDstPiPi_NR',
                'brDstPiPi_D1', 'brDstPiPi_D1st', 'brDstPiPi_D2st', 'brDstPiPi_D2Sst', 'brDstPiPi_D2Sst_D1Pi', 'brDstPiPi_D2Sst_D1stPi', 'brDstPiPi_D2Sst_D2stPi',
                'brDstPiPi_D2S_D2stPi',
                'brDstPiPi_D1Pi', 'brDstPiPi_D1stPi', 'brDstPiPi_D2stPi']:
        card += nnn + ' shape'  + aux*nCat + '\n'
    # card += 'Dst2S_width shape' + aux*nCat + '\n'
    # card += 'brDstst_DststPi shape' + aux*nCat + '\n'
    # card += 'brD2420_DstPiPi shape' + aux*nCat + '\n'
    # card += 'brD2430_DstPiPi shape' + aux*nCat + '\n'
    # card += 'brD2460_DstPiPi shape' + aux*nCat + '\n'
    # card += 'brDst2S_DstPiPi shape' + aux*nCat + '\n'

    # Hc mix composition
    def brShapeSys(relevantSamples=[], shapeNames=[]):
        aux = ''
        for p in processes:
            if p in relevantSamples: aux += ' 1.'
            else: aux += ' -'
        out = ''

        for n in shapeNames:
            # type = 'shapeU' if 'Kst' in n else 'shape'
            type = 'shape'
            out += 'br'+n+' ' + type + aux*nCat + '\n'
        return out

    card += brShapeSys(['Bd_DstDu'], ['Bd_DstDuK', 'Bd_DstDustK', 'Bd_DstDust', 'Bd_DstDuKst'])
    card += brShapeSys(['Bd_DstDd'], ['Bd_DstDdK','Bd_DstDdstK','Bd_DstDdst','Bd_DstDd','Bd_DstDdKst'])
    card += brShapeSys(['Bd_DstDs'], ['Bd_DstDs', 'Bd_DstDsst', 'Bd_DstDsst0'])
    card += brShapeSys(['Bu_DstDu'], ['Bu_DstDuK', 'Bu_DstDustK', 'Bu_DstDuKst', 'Bu_DstDu'])
    card += brShapeSys(['Bu_DstDd'], ['Bu_DstDdK', 'Bu_DstDdstK', 'Bu_DstDdKst'])
    card += brShapeSys(['Bs_DstDs'], ['Bs_DstDsKst']) # Threated all together in scale sys above



    card += 60*'-'+'\n'

    ######################################################
    ########## MC statistical uncertainties
    ######################################################

    if not args.noMCstats:
        for r in args.controlRegions:
            card += 'ctrl_'+r+' autoMCStats 0 1 1\n'

        if args.useMVA:
            card += 'MVA autoMCStats 2 1 1\n'
        else:
            if args.signalRegProj1D:
                if not args.noLowq2:
                    card += args.signalRegProj1D+'_q2bin0 autoMCStats 0 1 1\n'
                    card += args.signalRegProj1D+'_q2bin1 autoMCStats 0 1 1\n'
                if args.unblinded:
                    card += args.signalRegProj1D+'_q2bin2 autoMCStats 0 1 1\n'
                    card += args.signalRegProj1D+'_q2bin3 autoMCStats 0 1 1\n'
            else:
                if not args.noLowq2:
                    card += 'Unrolled_q2bin0 autoMCStats 0 1 1\n'
                    card += 'Unrolled_q2bin1 autoMCStats 0 1 1\n'
                if args.unblinded:
                    card += 'Unrolled_q2bin2 autoMCStats 0 1 1\n'
                    card += 'Unrolled_q2bin3 autoMCStats 0 1 1\n'

        card += 60*'-'+'\n'

    ######################################################
    ########## Scorrelate systematics
    ######################################################

    signalChannel = args.signalRegProj1D if args.signalRegProj1D else 'Unrolled'

    # card += 'nuisance edit drop * * B2Dst'+schemeFF +'.* ifexists\n'

    # Relax prior increasing width by a factor 2
    # card += 'nuisance edit add * * B2Dst'+schemeFF +'.* shape 0.5 overwrite\n'

    if args.decorrelateFFpars:
        for n in FreeParFF:
            if n == 'R0':
                continue
            parName = 'B2Dst'+schemeFF+n
            card += 'nuisance edit rename * ' + signalChannel+'_q2bin[01] ' + parName + ' ' + parName+'_ctrlReg'+category.name+'\n'
            card += 'nuisance edit rename * ctrl_.* ' + parName + ' ' + parName+'_ctrlReg'+category.name+'\n'
            card += 'nuisance edit rename * ' + signalChannel+'_q2bin[23] ' + parName + ' ' + parName+'_sigReg'+category.name+'\n'
            card += 'nuisance edit drop * * ' + parName +'\n'

        card += 60*'-'+'\n'

    return card

def createCombinationCard(fitRegionsOnly=False):
    clFull = card_location.replace('.txt', '_fitRegionsOnly.txt') if fitRegionsOnly else card_location
    cmd = 'cd '+os.path.dirname(clFull)+'; '
    cl = os.path.basename(clFull)
    cmd += 'combineCards.py'
    for c in categoriesToCombine:
        singleCatCardLoc = clFull.replace('comb', c)
        nWait = 0
        while not os.path.isfile(singleCatCardLoc):
            if nWait > 10:
                print '[ERROR] Waited too long...goodbye.'
                raise
            print 'Waiting for {} card to be produced'.format(c)
            time.sleep(60)
            nWait += 1
        cmd += ' {}={}'.format(c, cl.replace('comb', c))
    cmd += ' > ' + cl
    runCommandSafe(cmd)

    # Editing the nuisace renaiming
    cardStream = open(clFull, 'r')
    lines = cardStream.readlines()
    cardStream.close()

    cardStream = open(clFull, 'w')
    # Copy the whole file w/o the nuisance edit lines
    for line in lines:
        if not line.startswith('nuisance edit'):
            cardStream.write(line)
    # Re-write them down
    if args.decorrelateFFpars:
        signalChannel = args.signalRegProj1D if args.signalRegProj1D else 'Unrolled'
        for n in FreeParFF:
            if n == 'R0':
                continue
            parName = 'B2Dst'+schemeFF+n
            for c in categoriesToCombine:
                chName = c + '_' + signalChannel+'_q2bin[01]'
                cardStream.write('nuisance edit rename * ' + chName + ' ' + parName + ' ' + parName+'_ctrlReg'+c.capitalize()+'\n')
                chName = c + '_ctrl_.*'
                cardStream.write('nuisance edit rename * ' + chName + ' ' + parName + ' ' + parName+'_ctrlReg'+c.capitalize()+'\n')
                chName = c + '_' + signalChannel+'_q2bin[23]'
                cardStream.write('nuisance edit rename * ' + chName + ' ' + parName + ' ' + parName+'_sigReg'+c.capitalize()+'\n')
            cardStream.write('nuisance edit drop * * ' + parName +'\n')

    cardStream.close()


########################### -------- Create the workspace ------------------ #########################

def createWorkspace(cardLoc):
    print '-----> Creating workspace'
    print cardLoc
    cmd = 'text2workspace.py ' + cardLoc
    cmd += ' -o ' + cardLoc.replace('.txt', '.root')
    cmd += ' --no-b-only --verbose 1 --channel-masks'
    # cmd += ' --no-wrappers'
    output = runCommandSafe(cmd)

    text_file = open(cardLoc.replace('.txt', '_text2workspace.out'), 'w')
    text_file.write(output)
    text_file.close()

    text_file = open(webFolder + '/' + os.path.basename(cardLoc).replace('.txt', '_text2workspace.out'), 'w')
    text_file.write(output)
    text_file.close()

########################### -------- Bias studies ------------------ #########################

def biasToysScan(card, out, seed=1, nToys=10, rVal=SM_RDst, maskStr=''):

    if not args.asimov:
        inputSpace = 'higgsCombineBestfit.MultiDimFit.mH120.root'
        if not os.path.isfile(os.path.join(out, inputSpace)):
            print '-------- Best fit snap'
            cmd = 'cd ' + out + '; '
            cmd += 'combine -M MultiDimFit'
            cmd += ' -d ' + card.replace('.txt', '.root')
            cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
            cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
            cmd += ' --setParameters r={:.2f}'.format(rVal)
            if maskStr:
                cmd += ','+maskStr
            cmd += ' --setParameterRanges r=0.1,0.5'
            cmd += ' -n Bestfit'
            cmd += ' --saveWorkspace --verbose 0'
            runCommandSafe(cmd)
        arr = rtnp.root2array(os.path.join(out, inputSpace), treename='limit')
        rVal = arr['r'][0]
        print 'Using r best fit value {:.4f}'.format(rVal)
    else:
        inputSpace = card.replace('.txt', '.root')



    print '-----> Generating toys (seed: {})'.format(seed)
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M GenerateOnly'
    cmd += ' -d ' + inputSpace
    cmd += ' --seed ' + str(seed)
    cmd += ' --noMCbonly 1'
    if args.asimov:
        cmd += ' --setParameters r={} --freezeParameters r'.format(rVal)
    else:
        cmd += ' --snapshotName MultiDimFit'
    cmd += ' --toysFrequentist -t {} --saveToys'.format(nToys)
    cmd += ' -n Toys -m {:.0f}'.format(1000*rVal)
    runCommandSafe(cmd)

    print '-----> Running the toys scans'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit'
    # cmd += ' --algo grid --points=100 -n Scan'
    cmd += ' --algo singles -n Singles'
    cmd += ' --robustFit 1 --cminDefaultMinimizerStrategy 1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --seed ' + str(seed)
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --toysFrequentist --toysFile higgsCombineToys.GenerateOnly.mH{:.0f}.{}.root -t {}'.format(1000*rVal, seed, nToys)
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' --setParameterRanges r=0.1,0.50'
    cmd += ' --trackParameters rgx{.*}'
    # cmd += ' --trackErrors rgx{.*}'
    cmd += ' -m {:.0f}'.format(1000*rVal)
    runCommandSafe(cmd)

def collectBiasToysResults(scansLoc, rVal=SM_RDst):
    print '-----> Collectiong bias toys scans'
    if not scansLoc[-1] == '/': scansLoc += '/'
    fnames = glob(scansLoc + 'higgsCombineScan.MultiDimFit.mH{:.0f}.*.root'.format(1000*rVal))
    res = None
    # seedsScans = []
    # for fname in fnames:
    #     idx = fname.find('.mH')
    #     seed = int(fname[idx+7:-5])
    #     seedsScans.append(str(seed))
    #     if seed > 99 and seed < 199:
    #         continue
    #     # print seed
    #     if res is None:
    #         res = getUncertaintyFromLimitTree(fname, verbose=False)
    #     else:
    #         res = np.concatenate((res, getUncertaintyFromLimitTree(fname, verbose=False)), axis=0)
    # if not res is None:
    #     r = res[:,0]
    #     rLoErr = res[:,1]
    #     rHiErr = res[:,2]
    # print 'Scan: ', ' '.join(seedsScans)

    fnames = glob(scansLoc + 'higgsCombineSingles.MultiDimFit.mH{:.0f}.*.root'.format(1000*rVal))
    seedsSingles = []
    trackedParam = None
    for fname in fnames:
        idx = fname.find('.mH')
        seed = int(fname[idx+7:-5])
        if seed > 300:
            continue
        seedsSingles.append(str(seed))

        auxRes, auxTracked = getResultsFromMultiDimFitSingles(fname, verbose=False, getTrackedParam=True)
        if res is None:
            res = auxRes
            trackedParam = auxTracked
        else:
            res = np.concatenate((res, auxRes), axis=0)
            for n in trackedParam.keys():
                trackedParam[n] += auxTracked[n]

    r = res[:,0]
    rLoErr = res[:,2]
    rHiErr = res[:,3]
    print 'Singles: ', ' '.join(seedsSingles)

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(8,6))
    plt.errorbar(np.arange(1, 1+r.shape[0]), r, yerr=np.column_stack((rLoErr, rHiErr)).T, fmt='o', color='#1f77b4', label='Toys fit results')
    m = np.mean(r)
    sm = np.std(r)/np.sqrt(r.shape[0])
    x = [0, r.shape[0]+1]
    plt.fill_between(x, 2*[m-sm], 2*[m+sm], color='#ff7f0e', alpha=0.4)
    plt.plot(x, 2*[m], color='#d62728', lw=1, label='Toys mean')
    plt.plot(x, [rVal, rVal], 'm--', lw=2, label='Injected value')
    ymin, ymax = plt.ylim(np.min(r - 2*rLoErr), np.max(r + 2*rHiErr))
    xmin, xmax = plt.xlim()
    plt.text(xmin + 0.2*(xmax-xmin), ymin + 0.07*(ymax-ymin), 'Estimated bias: $({:.2f} \pm {:.2f}) \cdot 10^{{-2}}$ '.format(100*(m-rVal), 100*sm))
    plt.legend(loc='upper right', numpoints=1)
    plt.xlabel('Toy number')
    plt.ylabel(r'$R(D^*)$')
    plt.savefig(outdir + '/fig/biasStudy_toysResults.png')

    webFolderBias = webFolder + '/biasStudy'
    if not os.path.isdir(webFolderBias):
        os.system('mkdir -p '+webFolderBias)
        os.system('cp {d}/../index.php {d}'.format(d=webFolderBias))
    plt.savefig(webFolderBias + '/toysResults.png')

    z = (r - rVal)/(0.5*(rLoErr + rHiErr))
    h = create_TH1D(z, name='hZtest', binning=[int(2*np.sqrt(r.shape[0])), -4, 4], axis_title=['#hat{R(D*)} - R(D*) / #sigma', 'Number of toys'])
    h.Sumw2()
    h.Fit('gaus', 'ILQ')
    rt.gStyle.SetStatY(0.95)
    c = drawOnCMSCanvas(CMS_lumi, [h])
    c.SaveAs(outdir + '/fig/biasStudy_zTest.png')
    c.SaveAs(webFolderBias + '/zTest.png')

    if not trackedParam is None:
        for n, x in trackedParam.iteritems():
            h = create_TH1D(np.array(x), name='h'+n,
                            binning=[int(2*np.sqrt(r.shape[0])), np.min(x), np.max(x)],
                            axis_title=[n, 'Number of toys'])
            h.Sumw2()
            h.Fit('gaus', 'ILQ')
            rt.gStyle.SetStatY(0.95)
            c = drawOnCMSCanvas(CMS_lumi, [h])
            c.SaveAs(webFolderBias + '/bestFitDistribution_'+n+'.png')


########################### -------- Likelihood scan ------------------ #########################

def dumpNuisFromScan(tag, out):
    print 'Dumping nuisances'
    name = out+'higgsCombine{}.MultiDimFit.mH120.root'.format(tag)
    df = pd.DataFrame(rtnp.root2array(name, treename='limit'))
    aux = []
    for v in df.columns:
        if not v.startswith('trackedParam_'):
            continue
        if v == 'trackedParam_CMS_th1x':
            continue
        aux.append([v[13:], df[v].iloc[0] ])

    t = PrettyTable()
    t.field_names = ['Parameter', 'Best fit']
    for var, val in sorted(aux, key=lambda x: np.abs(x[1]), reverse=True):
        t.add_row([var, '{:.2f}'.format(val)])

    with open(webFolder+'/scanNuisanceOut_'+tag+'.txt', 'w') as dumpfile:
        dumpfile.write('{}\n'.format(t))

    # Create pulls distribution
    hNuisances = create_TH1D(np.array([x[1] for x in aux]), 'hNuisDiff', binning=[41, -4.5, 4.5],
                             axis_title=['Post-pre fit nuisance difference [#sigma]', '# Nuisance'])
    hNuisances.Sumw2()
    binWidth = hNuisances.GetBinWidth(1)
    fGaus = rt.TF1('fFit', '{}*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))/({:.4f}*[1])'.format(len(aux)*binWidth, np.sqrt(2*np.pi)), -5, 5)
    fGaus.SetParameters(0,1)
    fGaus.SetParNames('#mu','#sigma')
    fGaus.SetLineColor(rt.kRed-4)
    fGaus.SetLineStyle(7)
    fGaus.SetLineWidth(3)
    hNuisances.Fit(fGaus, 'QWL')
    cAux = drawOnCMSCanvas(CMS_lumi, [hNuisances], tag='nuisDiff')
    textPave = hNuisances.FindObject('stats')
    textPave.SetY1NDC(0.7)
    textPave.SetY2NDC(0.95)
    cAux.SaveAs(webFolder+'/scanNuisanceOut_'+tag+'_distribution.png')


def runScan(tag, card, out, catName, rVal=SM_RDst, rLimits=[0.1, 0.7], nPoints=30, maskStr='', strategy=1, freezePars=[], draw=True, dumpNuis=False):
    if not out[-1] == '/': out += '/'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit'
    cmd += ' --algo grid --points='+str(nPoints)
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy='+str(strategy)
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --X-rtd MINIMIZER_analytic'
    cmd += ' -d ' + card.replace('.txt', '.root')
    if len(rLimits) == 1:
        cmd += ' --centeredRange {:.2f}'.format(rLimits[0])
    else:
        cmd += ' --rMin={:.2f} --rMax={:.2f}'.format(*rLimits)
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    if maskStr:
        cmd += ','+maskStr
    if len(freezePars):
        cmd += ' --freezeParameters ' + ','.join(freezePars)
    cmd += ' --trackParameters rgx{.*}'
    cmd += ' -n ' + tag
    cmd += ' --verbose -1'
    output = runCommandSafe(cmd)
    if args.verbose:
        print output

    if dumpNuis:
        dumpNuisFromScan(tag, out)

    if draw:
        json.dump({'r': 'R(D*)'}, open(out+'renameDicLikelihoodScan.json', 'w'))

        cmd = 'cd ' + out + '; '
        cmd += '../../plot1DScan.py higgsCombine{t}.MultiDimFit.mH120.root -o scan{t}'.format(t=tag)
        cmd += ' --main-label "{} {}'.format('Obs.' if not args.asimov else 'Asimov', catName)
        if not args.unblinded: cmd += ' (blinded)'
        cmd += '"'
        cmd += ' --translate ' + out+'renameDicLikelihoodScan.json'
        runCommandSafe(cmd)
        cmd = 'cp {}scan{}.png {}/'.format(out, tag, webFolder)
        runCommandSafe(cmd)

    print 'Extracting new POI boundaries'
    res = getUncertaintyFromLimitTree(out+'higgsCombine{}.MultiDimFit.mH120.root'.format(tag))
    rLimitsOut = [res[0] - 3*res[1], res[0] + 3*res[2]]
    rValOur = res[0]
    with open(webFolder+'/scan_results.txt', 'a') as dumpFile:
        strRes = tag
        strRes += ' '*(35-len(strRes))
        strRes += '{:.3f} +{:.3f} / -{:.3f}'.format(res[0], res[2], res[1])
        dumpFile.write(strRes + '\n')
    if args.showPlots:
        display(Image(filename=out+'/scan'+tag+'.png'))
    return rValOur, rLimitsOut


########################### -------- Categories compatibility -------- #########################

def categoriesCompatibility(card, out, rVal=SM_RDst, rLimits=[0.1, 0.7]):
    fLog = open(webFolder + '/categoriesCompatibility.txt', 'w')
    print '----- Running nominal fit'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo singles'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --robustFit 1 --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    cmd += ' --setParameterRanges r={:.4f},{:.4f}'.format(*rLimits)
    cmd += ' -n _catCompNominal'
    runCommandSafe(cmd)
    arr = rtnp.root2array(out + '/higgsCombine_catCompNominal.MultiDimFit.mH120.root', treename='limit')
    rBestFit, rBestFitDown, rBestFitUp = arr['r']
    s = 'Using r best fit value {:.4f} +/- {:.4f}'.format(rBestFit, 0.5*(rBestFitUp-rBestFitDown))
    print s
    fLog.write(s + '\n')

    print '----- Creating workspace with independent signal strength'
    wsLoc = card.replace('.txt', '_rCat.root')
    if not card.endswith('.txt'):
        print 'categoriesCompatibility needs txt card input.'
        print card
        raise
    cmd = 'text2workspace.py ' + card
    cmd += ' -o ' + wsLoc
    cmd += ' -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel'
    for c in categoriesToCombine:
        cmd += ' --PO map=\'.*'+c+'.*/tau:r' + c.capitalize() + '[{},{},{}]\''.format(rVal, rLimits[0], rLimits[1])
    cmd += ' --no-b-only --verbose 1'
    output = runCommandSafe(cmd)
    text_file = open(card.replace('.txt', '_rCat_text2workspace.out'), "w")
    text_file.write(output)
    text_file.close()

    print '----- Running fit with r uncorrelated in each category'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo singles'
    cmd += ' -d ' + wsLoc
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -n _catComp_rCatIndep'
    runCommandSafe(cmd)

    arr = rtnp.root2array(out + '/higgsCombine_catComp_rCatIndep.MultiDimFit.mH120.root', treename='limit')
    rFit = []
    for c in categoriesToCombine:
        x = arr['r'+c.capitalize()]
        rFit.append( [x[0], 0.5*(np.max(x) - np.min(x))] )
        fLog.write(c + ': {:.4f} +/- {:.4f}'.format(*rFit[-1]) + '\n')
    rFit = np.array(rFit)
    rMean = np.sum(rFit[:,0]/np.square(rFit[:,1])) / np.sum(1./np.square(rFit[:,1]))
    print 'Average observed r: {:.4f}'.format(rMean)
    fLog.write('Average observed r: {:.4f}'.format(rMean) + '\n')

    chi2 = np.sum(np.square((rFit[:,0] - rMean)/rFit[:,1]))
    dof = rFit.shape[0] - 1
    pval = scipy_chi2.sf(chi2, dof)
    print 'Chi2 = {:.2f} ({:.1f}%)'.format(chi2, 100*pval)
    fLog.write('Chi2 = {:.2f} ({:.1f}%)'.format(chi2, 100*pval) + '\n')

    print '----- Running fit with r fixed to bestfit in each category'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo fixed'
    cmd += ' --fixedPointPOIs ' + ','.join(['r{}={:.4f}'.format(c.capitalize(), rBestFit) for c in categoriesToCombine])
    cmd += ' -d ' + wsLoc
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -n _catComp_rCatFixed'
    runCommandSafe(cmd)
    arr = rtnp.root2array(out + '/higgsCombine_catComp_rCatFixed.MultiDimFit.mH120.root', treename='limit')
    chi2 = 2*arr['deltaNLL'][1]
    pval = scipy_chi2.sf(chi2, dof)
    print 'Wilks Chi2 = {:.2f} ({:.1f}%)'.format(chi2, 100*pval)
    fLog.write('Wilks Chi2 = {:.2f} ({:.1f}%)'.format(chi2, 100*pval) + '\n')
    fLog.close()

    return

########################### -------- Fit Diagnostic ------------------ #########################
def defineChannelMasking(histo):
    visibleChannels = ['ctrl_'+r for r in args.controlRegions]
    fitVar = args.signalRegProj1D if args.signalRegProj1D else 'Unrolled'
    if not args.noLowq2:
        visibleChannels += [fitVar+'_q2bin0', fitVar+'_q2bin1']
    if args.unblinded:
        visibleChannels += [fitVar+'_q2bin2', fitVar+'_q2bin3']

    if args.category == 'comb':
        visibleChannels = ['rgx{mask_[hml].*_'+r+'}' for r in visibleChannels]
    else:
        visibleChannels = ['mask_'+r for r in visibleChannels]

    channelMaskingStr = 'rgx{mask_.*}=1,'+','.join([c+'=0' for c in visibleChannels])
    return channelMaskingStr

def runFitDiagnostic(tag, card, out, forceRDst=False, maskStr='', rVal=SM_RDst, rLimits=[0.1, 0.7], seed=6741, strategy=1):
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M FitDiagnostics'
    cmd += ' --robustFit 1 --robustHesse 1 --cminDefaultMinimizerStrategy '+str(strategy)+' --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --seed ' + str(seed)
    cmd += ' -d ' + card.replace('.txt', '.root')
    if forceRDst:
        cmd += ' --skipSBFit'
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
    cmd += ' --saveShapes --saveWithUncertainties --saveNormalizations  --saveWorkspace'
    cmd += ' --trackParameters rgx{.*}'
    cmd += ' --plots'
    cmd += ' --verbose -1'
    output = runCommandSafe(cmd)
    if rt.gROOT.IsBatch():
        print 50*'#'
        print 20*'#' + ' Fit Diag ' + 20*'#'
        print 50*'#'
        print output
    for line in output.split('\n'):
            if 'ERROR' in line: print line.replace('ERROR', '\033[1m\x1b[31mERROR\x1b[0m')
            if 'Error' in line: print line.replace('Error', '\033[1m\x1b[31mError\x1b[0m')
            if forceRDst:
                if 'customStartingPoint' in line: print line

    if not out[-1] == '/': out += '/'
    arr = rtnp.root2array(out + 'higgsCombine{}.FitDiagnostics.mH120.{}.root'.format(runName, seed), treename='limit')
    if forceRDst:
        if len(arr['limit']) > 1:
            print '[ERROR] Multiple values for R(D*):', arr['limit']
        else:
            print 'R(D*) value fixed to', arr['limit'][0]
    else:
        try:
            c, d, u, _ = arr['limit']
            print 'R(D*) = {:.3f} +{:.3f}/-{:.3f} [{:.1f} %]'.format(c, u-c, c-d, 100*(u-d)*0.5/c)
        except:
            print '[ERROR] Limit output format not recognized'
            print arr['limit']

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
    if args.category == 'comb':
        for catregName in [k.GetTitle() for k in fd.GetListOfKeys()]:
            c = catregName.split('_')[0]
            if not c in histo_postfit.keys():
                histo_postfit[c] = {}
            regName = catregName[len(c)+1:]
            histo_postfit[c][regName] = {}

            for n, h in histo_prefit[c][regName].iteritems():
                if '__' in n:
                    continue
                if n == 'data':
                    histo_postfit[c][regName]['data'] = h.Clone(h.GetName() + '_data')
                else:
                    h_post = h.Clone(h.GetName() + '_postfit')
                    h_fit = fd.Get(catregName+'/'+n)
                    if not h_fit:
                        print n+' missing from '+c+' '+regName
                        continue
                    for i in range(1, h_post.GetNbinsX()+1):
                        h_post.SetBinContent(i, h_fit.GetBinContent(i))
                        h_post.SetBinError(i, h_fit.GetBinError(i))

                    histo_postfit[c][regName][n] = h_post

            for k in histo_prefit[c].keys():
                if not k.startswith('h2D_q2bin'):
                    continue
                if k in histo_postfit[c].keys():
                    break
                histo_postfit[c][k] = {}
                for n in histo_prefit[c][k].keys():
                    histo_postfit[c][k][n] = histo_prefit[c][k][n].Clone()
                    histo_postfit[c][k][n].Reset()
    else:
        for regName in [k.GetTitle() for k in fd.GetListOfKeys()]:
            histo_postfit[regName] = {}

            for n, h in histo_prefit[regName].iteritems():
                if '__' in n:
                    continue
                if n == 'data':
                    histo_postfit[regName]['data'] = h.Clone(h.GetName() + '_data')
                else:
                    h_post = h.Clone(h.GetName() + '_postfit')
                    h_fit = fd.Get(regName+'/'+n)
                    if not h_fit:
                        print n+' missing from '+regName
                        continue
                    for i in range(1, h_post.GetNbinsX()+1):
                        h_post.SetBinContent(i, h_fit.GetBinContent(i))
                        h_post.SetBinError(i, h_fit.GetBinError(i))

                    histo_postfit[regName][n] = h_post

        for k in histo_prefit.keys():
            if not k.startswith('h2D_q2bin'): continue
            histo_postfit[k] = {}
            for n in histo_prefit[k].keys():
                histo_postfit[k][n] = histo_prefit[k][n].Clone()
                histo_postfit[k][n].Reset()

    h2 = fFitDiagnostics.Get('covariance_fit_' + ('b' if forceRDst else 's'))
    n = None
    nRateParamX = None
    for il, labObj in enumerate(h2.GetXaxis().GetLabels()):
        lab = labObj.GetName()
        if lab.startswith('prop_bin') and n is None:
            n = il
        elif not lab.startswith('prop_bin') and not (n is None) and nRateParamX is None:
            nRateParamX = il + 1
    if n is None:
        n = il

    nR = None
    nRateParamY = None
    for il, labObj in enumerate(reversed(h2.GetYaxis().GetLabels())):
        lab = labObj.GetName()
        if lab == 'r':
            nR = il+1
        elif lab.startswith('prop_bin'):
            nRateParamY = il
            break

    h2.Scale(100.)
    rt.gStyle.SetPaintTextFormat('.0f')
    N = h2.GetNbinsX()
    # n=80
    h2.LabelsOption("v")
    gSF = 70/float(n)

    if not nR is None:
        h2.SetMarkerSize(2.0*gSF)
        h2.GetXaxis().SetLabelSize(0.07*gSF)
        h2.GetYaxis().SetLabelSize(0.07*gSF)
        h2.GetXaxis().SetRange(1, n)
        h2.GetYaxis().SetRangeUser(nR-1, nR)
        h2.GetZaxis().SetRangeUser(-100, 100)
        h2.GetZaxis().SetNdivisions(-304)
        CC1 = drawOnCMSCanvas(CMS_lumi, [h2, h2], ['colz', 'text same'], size=(1200, 300), tag='tl1', mL=0.03, mR=0.08, mB=0.65, mT=0.1)
        CC1.SaveAs(out+'fig/correlationR'+ ('_RDstFixed' if forceRDst else '')+'.png')
        CC1.SaveAs(webFolder+'/correlationR'+ ('_RDstFixed' if forceRDst else '')+'.png')
        CC1.SaveAs(webFolder+'/correlationR'+ ('_RDstFixed' if forceRDst else '')+'.pdf')

    h2.SetMarkerSize(.5*gSF)
    h2.GetXaxis().SetLabelSize(0.02*gSF)
    h2.GetYaxis().SetLabelSize(0.02*gSF)
    h2.GetXaxis().SetRange(1, n)
    h2.GetYaxis().SetRangeUser(N-n, N)
    h2.GetZaxis().SetRangeUser(-100, 100)
    h2.GetZaxis().SetNdivisions(510)
    CC = drawOnCMSCanvas(CMS_lumi, [h2, h2], ['colz', 'text same'], size=(900, 700), tag='tl', mL=0.12, mR=0.135, mB=0.16)
    CC.SaveAs(out+'fig/covariance_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')
    CC.SaveAs(webFolder+'/covariance_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')
    CC.SaveAs(webFolder+'/covariance_zoom'+ ('_RDstFixed' if forceRDst else '')+'.pdf')

    if nRateParamY>1:
        gSF = 5/float(nRateParamY)
        h2.SetMarkerSize(1.2*gSF)
        h2.GetXaxis().SetLabelSize(0.035*gSF)
        h2.GetYaxis().SetLabelSize(0.035*gSF)
        h2.GetXaxis().SetRange(nRateParamX, N)
        h2.GetYaxis().SetRangeUser(0, nRateParamY)
        h2.GetZaxis().SetRangeUser(-100, 100)
        h2.GetZaxis().SetNdivisions(510)
        CC = drawOnCMSCanvas(CMS_lumi, [h2, h2], ['colz', 'text same'], size=(900, 700), tag='tl', mL=0.2, mR=0.135, mB=0.25)
        CC.SaveAs(out+'fig/covariance_rateParam_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')
        CC.SaveAs(webFolder+'/covariance_rateParam_zoom'+ ('_RDstFixed' if forceRDst else '')+'.png')
        CC.SaveAs(webFolder+'/covariance_rateParam_zoom'+ ('_RDstFixed' if forceRDst else '')+'.pdf')

    return histo_postfit, CC, fFitDiagnostics

def extactCovarianceBlock(tag, out, parameters=['r', 'B2DstCLNeig1'], forceRDst=False):
    print 'Foction work in progress! To be compleated!!!!'
    exit()
    return

def nuisancesDiff(tag, out, forceRDst):
    runName = tag + ('_RDstFixed' if forceRDst else '')
    cmd = 'python diffNuisances.py ' + out + '/fitDiagnostics{}.root'.format(runName)
    if forceRDst:
        cmd += ' --skipFitSB'
    else:
        cmd += ' --skipFitB'
    cmd += ' --all'
    cmd += ' --abs'
    cmd += ' -g {}/nuisance_difference'.format(out) + runName + '.root'
    output = runCommandSafe(cmd)
    print 'Done'
    nName, nValPost, nSigma = dumpDiffNuisances(output, out, tag='RDstFixed' if forceRDst else '',
                      useBonlyResults=forceRDst, parsToPrint=100)
    cmd = 'cp {}/nuisance_difference{}*txt {}/'.format(out, '_RDstFixed' if forceRDst else '', webFolder)
    runCommandSafe(cmd)

    print 'Crating nuisances difference distribution'
    hNuisances = create_TH1D(nSigma, 'hNuisDiff', binning=[41, -4.5, 4.5],
                             axis_title=['Post-pre fit nuisance difference [#sigma]', '# Nuisance'])
    hNuisances.Sumw2()
    binWidth = hNuisances.GetBinWidth(1)
    fGaus = rt.TF1('fFit', '{}*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))/({:.4f}*[1])'.format(nSigma.shape[0]*binWidth, np.sqrt(2*np.pi)), -5, 5)
    fGaus.SetParameters(0,1)
    fGaus.SetParNames('#mu','#sigma')
    fGaus.SetLineColor(rt.kRed-4)
    fGaus.SetLineStyle(7)
    fGaus.SetLineWidth(3)
    hNuisances.Fit(fGaus, 'QWL')
    cAux = drawOnCMSCanvas(CMS_lumi, [hNuisances], tag='nuisDiff')
    textPave = hNuisances.FindObject('stats')
    textPave.SetY1NDC(0.7)
    textPave.SetY2NDC(0.95)
    cAux.SaveAs(outdir+'/fig/nuisanceDifferenceDistribution.png')
    cAux.SaveAs(webFolder+'/nuisanceDifferenceDistribution.png')

    return


########################### -------- Uncertainty breakdown ------------------ #########################

def runUncertaintyBreakDownScan(card, out, catName, rVal=SM_RDst, rLimits=[0.1, 0.7], maskStr=''):
    if not out[-1] == '/': out += '/'
    print '--------> Running uncertainty breakdown <--------------'
    print '--------> Nominal scan'
    rValOut, rLimitsOut = runScan('Nominal', card, out, catName, rVal, rLimits, nPoints=150, maskStr=maskStr, strategy=args.scanStrategy, draw=False)
    sig = (rLimitsOut[1] - rValOut)/3
    rLimitsTight = [rValOut - 2*sig, rValOut + 2*sig]

    print '--------> Best fit snap'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --setParameters r={:.2f}'.format(rValOut)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' --setParameterRanges r={:.3f},{:.3f}'.format(*rLimitsTight)
    cmd += ' -n Bestfit'
    cmd += ' --saveWorkspace --verbose -1'
    runCommandSafe(cmd)

    print '--------> Statistical uncertanty only'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=100'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -d higgsCombineBestfit.MultiDimFit.mH120.root'
    cmd += ' --snapshotName MultiDimFit'
    cmd += ' --rMin={:.3f} --rMax={:.3f}'.format(*rLimitsTight)
    cmd += ' -n StatOnly'
    cmd += ' --freezeParameters allConstrainedNuisances'
    if maskStr: cmd += ' --setParameters ' + maskStr
    cmd += ' --fastScan' # To be added if there are no free parameters otherwise
    cmd += ' --verbose -1'
    runCommandSafe(cmd)
    getUncertaintyFromLimitTree(out + 'higgsCombineStatOnly.MultiDimFit.mH120.root')

    print '--------> MC stats and Statistical uncertanty only'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=100'
    cmd += ' --cminDefaultMinimizerStrategy=1 --robustFit 1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -d higgsCombineBestfit.MultiDimFit.mH120.root'
    cmd += ' --rMin={:.4f} --rMax={:.4f}'.format(*rLimitsTight)
    cmd += ' -n MCstat'
    cmd += ' --snapshotName MultiDimFit'
    if maskStr: cmd += ' --setParameters ' + maskStr
    cmd += ' --freezeNuisanceGroups=autoMCStats'
    cmd += ' --verbose -1'
    runCommandSafe(cmd)
    getUncertaintyFromLimitTree(out + 'higgsCombineMCstat.MultiDimFit.mH120.root')


    json.dump({'r': 'R(D*)'}, open(out+'renameDicLikelihoodScan.json', 'w'))

    cmd = 'cd ' + out + '; '
    cmd += 'python ../../plot1DScan.py higgsCombineNominal.MultiDimFit.mH120.root'
    cmd += ' --main-label "{} {}'.format('Obs.' if not args.asimov else 'Asimov', catName)
    if not args.unblinded: cmd += ' (blinded)'
    cmd += '"'
    cmd += ' --others'
    cmd += ' "higgsCombineMCstat.MultiDimFit.mH120.root:Stat. + Syst.:4"'
    cmd += ' "higgsCombineStatOnly.MultiDimFit.mH120.root:Stat. only:2"'
    cmd += ' --breakdown "MC stat.,syst.,stat."'
    cmd += ' --translate ' + out+'renameDicLikelihoodScan.json'
    cmd += ' -o scanBreakdown'
    runCommandSafe(cmd)
    cmd = 'cp {}scanBreakdown.png {}/'.format(out, webFolder)
    cmd += '; cp {}scanBreakdown.pdf {}/'.format(out, webFolder)
    status, output = commands.getstatusoutput(cmd)
    runCommandSafe(cmd)
    if args.showPlots:
        display(Image(filename=out+'scanBreakdown.png'))



def runUncertaintyBreakDownTable(card, out, catName, rVal=SM_RDst, rLimits=[0.1, 0.7], maskStr=''):
    print '--------> Running uncertainty breakdown <--------------'
    if not out[-1] == '/':
        out += '/'
    out += 'uncertaintyBreakDownTable/'
    if not os.path.isdir(out):
        os.makedirs(out)


    uncRemaining = []
    uncAss = []
    uncNames = []

    fLog = open(webFolder + '/uncertaintyBreakDownTable_log.txt', 'w')
    print '----- Running nominal fit'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=200'
    cmd += ' -d ' + card.replace('.txt', '.root')
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --setParameters r={:.2f}'.format(rVal)
    if maskStr:
        cmd += ','+maskStr
    cmd += ' --setParameterRanges r={:.4f},{:.4f}'.format(*rLimits)
    cmd += ' -n _total'
    runCommandSafe(cmd)
    res = getUncertaintyFromLimitTree(out + 'higgsCombine_total.MultiDimFit.mH120.root', verbose=False, drawPlot=False)
    r, rDown, rUp = res[0], res[-3], res[-2]
    uncRemaining.append(0.5*(rUp-rDown))
    uncNames.append('total')
    s = 'Nominal fit: {:.4f} +/- {:.4f}'.format(r, uncRemaining[-1])
    print s
    fLog.write(s + '\n')
    fLog.close()

    rDown = r - 2.1*(r - rDown)
    rUp = r + 2.1*(rUp - r)
    idx = cmd.find('--setParameterRanges') + len('--setParameterRanges ')
    oldRange = cmd[idx: idx + 15]
    newRange = 'r={:.4f},{:.4f}'.format(rDown, rUp)
    cmd = cmd.replace(oldRange, newRange)

    print '-------- Best fit snap'
    cmdBestFit = cmd.replace(' --algo grid --points=200', '')
    cmdBestFit = cmdBestFit.replace(' -n _total', ' -n Bestfit')
    cmdBestFit += ' --saveWorkspace --verbose 0'
    runCommandSafe(cmdBestFit)

    print '-------- Breaking the uncertainty'
    def extractUncertainty(tag, type='scan'):
        if type=='fit':
            arr = rtnp.root2array(out + 'higgsCombine_'+tag+'.MultiDimFit.mH120.root', treename='limit')
            r, rDown, rUp = arr['r']
        elif type=='scan':
            res = getUncertaintyFromLimitTree(out + 'higgsCombine_'+tag+'.MultiDimFit.mH120.root', verbose=False, drawPlot=False)
            r, rDown, rUp = res[0], res[-3], res[-2]
        else:
            print 'Type not recognised'
            raise
        uncRemaining.append(0.5*(rUp-rDown))
        uncNames.append(tag)
        uncAss.append(np.sqrt(uncRemaining[-2]**2 - uncRemaining[-1]**2))
        s = tag + ': {:.4f} +/- {:.4f} ({:.4f})'.format(r, uncRemaining[-1], uncAss[-1])
        print s
        fLog = open(webFolder + '/uncertaintyBreakDownTable_log.txt', 'a')
        fLog.write(s + '\n')
        fLog.close()
        return 0.5*(rUp-rDown), [r - 2.1*(r - rDown), r + 2.1*(rUp - r)]

    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo grid --points=100'
    cmd += ' -d higgsCombineBestfit.MultiDimFit.mH120.root'
    cmd += ' --robustFit 1'
    cmd += ' --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    if maskStr:
        cmd += ' --setParameters '+maskStr
    cmd += ' --setParameterRanges ' + newRange
    cmd += ' --snapshotName MultiDimFit'
    cmd += ' --freezeNuisanceGroups=autoMCStats'

    # MC statistics uncertanty
    print '----> Freezing MC stat'
    cmdAux = cmd + ' -n _MCstat'
    runCommandSafe(cmdAux)
    dr, rLims = extractUncertainty('MCstat')
    idx = cmd.find('--setParameterRanges') + len('--setParameterRanges ')
    cmd = cmd.replace(cmd[idx: idx + 15], 'r={:.4f},{:.4f}'.format(*rLims))
    print ' '

    groupsDefFile = os.path.join(basedir, 'Combine/uncertaintyBreakdownTableGroups.yml')
    groups = yaml.load(open(groupsDefFile, 'r'))
    frozenNuisance = []
    firstToRun = ''
    for ig, group in enumerate(groups):
        if args.freezeFF and group['tag'] == 'formFactors':
            continue
        print '----> Freezing ' + group['tag']
        frozenNuisance += group['nuisance']
        cmdAux = cmd + ' -n _' + group['tag']
        cmdAux += ' --freezeParameters ' + ','.join(frozenNuisance)
        if not firstToRun or firstToRun == group['tag']:
            runCommandSafe(cmdAux)
            firstToRun = ''
        dr, rLims = extractUncertainty(group['tag'])
        idx = cmd.find('--setParameterRanges') + len('--setParameterRanges ')
        cmd = cmd.replace(cmd[idx: idx + 15], 'r={:.4f},{:.4f}'.format(*rLims))

        print ' '

    # Remove all of them (statistical)
    print '----> Freezing all nuisances'
    cmd = cmd.replace(' --freezeNuisanceGroups=autoMCStats', ' --freezeParameters allConstrainedNuisances')
    # cmd = cmd.replace(' --freezeNuisanceGroups=autoMCStats', ' --freezeParameters rgx{.*}')
    cmd += ' -n _remainingNuis'
    cmd += ' --fastScan'
    runCommandSafe(cmd)
    uncStat, _ = extractUncertainty('remainingNuis', type='scan')
    uncNames.append('stat')
    uncAss.append(uncStat)
    s = 'Statistics: XX +/- {:.4f} ({:.4f})'.format(uncStat, uncStat)
    print s
    fLog = open(webFolder + '/uncertaintyBreakDownTable_log.txt', 'a')
    fLog.write(s + '\n')

    dicDump = {'uncStat': uncStat, 'uncNames': uncNames, 'uncAss': uncAss}
    pickle.dump(dicDump, open(out+'/results.pkl', 'wb'))

    totExt = np.sqrt(np.sum(np.square(uncAss)))
    s = '\n\nUncertianty sum check: {:.4f}, Delta = {:1.2e}'.format(totExt, totExt - uncRemaining[0])
    print s
    fLog.write(s + '\n')
    fLog.close()

    print '\n\n-----------> Creating latex table'
    supplementDic = {'MCstat'       : 'Finite MC sample size',
                     'stat'         : 'Statistical',
                     'remainingNuis': 'Others'
                    }
    fTable = open(webFolder + '/uncertaintyBreakDownTable_latex.txt', 'w')
    s = r'\begin{tabular}{|lr|}' + '\n'
    s += r' \hline' + '\n'
    s += r' Uncertianty Source & Size [$10^{-2}$] \\' + '\n'
    s += r' \hline' + '\n'
    s += r' \hline' + '\n'

    uncSys = np.sqrt(uncRemaining[0]**2 - uncRemaining[-1]**2)
    s += r' Systematics & ' + '{:.2f}'.format(100*uncSys) + r' \\' + '\n'
    s += r' \hline' + '\n'
    fTable.write(s)

    for i in range(1, len(uncNames)):
        n = uncNames[i]
        s = ' '
        if n in supplementDic.keys():
            title = supplementDic[n]
        else:
            title = groups[i-2]['title']

        if n in ['stat']:
            fTable.write(r' \hline' + '\n')
        else:
            title = '\hspace{3mm} ' + title

        s += title + ' & '
        if uncAss[i-1] > 0.0001:
            s += '{:.2f}'.format(100*uncAss[i-1])
        else:
            s += '< 0.01'
        s += r' \\' + '\n'
        fTable.write(s)

    s = r' \hline' + '\n'
    s += r' Total & ' + '{:.2f}'.format(100*uncRemaining[0]) + r' \\' + '\n'
    s += r' \hline' + '\n'
    s += r'\end{tabular}' + '\n'
    fTable.write(s)
    fTable.close()

    return


########################### -------- Externalize parameters ------------------ #########################

def externalizeUncertainty(card, out, parameters=['B2DstCLNeig1', 'B2DstCLNeig2', 'B2DstCLNeig3'], center='preFit', sigma=1, tag='FF', rVal=SM_RDst, rLimits=[0.1, 0.7]):
    tag += '_center' + center

    if not out[-1] == '/':
        out += '/'
    out += 'externalization_' + tag
    if not os.path.isdir(out):
        os.makedirs(out)

    fLog = open(webFolder + '/externalize_'+tag+'.txt', 'w')
    s = 'Externalizing paramters:' + ' '.join(parameters) + '\nCenter: ' +  center
    print s
    fLog.write(s + '\n')

    if center == 'postFit':
        # Creating bestfit snapshop
        cmd = 'cd ' + out + '; '
        cmd += 'combine -M MultiDimFit'
        cmd += ' -d ' + card.replace('.txt', '.root')
        cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
        cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
        cmd += ' --setParameters r={:.2f}'.format(rVal)
        cmd += ' --setParameterRanges r={:.4f},{:.4f}'.format(*rLimits)
        cmd += ' -n Bestfit'
        cmd += ' --saveWorkspace --verbose 0'
        cmd += ' --trackParameters ' + ','.join(parameters)
        runCommandSafe(cmd)

        inputSpace = 'higgsCombineBestfit.MultiDimFit.mH120.root'
        arr = rtnp.root2array(out + '/higgsCombineBestfit.MultiDimFit.mH120.root', treename='limit')
        centralValue = {}
        for n in parameters:
            centralValue[n] = arr['trackedParam_'+n][0]
    elif center == 'preFit':
        inputSpace = card.replace('.txt', '.root')
        centralValue = dict.fromkeys(parameters, 0)
    else:
        print 'Center option', center, 'not recognized'

    for n in parameters:
        s = n + ' central value: {:.2f}'.format(centralValue[n])
        print s
        fLog.write(s + '\n')

    print '----- Running central fit'
    cmd = 'cd ' + out + '; '
    cmd += 'combine -M MultiDimFit --algo singles'
    cmd += ' -d ' + inputSpace
    cmd += ' --robustFit 1  --cminDefaultMinimizerStrategy=1 --X-rtd MINIMIZER_analytic'
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' --rMin={:.3f} --rMax={:.3f}'.format(*rLimits)
    cmd += ' --setParameters r={:.3f},'.format(rVal) + ','.join([n+'={:.2f}'.format(centralValue[n]) for n in parameters])
    cmd += ' --freezeParameters ' + ','.join(parameters)
    cmd += ' -n _ext'+tag+'_central'
    runCommandSafe(cmd)
    arr = rtnp.root2array(out + '/higgsCombine_ext'+tag+'_central.MultiDimFit.mH120.root', treename='limit')
    rCentral, d, u = arr['r']
    drCentral = (u-d)*0.5
    s = 'Central value {:.4f} +/- {:.4f}\n'.format(rCentral, drCentral)
    print s
    fLog.write(s + '\n')

    deltaR = np.zeros((len(parameters), 2))
    for ip, p in enumerate(parameters):
        rMod = []
        for mod in [sigma, -sigma]:
            newValueStr = p+'={:.2f}'.format(centralValue[p] + mod)
            cmdAux = cmd.replace(p+'={:.2f}'.format(centralValue[p]), newValueStr)
            cmdAux = cmdAux.replace(tag+'_central', tag+'_variation')
            runCommandSafe(cmdAux)
            arr = rtnp.root2array(out + '/higgsCombine_ext'+tag+'_variation.MultiDimFit.mH120.root', treename='limit')
            r, d, u = arr['r']
            dr = (u-d)*0.5
            s = '{} {:+.1f} ({:+.1f}): {:.4f} +/- {:.4f}'.format(p, mod, centralValue[p] + mod, r, dr)
            print s
            fLog.write(s + '\n')
            rMod.append(r)
        rMod = np.sort(rMod)
        s = p+' variation: {:.4f} {:+.4f}/{:+.4f}\n'.format(rCentral, rMod[1]-rCentral, rMod[0]-rCentral)
        print s
        fLog.write(s + '\n')
        deltaR[ip] = rMod

    extUnc = np.sqrt(np.sum(np.square(0.5*(deltaR[:,1] - deltaR[:,0]))))
    totUnc = np.hypot(drCentral, extUnc)

    s = 'Externalized result: {:.4f} +/- {:+.4f} (prof.) +/- {:+.4f} (ext.) = {:.4f} +/- {:+.4f}'.format(rCentral, drCentral, extUnc, rCentral, totUnc)
    print s
    fLog.write(s + '\n')
    fLog.close()
    return


########################### -------- Nuisances impact ------------------ #########################

def runNuisanceImpacts(card, out, catName, maskStr='', rVal=SM_RDst, submit=True, collect=True):
    if not out[-1] == '/': out += '/'

    if submit:
        if os.path.isdir(out+'impactPlots'):
            os.system('rm -rf '+out+'impactPlots')
        os.mkdir(out+'impactPlots')

        print '----- Running initial fit'
        cmd = 'cd {}impactPlots; '.format(out)
        cmd += ' combineTool.py -M Impacts --doInitialFit -m 0'
        cmd += ' --robustFit 1 --X-rtd MINIMIZER_analytic'
        cmd += ' --cminDefaultMinimizerStrategy=1'
        cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
        cmd += ' -d ' + card.replace('.txt', '.root')
        cmd += ' --setParameters r={:.2f}'.format(rVal)
        if maskStr:
            cmd += ','+maskStr
        cmd += ' --setParameterRanges r=0,0.6'
        cmd += ' --verbose 1'
        runCommandSafe(cmd)

        # If running on Tier2 condor remmeber to add this line to CombineToolBase.py ln 11
        # ``source /cvmfs/cms.cern.ch/cmsset_default.sh``
        print '----- Running all the fits'
        cmd = 'cd {}impactPlots;'.format(out)
        cmd += ' combineTool.py -M Impacts --doFits -m 0'
        cmd += ' --robustFit 1 --X-rtd MINIMIZER_analytic'
        cmd += ' --cminDefaultMinimizerStrategy=1'
        cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
        cmd += ' --parallel 100 --job-mode condor --task-name combineImpacts_'+os.path.basename(card).replace('.txt', '')
        cmd += ' --sub-opts "{}"'.format(stringJubCustomizationCaltechT2.replace('"', '\\\"').replace('$', '\$'))
        cmd += ' -d ' + card.replace('.txt', '.root')
        cmd += ' --setParameters r={:.2f}'.format(rVal)
        if maskStr:
            cmd += ','+maskStr
        cmd += ' --verbose -1'
        runCommandSafe(cmd)

    if collect:
        status, output = commands.getstatusoutput('condor_q')
        while 'combineImpacts_'+os.path.basename(card).replace('.txt', '') in output:
            time.sleep(30)
            status, output = commands.getstatusoutput('condor_q')
            for l in output.split('\n'):
                if 'combineImpacts_'+os.path.basename(card).replace('.txt', '') in l:
                    print l
                    sys.stdout.flush()
        cmd = 'cd {}impactPlots;'.format(out)
        cmd += ' combineTool.py -M Impacts -o impacts.json -m 0'
        cmd += ' -d ' + card.replace('.txt', '.root')
        runCommandSafe(cmd)


        rename = {
        'r': 'R(D*)',
        'mutauNorm': 'Norm B#rightarrow D*l#nu',
        'BdpT': 'B^{0} p_{T} spectrum',
        'BuPt': 'B^{+} p_{T} spectrum',
        'B2DstCLNR0':'R_{0} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig1':'#lambda_{1} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig2':'#lambda_{2} (CLN B#rightarrow D*l#nu)',
        'B2DstCLNeig3':'#lambda_{3} (CLN B#rightarrow D*l#nu)',
        'trgSF': 'Trigger scale factor',
        'trkEff': 'Tracking efficiency (control to signal region transfer factor)',
        'B2DstHcTransferFactor': 'Control to signal region transfer factor due to charged D decays',
        'softTrkEff': 'Soft tracks tracking efficiency',
        'softTrkEff_w': 'Soft tracks tracking efficiency (w)',
        'softTrkEff_s': 'Soft tracks tracking efficiency (s)',
        'tkPVfrac': 'Additional tracks origin',
        'overallNorm': 'Overall norm',
        'overallNormMu7_IP4': 'Overall norm (Low)',
        'overallNormMu9_IP6': 'Overall norm (Mid)',
        'overallNormMu12_IP6': 'Overall norm (High)',
        'fDststWide': 'D^{**} wide resonances fraction',
        'D2420_width': 'Width D**(2420)',
        'D2430_width': 'Width D**(2430)',
        'D2460_width': 'Width D**(2460)',
        'RDs_stst': 'R(D**_{(s)})',
        'Dst2S_width': 'Width D*(2S)',

        'brBd_DstDs' : 'Br. frac. B^{0}#rightarrow D*D_{s}',
        'brBd_DstDsst' : 'Br. frac. B^{0}#rightarrow D*D_{s}*',
        'brBd_DstDsst0' : 'Br. frac. B^{0}#rightarrow D*D_{s0}*',
        }
        for c in ['High', 'Mid', 'Low']:
            rename['tkPVfrac'+c] = 'Additional tracks origin ({})'.format(c)
            for i in range(1,5):
                s = str(i)
                rename['BdpT'+c+'_lam'+s] = 'B^{0} p_{T} '+c+' #lambda_{'+s+'}'
                rename['BuPt'+c+'_lam'+s] = 'B^{+} p_{T} '+c+' #lambda_{'+s+'}'

        procName_dic = {
        'mu'        : 'B^{0}#rightarrow D*#mu#nu',
        'tau'       : 'B^{0}#rightarrow D*#tau#nu',
        'DstPi'     : 'B#rightarrow D**(D*#pi)l#nu',
        'DstPiPi'   : 'B#rightarrow D**(D*#pi#pi)l#nu',
        'DstK'      : 'B_{s}#rightarrow D**_{s}(D*K)l#nu',
        'DuMu'      : 'D^{0}#rightarrow #muX',
        'DdMu'      : 'D^{+}#rightarrow #muX',
        'DsMu'      : 'D_{s}#rightarrow #muX',
        }
        for n in procName_dic: rename[n+'Br'] = 'Br. frac. ' + procName_dic[n]

        d = json.load(open(out+'impactPlots/impacts.json', 'r'))
        for par in d['params']:
            name = str(par['name'])
            if name.startswith('prop_bin'):
                label = name.replace('prop_bin', 'MC stat. ')
                label = label.replace('M2_miss_', 'M^{2}_{miss} ')
                label = label.replace('Est_mu_', 'E*_{#mu} ')
                label = label.replace('q2bin', '[b_{q^{2}}=')
                label = label.replace('_bin', '] ')
                rename[name] = label + 10*' '
            elif re.match(r'trgMu[0-9]+_IP[0-9]+SF_pt[0-9]+', name):
                label = 'Trigger scale factors ' + re.search(r'Mu[0-9]+_IP[0-9]', name).group(0)
                idx = name.find('_pt')
                label += ' bin ' + name[idx+3:]
                rename[name] = label
        json.dump(rename, open(out+'impactPlots/rename.json', 'w'))

        cmd = 'cd {};'.format(out)
        cmd += 'plotImpacts.py -i impactPlots/impacts.json -o impacts -t impactPlots/rename.json --max-pages 1'
        runCommandSafe(cmd)
        cmd = 'cp {}impacts.pdf {}/'.format(out, webFolder)
        runCommandSafe(cmd)

        cmd = 'cd {};'.format(out)
        cmd += 'plotImpacts.py -i impactPlots/impacts.json -o impacts_full -t impactPlots/rename.json'
        runCommandSafe(cmd)
        cmd = 'cp {}impacts_full.pdf {}/'.format(out, webFolder)
        runCommandSafe(cmd)

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
    cmd += ' --cminFallbackAlgo Minuit2,Migrad,0'
    cmd += ' -d ' + card.replace('.txt', '.root')
    if fixRDst:
        cmd += ' --freezeParameters r --setParameters r={:.3f}'.format(rVal)
    if maskEvalGoF:
        cmd += ' --setParametersForEval ' + maskEvalGoF
    cmd += ' -n Obs'+tag
    cmd += ' -t 0 -s 100'
    cmd += ' --verbose -1'
    runCommandSafe(cmd)
    arr = rtnp.root2array(gofOutdir+'/higgsCombineObs'+tag+'.GoodnessOfFit.mH120.100.root', treename='limit')
    s_obs = arr['limit'][0]
    print 'Observed test statistics: {:.1f}'.format(s_obs)


    # Run the test stat toy distribution
    cmdToys = cmd.replace('-n Obs', '-n Toys')
    nToysPerRep = 40 if args.runInJob else 10
    cmdToys = cmdToys.replace('-t 0 -s 100', '-t {} -s -1'.format(nToysPerRep))
    print cmdToys

    nRep = 10 if args.runInJob else 20
    p = Pool(min(20,nRep))
    outputs = p.map(runCommand, nRep*[cmdToys])
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
    nDigits = '3' if s_obs < 1 else '1'
    plt.plot([s_obs, s_obs], [0, np.max(content)], 'm--', label=('Obs: {:.'+nDigits+'f}\np-val {:.1f}%').format(s_obs, 100*p_val))
    plt.legend(loc='upper right')
    plt.xlabel('Test statistic (algo {})'.format(algo))
    plt.ylabel('Probability / {:.1f}'.format(0.5*(center[2]-center[1])))
    text = 'Cat. '+args.category.capitalize()
    if fixRDst:
        text += ' , R(D*) fixed'
    plt.text(0.04, 0.93, text, transform = plt.gca().transAxes)
    plt.savefig(out + 'fig/resultsGoF'+tag+'.png')
    plt.savefig(webFolder + '/resultsGoF'+tag+'.png')

    strRes = tag
    strRes += ' '*(45-len(strRes))
    strRes += '{:.3f} ({:.1f}%)'.format(s_obs, 100*p_val)
    strRes += ' '*(70-len(strRes))
    strRes += '{:.3f} {:.3f} {:.3f}'.format(np.percentile(s_toys, 50), np.percentile(s_toys, 95), np.percentile(s_toys, 99))
    print strRes
    with open(webFolder+'/GoF_results.txt', 'a') as dumpFile:
        dumpFile.write(strRes + '\n')
    os.system('echo "{}" >> {}GoF_results.txt'.format(strRes, out));


########################### -------- Condor submissions ------------------ #########################
# Check for the right singularity using: ll /cvmfs/singularity.opensciencegrid.org/cmssw/
mem = '7500'
if 'histos' in args.step:
    mem = '12000'
elif args.category == 'comb' and 'fitDiag' in args.step:
    mem = '15000'
jdlTemplate = '\n'.join([
              'executable        = '+basedir+'Combine/condorJob.sh',
              'arguments         = {arguments}',
              'output            = {outdir}/job_{jN}_$(ClusterId).out',
              'error             = {outdir}/job_{jN}_$(ClusterId).err',
              'log               = {outdir}/job_{jN}_$(ClusterId).log',
              'JobPrio           = -1',
              'WHEN_TO_TRANSFER_OUTPUT = ON_EXIT_OR_EVICT',
              '+MaxRuntime       = 1200',
              '+JobQueue         = ' + ('"Normal"' if ('fitDiag' in args.step) else '"Short"'),
              '+RunAsOwner       = True',
              '+InteractiveUser  = True',
              '+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"',
              '+SingularityBindCVMFS = True',
              'run_as_owner      = True',
              'RequestDisk       = 15000000', #KB
              'RequestMemory     = ' + mem, #MB
              'RequestCpus       = {:.0f}'.format(np.ceil(float(mem)/4000)),
              'x509userproxy     = $ENV(X509_USER_PROXY)',
              'on_exit_remove    = (ExitBySignal == False) && (ExitCode == 0)',
              'on_exit_hold      = (ExitBySignal == True) || (ExitCode != 0)',
              # 'periodic_release  =  (NumJobStarts < 2) && ((CurrentTime - EnteredCurrentStatus) > (60*5))',
              '+PeriodicRemove   = ((JobStatus =?= 2) && ((MemoryUsage =!= UNDEFINED && MemoryUsage > 5*RequestMemory)))',
              # 'max_retries       = 3',
              'requirements      = Machine =!= LastRemoteHost',
              'universe          = vanilla',
              'queue 1',''])

def submitRunToCondor():
    jobDir = outdir + '/condor'
    if not os.path.exists(jobDir):
        os.makedirs(jobDir)

    jN = 0
    for name in glob(jobDir + '/job_*.jdl'):
        n = os.path.basename(name)
        jN = max(jN, int(n[4:-4])+1)
    jN = str(jN)

    arguments = os.environ['PWD'] + ' '
    arguments += ' '.join(sys.argv).replace(' --submit', ' --runInJob')

    job = jdlTemplate.format(outdir=jobDir, arguments=arguments, jN=jN)
    with open(jobDir+'/job_'+jN+'.jdl', 'w') as fsub:
        fsub.write(job)

    cmd = 'cd '+jobDir+'; condor_submit job_'+jN+'.jdl'
    cmd += ' -batch-name ' + card_name + '_jN'+jN
    runCommandSafe(cmd)
    print 'Jobs submitted'


######################################################################################################


if __name__ == "__main__":
    if args.submit:
        submitRunToCondor()
        exit()

    if 'clean' in args.step:
        print '-----> Cleaning previous results'
        cleanPreviousResults()
        args.step.remove('clean')

    if 'histos' in args.step:
        if args.category == 'comb':
            print 'Histo should be created ahead running single categories'
        else:
            createHistograms(categoriesDef[args.category])
        args.step.remove('histos')

    if not args.step: exit()
    print '-----> Loading histograms'
    histo = None
    if args.category == 'comb':
        histo = {}
        for c in categoriesToCombine:
            print '---- Loading', c
            present = False
            while not present:
                n = len(glob(os.path.join(histo_file_dir, card_name.replace('comb', c)) + '_*.root'))
                if n>50:
                    present = True
                else:
                    print 'Waiting for ' + c
                    time.sleep(60*10)
            histo[c] = loadHisto4CombineFromRoot(histo_file_dir, card_name.replace('comb', c))
    else:
        loadShapeVar = 'card' in args.step
        histo = loadHisto4CombineFromRoot(histo_file_dir, card_name, loadShapeVar=loadShapeVar, verbose=False)

    if 'preFitPlots' in args.step:
        if args.category == 'comb':
            cPre = {}
            for c in categoriesToCombine:
                cPre[c] = drawPlots('prefit'+c.capitalize(), histo[c], c.capitalize(), scale_dic={'tau': SM_RDst})
        else:
            cPre = drawPlots('prefit', histo, args.category.capitalize(), scale_dic={'tau': SM_RDst})
        args.step.remove('preFitPlots')

    if 'shapeVarPlots' in args.step:
        if args.category == 'comb':
            cPre = {}
            for c in categoriesToCombine:
                drawShapeVarPlots(card=card_name.replace('comb', c), tag='_'+c)
        else:
            drawShapeVarPlots(card_name)
        args.step.remove('shapeVarPlots')

    if 'card' in args.step:
        for fitRegionsOnly in [False, True]:
            print '-----> Creating the datacard' + (' (fit regions only)' if fitRegionsOnly else '')
            cl = card_location.replace('.txt', '_fitRegionsOnly.txt') if fitRegionsOnly else card_location
            if args.category == 'comb':
                createCombinationCard(fitRegionsOnly)
            else:
                card = createSingleCard(histo, categoriesDef[args.category], fitRegionsOnly)
                if args.showCard:
                    print 3*'\n'
                    print '--------> Dumping Combine card (fitRegionsOnly {})'.format(fitRegionsOnly)
                    print card
                    print 3*'\n'
                print 'Card location:', cl
                fc = open(cl, 'w')
                fc.write(card)
                fc.close()

            cmd = 'cp '+cl+' '+webFolder+'/'+os.path.basename(cl)
            runCommandSafe(cmd)

        if args.validateCard:
            cl = card_location.replace('.txt', '_fitRegionsOnly.txt')
            cmd = 'cd ' + os.path.dirname(cl) + '; '
            cmd += 'ValidateDatacards.py ' + os.path.basename(cl) + ' -p 3 -c 0.2'
            print cmd
            status, output = commands.getstatusoutput(cmd)
            if ('ERROR' in output) or ('There were  ' in output):
                print output

        print '\n'
        args.step.remove('card')

    if not args.step: exit()

    if 'workspace' in args.step:
        createWorkspace(card_location)
        createWorkspace(card_location.replace('.txt', '_fitRegionsOnly.txt'))
        print '\n'

    if 'bias' in args.step:
        biasOut = outdir + '/biasStudyToys'
        if not os.path.isdir(biasOut):
            os.makedirs(biasOut)

        both = args.runBiasAnalysis and args.runBiasToys
        if (not args.runBiasAnalysis) or both:
            biasToysScan(card_location.replace('.txt', '_fitRegionsOnly.txt'), biasOut, args.seed, args.nToys)
        if (not args.runBiasToys) or both:
            collectBiasToysResults(biasOut, args.toysRDst)

    if 'scan' in args.step:
        print '-----> Running likelyhood scan'
        if len(args.RDstLims) > 0:
            rLimits = args.RDstLims
        elif args.unblinded:
            rLimits = [0.08]
        elif not args.unblinded:
            if args.category == 'comb':
                rLimits = [0.15]
            else:
                rLimits = [0.25]

        if args.maskScan:
            maskList = []
            for n in args.maskScan:
                knList = histo.keys()
                if args.category == 'comb':
                    knList = []
                    for c, hDic in histo.iteritems():
                        for n in hDic.keys():
                            knList.append(c+'_'+n)
                for kn in knList:
                    if not re.match(n, kn) is None:
                        print 'Masking', kn
                        maskList.append('mask_'+kn+'=1')
            maskStr = ','.join(maskList)
            fit_RDst, rDst_postFitRegion = runScan(args.scanTag, card_location, outdir,
                                                   args.category.capitalize(),
                                                   maskStr=maskStr,
                                                   rLimits=rLimits, strategy=0, draw=True, dumpNuis=True)
        elif args.scanTag:
            fit_RDst, rDst_postFitRegion = runScan(args.scanTag, card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                                                   args.category.capitalize(),
                                                   freezePars=args.freezeParsScan,
                                                   rLimits=rLimits,
                                                   strategy=0, draw=True, dumpNuis=True)
        else:
            fit_RDst, rDst_postFitRegion = runScan('Base', card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                                                   args.category.capitalize(),
                                                   rLimits=rLimits,
                                                   strategy=0, draw=True, dumpNuis=True)
    else:
        rDst_postFitRegion = args.RDstLims if len(args.RDstLims) == 2 else [0.15, 0.45]
        fit_RDst = SM_RDst
    rDst_postFitRegion[0] = max(0, rDst_postFitRegion[0])
    print '[INFO] R(D*) central value set to {:.3f}'.format(fit_RDst)
    print '[INFO] R(D*) boundaries set to [{:.3f}, {:.3f}]'.format(*rDst_postFitRegion)

    if 'catComp' in args.step:
        if not args.category == 'comb':
            print '[WARNING]: Skipping catComp because not defined for category='+args.category
        else:
            print '-----> Running categories compatibility'
            categoriesCompatibility(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                                    rVal=SM_RDst, rLimits=rDst_postFitRegion
                                    )

    if 'uncBreakdownTable' in args.step:
        runUncertaintyBreakDownTable(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                                args.category.capitalize(),
                                rVal=fit_RDst, rLimits=rDst_postFitRegion)

    if 'externalize' in args.step:
        print '-----> Running externalization'
        externalizeUncertainty(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                               parameters=args.externPars,
                               center=args.externCenter,
                               sigma=args.externSigma,
                               tag=args.externTag,
                               rVal=SM_RDst, rLimits=[0.1, 0.45]
                               )

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

    if 'fitDiag' in args.step:
        print '-----> Running fit diagnostic'
        globalChannelMaskingStr = defineChannelMasking(histo)
        runFitDiagnostic(args.cardTag, card_location, outdir,
                         strategy=0 if args.category == 'comb' else 1,
                         forceRDst=args.forceRDst, maskStr=globalChannelMaskingStr,
                         rVal=fit_RDst, rLimits=rDst_postFitRegion)

    if 'postFitPlots' in args.step:
        print '-----> Getting postfit results'
        histo_post, _, _ = getPostfitHistos(args.cardTag, outdir, forceRDst=args.forceRDst, histo_prefit=histo)
        if args.category == 'comb':
            cPost = {}
            cPrePostComp = {}
            for c in categoriesToCombine:
                tag = 'postfit'+c.capitalize() + ('_RDstFixed' if args.forceRDst else '')
                cPost[c] = drawPlots(tag, histo_post[c], c.capitalize())
                cPrePostComp[c] = drawPrePostFitComparison(histo[c], histo_post[c], tag=c + ('_RDstFixed' if args.forceRDst else ''))
        else:
            tag = 'postfit'+ ('_RDstFixed' if args.forceRDst else '')
            cPost = drawPlots(tag, histo_post, args.category.capitalize())
            cPrePostComp = drawPrePostFitComparison(histo, histo_post, tag=('_RDstFixed' if args.forceRDst else ''))
        nuisancesDiff(args.cardTag, outdir, args.forceRDst)
        print '\n'

    if 'uncBreakdownScan' in args.step:
        runUncertaintyBreakDownScan(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                                args.category.capitalize(),
                                rVal=fit_RDst, rLimits=rDst_postFitRegion)

    if 'impacts' in args.step:
        print '-----> Running impact plots'
        submit, collect = True, True
        if args.collectImpacts and not args.subOnlyImpacts:
            submit, collect = False, True
        elif not args.collectImpacts and args.subOnlyImpacts:
            submit, collect = True, False
        runNuisanceImpacts(card_location.replace('.txt', '_fitRegionsOnly.txt'), outdir,
                           args.category.capitalize(),
                           rVal=fit_RDst, submit=submit, collect=collect)
