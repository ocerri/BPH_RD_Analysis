#!/usr/bin/env python

import sys, os, pickle, yaml, time, argparse
sys.path.append('../lib')
if os.environ['CMSSW_VERSION'] != 'CMSSW_10_2_3':
    raise

from glob import glob
from array import array

import numpy as np
import scipy.stats as sps
import pandas as pd

import matplotlib.pyplot as plt
import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp


import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "     Internal"

from histo_utilities import create_TH1D, create_TH2D, std_color_list, SetMaxToMaxHist, make_ratio_plot
from progressBar import ProgressBar
from lumi_utilities import getLumiByTrigger
from pileup_utilities import pileupReweighter
from beamSpot_calibration import getBeamSpotCorrectionWeights

from categoriesDef import categories
from analysis_utilities import drawOnCMSCanvas, DSetLoader
from pT_calibration_reader import pTCalReader as calibrationReader

parser = argparse.ArgumentParser(description='Script used to calibrate Bd kin',
                                 epilog='Test example: ./kinematicCalibration_Bd_JpsiKst.py',
                                 add_help=True
                                 )
parser.add_argument ('--category', '-c', type=str, default='high', choices=['low', 'mid', 'high'], help='Category.')
parser.add_argument ('--version', '-v', default='test', help='Version.')
parser.add_argument ('--showPlots', default=False, action='store_true', help='Show plots by setting ROOT batch mode OFF (default ON)')
parser.add_argument ('--draw_precal', default=False, action='store_true', help='Draw also precal plots')
parser.add_argument ('--verbose', default=False, action='store_true', help='Verbose')
args = parser.parse_args()

print '\n'+40*'#'
print 'Kinematic calibration of Bd -> JpsiK*'
print 'Category:', args.category
print 40*'#'

# Example with loop: for cat in "low" "mid" "high"; do ./kinematicCalibration_Bd_JpsiKst.py -c $cat; done;

if not args.showPlots:
    rt.gROOT.SetBatch(True)
    plt.ioff()
    plt.switch_backend('Agg')

cat = categories[args.category]
version = args.version

mB_var = 'mass_candB'
m_B0 = 5.27963 #1e-3*Particle.from_string('B0').mass

cl = rt.TLine()
cl.SetLineColor(6)
cl.SetLineStyle(9)
cl.SetLineWidth(2)

catText = rt.TLatex()
catText.SetTextAlign(31)
catText.SetTextSize(0.06)
catText.SetTextFont(42)
catText.SetTextSize(0.05)


webFolder = os.environ['HOME']+'/public_html/BPH_RDst/kinematicCalibration_Bd_JpsiKst/'+version+'/'
if not os.path.isdir(webFolder):
    os.makedirs(webFolder)
    os.system('cp {v}../index.php {v}'.format(v=webFolder))

dataLoc = '/storage/af/group/rdst_analysis/BPhysics/data/'

def getPolyCorrection(hNum, hDen, deg, tag, verbose=False):
    y = []
    yerr = []
    x = []
    hAux = hNum.Clone('hAux')
    hAux.Divide(hDen)
    for i in range(1, hNum.GetNbinsX()+1):
        if np.abs(hAux.GetBinContent(i))<1e-3:
            continue
        x.append(hAux.GetBinCenter(i))
        y.append(hAux.GetBinContent(i))
        yerr.append(hAux.GetBinError(i))

    beta, covBeta = np.polyfit(x, y, deg=deg, full=False, w=1./np.array(yerr), cov='unscaled')
    eigVal, eigVec = np.linalg.eig(covBeta)
    eigSig = np.sqrt(eigVal)

    if verbose:
        print tag
        print 'Beta:  '+', '.join(beta.shape[0]*['{:1.2e}']).format(*beta)
        print 'Error: '+', '.join(beta.shape[0]*['{:1.2e}']).format(*np.sqrt(np.diag(covBeta)))

    betaVar = []
    for i in range(eigSig.shape[0]):
        betaVar.append(eigSig[i]*eigVec[:, i])
        if verbose:
            print '\n'
            print 'eigSigma: {:1.2e}'.format(eigSig[i])
            print 'eigVect: '+', '.join(beta.shape[0]*['{:.2f}']).format(*eigVec[:, i])
            print 'Variation: '+', '.join(beta.shape[0]*['{:1.2e}']).format(*betaVar[i])

    plt.figure(figsize=(8,6))
    plt.errorbar(x, y, yerr=yerr, fmt='k.', label='data/MC')
    plt.plot(x, np.polyval(beta, x), '-', label='Central')
    colors = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    for i in range(len(betaVar)):
        if i > 2:
            break
        yP = np.polyval(beta+betaVar[i], x)
        yM = np.polyval(beta-betaVar[i], x)
        plt.plot(x, yP, '--', color=colors[i], label='w, $\lambda_{} \pm 1\sigma$'.format(i))
        plt.plot(x, yM, '--', color=colors[i])

    plt.grid()
    ymin, ymax = plt.ylim()
    plt.ylim(max(0, ymin), min(8, ymax))
    plt.ylabel('data/MC')
    plt.xlabel('B '+tag.split('_')[0])
    plt.legend(loc='best', numpoints=1)
    plt.savefig(webFolder+tag+'_' + cat.name + '.png', bbox_inches='tight')

    dOut = {'beta': beta, 'betaVar' : betaVar}
    pickle.dump(dOut, open(dataLoc+'calibration/kinematicCalibration_Bd/'+tag+'_{}_{}.pkl'.format(cat.name, version), 'wb'))


# # Load MC
mcSample = DSetLoader('Bd_JpsiKst_General', candDir='ntuples_Bd2JpsiKst_220228')
dsetMC_loc = mcSample.skimmed_dir + '/{}_corr.root'.format(cat.name)
# dsetMC_loc = mcSample.skimmed_dir + '/{}_bare.root'.format(cat.name)
dfMC = pd.DataFrame(rtnp.root2array(dsetMC_loc))


effMCgen = mcSample.effMCgen
brFileLoc = dataLoc+'forcedDecayChannelsFactors_v2.pickle'
decayBR = pickle.load(open(brFileLoc, 'rb'))['JPsiKst']
effCAND = mcSample.effCand['effCAND']
effSkim = mcSample.getSkimEff(cat.name+'_bare')

xsec_eff = 1
dxsec = 0
for f, df in [effMCgen['xsec'], effMCgen['effGEN'], decayBR, effCAND, effSkim]:
    xsec_eff *= f
    dxsec += np.square(df/f)
dxsec = xsec_eff * np.sqrt(dxsec)
print 'Expected evts/fb: {:.0f} +/- {:.0f}'.format(xsec_eff, dxsec)

puRew = pileupReweighter(dsetMC_loc, 'hAllNTrueIntMC', trg=cat.trg)
dfMC['wPU'] = puRew.getPileupWeights(dfMC['MC_nInteractions'])

beamSpotCalLoc = '/storage/af/group/rdst_analysis/BPhysics/data/calibration/beamSpot/crystalball_calibration_v2_bs_'+args.category.capitalize()+'.yaml'
paramBeamSpotCorr = yaml.load(open(beamSpotCalLoc, 'r'))
dfMC['wBeamSpot'] = getBeamSpotCorrectionWeights(dfMC, paramBeamSpotCorr, ref='bs')

loc = dataLoc+'calibration/triggerScaleFactors/'
fTriggerSF = rt.TFile.Open(loc + 'HLT_' + cat.trg + '_SF_v22_count.root', 'READ')
hTriggerSF = fTriggerSF.Get('hSF_HLT_' + cat.trg)

ptmax = hTriggerSF.GetXaxis().GetXmax() - 0.01
ipmax = hTriggerSF.GetYaxis().GetXmax() - 0.01
etamax = hTriggerSF.GetZaxis().GetXmax() - 0.01

nX = hTriggerSF.GetNbinsX()
ptWeight = [[] for k in range(nX+2)]

dfMC['trgSF'] = np.ones(dfMC.shape[0])
for i, (pt, eta, ip) in enumerate(dfMC[['trgMu_pt', 'trgMu_eta', 'trgMu_sigdxy']].values):
    ix = hTriggerSF.GetXaxis().FindBin(min(ptmax, pt))
    iy = hTriggerSF.GetYaxis().FindBin(min(ipmax, ip))
    iz = hTriggerSF.GetZaxis().FindBin(min(etamax, np.abs(eta)))
    dfMC.at[i, 'trgSF'] = hTriggerSF.GetBinContent(ix, iy, iz)
#     if np.abs(dfMC.at[i, 'trgSF'] - 1) > 0.1:
#         print (4*'{:.2f} ').format(pt, eta, ip, hTriggerSF.GetBinContent(ix, iy, iz))
    ptWeight[ix].append(hTriggerSF.GetBinContent(ix, iy, iz))

# Muon ID scale factor
loc = dataLoc+'calibration/muonIDscaleFactors/Run2018ABCD_SF_MuonID_Jpsi.root'
fMuonIDSF = rt.TFile.Open(loc, 'READ')
hMuonIDSF = fMuonIDSF.Get('NUM_MediumID_DEN_genTracks_pt_abseta')

dfMC['muonSF'] = np.ones(dfMC.shape[0])
for i, (ptp, etap, ptm, etam) in enumerate(dfMC[['MC_mup_pt', 'MC_mup_eta', 'MC_mum_pt', 'MC_mum_eta']].values):
    ix = hMuonIDSF.GetXaxis().FindBin(min(39.9,ptp))
    iy = hMuonIDSF.GetYaxis().FindBin(np.abs(etap))
    wp = hMuonIDSF.GetBinContent(ix, iy)
    ix = hMuonIDSF.GetXaxis().FindBin(min(39.9,ptm))
    iy = hMuonIDSF.GetYaxis().FindBin(np.abs(etam))
    wm = hMuonIDSF.GetBinContent(ix, iy)
    dfMC.at[i, 'muonSF'] = wp * wm

dfMC['w'] = dfMC['wPU']*dfMC['wBeamSpot']*dfMC['muonSF']*dfMC['trgSF']


dpt_rel = np.abs(dfMC['B_pt']/dfMC['MC_B_pt'] - 1)
deta = np.abs(dfMC['B_eta'] - dfMC['MC_B_eta'])
dphi = np.abs(dfMC['B_phi'] - dfMC['MC_B_phi'])
dphi = np.abs(dphi - 2*np.pi*(dphi > np.pi).astype(np.float))

safe = np.logical_and(np.logical_and(dpt_rel < 0.1, deta < 0.05), dphi < 0.05)
print 'MC purity (kin): {:.1f}%'.format(100*np.sum(safe)/float(safe.shape[0]))
print 'MC purity (idx match): {:.1f}%'.format(100*np.sum(dfMC['MC_idxMatch'] == 1)/float(dfMC['MC_idxMatch'].shape[0]))


# # Load data
# datasets_loc = glob(dataLoc + 'cmsRD/ParkingBPH*/*2018*B2JpsiKst_210501*')
# lumi_tot = getLumiByTrigger(datasets_loc, cat.trg, verbose=True)
# if not lumi_tot:
expectedLumi = {'Low':6.4, 'Mid':20., 'High':26.} #fb^-1
lumi_tot = expectedLumi[cat.name]
print 'Total lumi (estimated): {:.1f} fb^-1'.format(lumi_tot)
CMS_lumi.integrated_lumi = lumi_tot



dsetRD_loc = dataLoc+'cmsRD/skimmed/B2JpsiKst_220228_{}_corr.root'.format(cat.name)
dfRD = pd.DataFrame(rtnp.root2array(dsetRD_loc))
N_sel_per_fb = float(dfRD.shape[0])/lumi_tot
print 'Selected events per fb: {:.0f}'.format(N_sel_per_fb)

# # Clean sets
cuts = [
    ['B_eta', [-0.8, 0.8]],
    [mB_var, [5.24, 5.34]],
]

fout = open(webFolder + 'additionalSelection.txt', 'w')
print 'Including cuts:'
for i, d in enumerate([dfMC, dfRD]):
    sel = np.ones_like(d['mum_pt']).astype(np.bool)
    for v, [ll, ul] in cuts:
        sel = np.logical_and(sel, np.logical_and(d[v] > ll, d[v] < ul))
        if i == 0:
            fout.write('{:.3f} < {} < {:.3f}'.format(ll, v, ul))
            print '{:.3f} < {} < {:.3f}'.format(ll, v, ul)
    if i == 0:
        dfMC = d[sel]
    else:
        dfRD = d[sel]
fout.close()

# # Compare data/MC
def makePlot(var, binning=None, axis_title=None, wMC=[None], wRD=None, tag='', legMC=['MC'], legRD='Data', saveFig=''):
    h = create_TH1D(dfRD[var], name='hRD', axis_title=axis_title,
                    weights=wRD, binning=binning, scale_histo='norm')
    h.SetMarkerStyle(15)
    hList = [h]
    leg = rt.TLegend(0.8, 0.7, 0.95,0.8)
    leg.SetBorderSize(0)
    leg.AddEntry(h, legRD, 'lep')

    for i, (w, l) in enumerate(zip(wMC, legMC)):
        h = create_TH1D(dfMC[var], name='hMC'+str(i), binning=hList[0].binning,
                        weights=w, scale_histo='norm', color=i)

        leg.AddEntry(h, l, 'le')
        hList.append(h)

    CMS_lumi.integrated_lumi = lumi_tot
    m = SetMaxToMaxHist(hList)
    c = drawOnCMSCanvas(CMS_lumi, hList, 'same')
    leg.Draw()
    c.leg = leg
    catText.DrawLatexNDC(0.95, 0.85, 'Category: {}'.format(cat.name))
    if saveFig:
        c.SaveAs(webFolder+saveFig)

makePlot('N_vtx', binning=[70, 0.5, 70.5], wMC=[dfMC['w']],
         axis_title=['Number of vertexes', 'Normalized entries'], saveFig='reconstructedVertexes_'+cat.name+'.png')

if args.draw_precal:
    makePlot(mB_var, binning=[75, 5.24, 5.34], wMC=[dfMC['w']],
             axis_title=['mass(#mu#mu#piK)', 'Normalized entries'], saveFig='mass_mumupiK_precal_'+cat.name+'.png')
    makePlot('mass_mumu', binning=[70, 3.0, 3.2], wMC=[dfMC['w']],
             axis_title=['mass(#mu#mu)', 'Normalized entries'], saveFig='mass_mumu_precal_'+cat.name+'.png')
    makePlot('mass_candKst', binning=[70, 0.82, 0.97], wMC=[dfMC['w']],
             axis_title=['mass(#piK)', 'Normalized entries'], saveFig='mass_Kst_precal_'+cat.name+'.png')


    b = {'Low': array('d', list(np.arange(7, 9.4, 0.1)) + [9.4]),
         'Mid': array('d', list(np.arange(9, 12.2, 0.1)) +[12.2]),
         'High': array('d', list(10+np.logspace(np.log10(12.2-10), np.log10(50), 30)))
        }
    makePlot('trgMu_pt', binning=b[cat.name],  wMC=[dfMC['w']],
             axis_title=['trigger muon p_{T}', 'Normalized entries'], saveFig='trgMu_pt_precal_'+cat.name+'.png')

    for n in ['otherMu', 'K', 'pi', 'Jpsi', 'Kst']:
        makePlot(n+'_pt', binning=3*[None], axis_title=[n+' p_{T}', 'Normalized entries'],
                  wMC=[dfMC['w']], saveFig=n+'_pt_precal_'+cat.name+'.png')


    for n in ['trgMu', 'otherMu', 'K', 'pi', 'Jpsi', 'Kst']:
        makePlot(n+'_eta', binning=[50, -2, 2], axis_title=[n+' #eta', 'Normalized entries'],
                  wMC=[dfMC['w']], saveFig=n+'_eta_precal_'+cat.name+'.png')


# # Eta Calibration

b=[30,-1.8,1.8]
binWdith = (b[2] - b[1])/float(b[0])
hRD = create_TH1D(dfRD['B_eta'], name='hRD', title='data',
                  axis_title=['B #eta (reco)',
                              '1/#sigma d#sigma/d#eta / '+'({:.2f})'.format(binWdith)],
                  binning=b, scale_histo='norm', opt='overflow+underflow'
                 )
hRD.SetMarkerStyle(15)

hMC = create_TH1D(dfMC['B_eta'], name='hMC', weights=dfMC['w'],
                  axis_title=['B #eta (reco)',
                              '1/#sigma d#sigma/d#eta / '+'({:.2f})'.format(binWdith)],
                  title = 'B#rightarrow J/#psi K*',
                  scale_histo='norm', color=0,
                  h2clone=hRD, opt='overflow+underflow')

CMS_lumi.extraText = '      Internal'
c = make_ratio_plot([hMC, hRD], ratio_bounds=[0.8, 1.2], draw_opt='E1')

CMS_lumi.CMS_lumi(c, -1, 0)
c.pad1.SetTopMargin(0.07)
c.pad1.SetRightMargin(0.035)
c.pad2.SetRightMargin(0.035)
# c.pad2.SetLogy()
c.leg.SetY1(0.3)
c.leg.SetY2(0.5)
c.leg.SetX1(0.35)
c.leg.SetX2(0.7)
c.Draw()

c.pad1.cd()
catText.SetTextSize(0.04)
catText.DrawLatexNDC(0.9, 0.8, 'Category: {}'.format(cat.name))
c.SaveAs(webFolder+'B_eta_' + cat.name + '_precal.png')

getPolyCorrection(hRD, hMC, 5, 'eta_polyCoeff', args.verbose)
cal_eta = calibrationReader(calibration_file=dataLoc+'calibration/kinematicCalibration_Bd/eta_polyCoeff_{}_{}.pkl'.format(cat.name, version))
dfMC['wEta'] = cal_eta.getWeights(dfMC['B_eta'])



b=[30,-1.8,1.8]
binWdith = (b[2] - b[1])/float(b[0])
hRD = create_TH1D(dfRD['B_eta'], name='hRD',
                  title='data',
                  axis_title=['B #eta (reco)',
                              '1/#sigma d#sigma/d#eta / '+'({:.2f})'.format(binWdith)],
                  binning=b,
                  scale_histo='norm',
                  opt='overflow+underflow'
                 )
hRD.SetMarkerStyle(15)

hMC1 = create_TH1D(dfMC['B_eta'], name='hMC',
                  weights=dfMC['w'],
                  axis_title=['B #eta (reco)',
                              '1/#sigma d#sigma/d#eta / '+'({:.2f})'.format(binWdith)],
                  title = 'B#rightarrow J/#psi K*',
                  scale_histo='norm', color=0,
                  h2clone=hRD, opt='overflow+underflow')

hMC2 = create_TH1D(dfMC['B_eta'], name='hMC',
                  weights=dfMC['w']*dfMC['wEta'],
                  axis_title=['B #eta (reco)',
                              '1/#sigma d#sigma/d#eta / '+'({:.2f})'.format(binWdith)],
                  title = 'B#rightarrow J/#psi K* (#eta corr)',
                  scale_histo='norm', color=1,
                  h2clone=hRD, opt='overflow+underflow')


CMS_lumi.extraText = '      Internal'
c = make_ratio_plot([hMC2, hMC1, hRD], ratio_bounds=[0.8, 1.2], draw_opt='E1')

CMS_lumi.CMS_lumi(c, -1, 0)
c.pad1.SetTopMargin(0.07)
c.pad1.SetRightMargin(0.035)
c.pad2.SetRightMargin(0.035)
# c.pad2.SetLogy()
c.leg.SetY1(0.3)
c.leg.SetY2(0.5)
c.leg.SetX1(0.35)
c.leg.SetX2(0.7)
c.Draw()

c.pad1.cd()
catText.SetTextSize(0.04)
catText.DrawLatexNDC(0.9, 0.8, 'Category: {}'.format(cat.name))
c.SaveAs(webFolder+'B_eta_' + cat.name + '_etaOnlycal.png')


# # pT Calibration


lowOff = 5
midOff = 10
highOff = 5
b = {'Low': array('d', list( lowOff+np.logspace(np.log10(12-lowOff), np.log10(90-lowOff), 40) )),
     'Mid': array('d', list( midOff+np.logspace(np.log10(14-midOff), np.log10(100-midOff), 45) )),
     'High': array('d', list( highOff+np.logspace(np.log10(18-highOff), np.log10(125-highOff), 50) )),
    }

binWdith = b[cat.name][1] - b[cat.name][0]
hRD = create_TH1D(dfRD['B_pt'], name='hRD',
                  title='data',
                  axis_title=['B p_{T} (reco) [GeV]',
                              '1/#sigma d#sigma/dp_{T}'],
                  binning=b[cat.name],
                  scale_histo='norm',
                  widthNorm=True,
                  opt='overflow+underflow'
                 )
hRD.SetMarkerStyle(15)

hMCb = create_TH1D(dfMC['B_pt'], name='hMCb',
                  title = 'B#rightarrow J/#psi K* (no #eta cal)',
                  weights=dfMC['w'],
                  scale_histo='norm', color=1,
                  widthNorm=True,
                  binning=hRD.binning, opt='overflow+underflow')

hMC = create_TH1D(dfMC['B_pt'], name='hMC',
                  weights=dfMC['w']*dfMC['wEta'],
                  title = 'B#rightarrow J/#psi K*',
                  scale_histo='norm', color=0,
                  widthNorm=True,
                  binning=hRD.binning, opt='overflow+underflow')

CMS_lumi.extraText = '      Internal'
cr = make_ratio_plot([hMC, hRD, hMCb], ratio_bounds=[0.5, 2.5], draw_opt='E1')
# c = make_ratio_plot([hRD, hMCb, hMC], ratio_bounds=[0.5, 10], draw_opt='E1')
CMS_lumi.CMS_lumi(cr, -1, 0)
cr.pad1.SetTopMargin(0.07)
cr.pad1.SetRightMargin(0.035)
cr.pad2.SetRightMargin(0.035)
cr.pad2.SetLogy()
cr.leg.SetY2(0.9)
cr.leg.SetY1(0.6)
cr.leg.SetX1(0.6)
cr.Draw()

cr.pad1.cd()
catText.SetTextSize(0.04)
catText.DrawLatexNDC(0.9, 0.5, 'Category: {}'.format(cat.name))

cr.SaveAs(webFolder+'BpT_preCal_' + cat.name + '.png')


getPolyCorrection(hRD, hMC, 3, 'pt_polyCoeff', args.verbose)
cal_pT = calibrationReader(calibration_file=dataLoc+'calibration/kinematicCalibration_Bd/pt_polyCoeff_{}_{}.pkl'.format(cat.name, version))
dfMC['wPt'] = cal_pT.getWeights(dfMC['MC_B_pt'])



h_var = {}
h_var['C'] = create_TH1D(dfMC['B_pt'], name='h_var_C', binning=hRD.binning, opt='underflow+overflow',
                         weights=dfMC['w']*dfMC['wEta']*dfMC['wPt'])
norm = float(h_var['C'].Integral())
h_var['C'].Scale(1./norm, 'width')

mVar = np.array(cal_pT.betaVar)
nIter = 1000
binContent = np.zeros((nIter, hRD.GetNbinsX()))
pb = ProgressBar(nIter)
for itx in range(nIter):
    pb.show(itx)
    p = cal_pT.beta + np.sum(np.random.normal(size=(cal_pT.beta.shape[0],1))*mVar, axis=0)

    wPt = np.polyval(p, dfMC['MC_B_pt'])
    hAux = create_TH1D(dfMC['B_pt'], name='hAux', binning=hRD.binning, opt='underflow+overflow',
                       weights=dfMC['w']*dfMC['wEta']*wPt)
    hAux.Scale(1./norm, 'width')

    for i in range(1, binContent.shape[1]+1):
        binContent[itx, i-1] = hAux.GetBinContent(i)

for n, var in {'Down':-1, 'Up':1}.iteritems():
    h_var[n] = h_var['C'].Clone('h_var_'+n)
    content = np.percentile(binContent, q=100*sps.norm.cdf(var), axis=0)
    for i in range(1, content.shape[0]+1):
        h_var[n].SetBinContent(i, content[i-1])



gr_stat = rt.TGraphAsymmErrors()
gr_sys = rt.TGraphAsymmErrors()
h_dr = hRD.Clone('h_aux_dataratio')
h_mr = hMC.Clone('h_aux_MCratio')
g_up = rt.TGraph()
g_up.SetPoint(0, hMC.GetBinCenter(1)-0.5*hMC.GetBinWidth(1), 1)
g_down = rt.TGraph()
g_down.SetPoint(0, hMC.GetBinCenter(1)-0.5*hMC.GetBinWidth(1), 1)
for ib in range(1, hRD.GetNbinsX()+1):
    x = hRD.GetBinCenter(ib)
    y = h_var['C'].GetBinContent(ib)
    c = h_dr.GetBinContent(ib)
    e = h_dr.GetBinError(ib)
    if y == 0:
        y = 1e-6
    h_dr.SetBinContent(ib, c/y)
    h_dr.SetBinError(ib, e/y)
    c = h_mr.GetBinContent(ib)
    e = h_mr.GetBinError(ib)
    h_mr.SetBinContent(ib, c/y)
    h_mr.SetBinError(ib, e/y)
    gr_stat.SetPoint(ib-1, x, y)
    dx = 0.5*hMC.GetBinWidth(ib)
    dy = h_var['C'].GetBinError(ib)
    gr_stat.SetPointError(ib-1, dx, dx, dy, dy)

    dy_low = max(y-h_var['Up'].GetBinContent(ib), y-h_var['Down'].GetBinContent(ib))
    dy_up = max(h_var['Up'].GetBinContent(ib)-y, h_var['Down'].GetBinContent(ib)-y)
    gr_sys.SetPoint(ib-1, x, y)
    gr_sys.SetPointError(ib-1, dx, dx, dy_low, dy_up)

    x_low = h_dr.GetBinCenter(ib) - 0.5*h_dr.GetBinWidth(ib)
    x_up = h_dr.GetBinCenter(ib) + 0.5*h_dr.GetBinWidth(ib)
    g_up.SetPoint(2*ib-1, x_low, (y+dy_up)/y)
    g_up.SetPoint(2*ib, x_up, (y+dy_up)/y)
    g_down.SetPoint(2*ib-1, x_low, (y-dy_low)/y)
    g_down.SetPoint(2*ib, x_up, (y-dy_low)/y)
g_up.SetPoint(2*ib+1, x_up, 1)
g_down.SetPoint(2*ib+1, x_up, 1)
gr_stat.SetLineColor(rt.kRed-4)
gr_stat.SetLineWidth(2)
gr_stat.SetMarkerColor(rt.kRed-4)
gr_sys.SetFillColor(rt.kRed-4)
gr_sys.SetFillStyle(3005)
gr_sys.SetLineWidth(0)


leg = rt.TLegend(0.5, 0.3, 0.93,0.55)
leg.SetBorderSize(0)
leg.AddEntry(hRD, hRD.GetTitle(), 'lep')
leg.AddEntry(hMC, 'B^{0}#rightarrow J/#psi K*', 'le')
leg.AddEntry(gr_stat, 'B^{0}#rightarrow J/#psi K* (p_{T} weights)', 'lep')
leg.AddEntry(gr_sys, 'Weight systematics', 'f')

SetMaxToMaxHist([hRD, hMC])
c = rt.TCanvas('c', 'c', 50, 50, 800, 700)
c.SetTickx(0)
c.SetTicky(0)

pad = rt.TPad('pmain', 'pmain', 0, 0.25, 1, 1)
pad.SetBottomMargin(0.015)
pad.SetTopMargin(0.07)
pad.SetRightMargin(0.05)
pad.SetLeftMargin(0.15)
pad.Draw()
pad.cd()
hRD.Draw()
hMC.Draw('same')
gr_stat.Draw('p')
gr_sys.Draw('2')
leg.Draw()
catText.DrawLatexNDC(0.9, 0.6, 'Category: {}'.format(cat.name))
CMS_lumi.extraText = '     Internal'
CMS_lumi.integrated_lumi = lumi_tot
CMS_lumi.CMS_lumi(pad, -1, 33, cmsTextSize=0.75*1.2, lumiTextSize=0.6*1.2)

c.cd()
pad = rt.TPad('ppull', 'ppull', 0, 0, 1, 0.25)
pad.SetBottomMargin(0.5)
pad.SetTopMargin(0.03)
pad.SetRightMargin(0.05)
pad.SetLeftMargin(0.15)
pad.Draw('same')
pad.cd()
h_dr.GetYaxis().SetTitle('RD/MC')
t = 0.15
h_dr.GetYaxis().SetRangeUser(1 - 3*t, 1 + 3*t)
h_dr.GetYaxis().SetTitleOffset(0.5)
h_dr.GetYaxis().SetTitleSize(0.14)
h_dr.GetYaxis().SetLabelSize(0.15)
h_dr.GetYaxis().SetNdivisions(-203)
h_dr.GetXaxis().SetTitleOffset(0.95)
h_dr.GetXaxis().SetTitleSize(0.2)
h_dr.GetXaxis().SetLabelSize(0.18)
h_dr.GetXaxis().SetTickSize(0.07)
h_dr.Draw('E0')
h_mr.Draw('sameE0')
g_up.SetFillColor(rt.kRed-4)
g_up.SetFillStyle(3005)
g_up.Draw('F')
g_down.SetFillColor(rt.kRed-4)
g_down.SetFillStyle(3005)
g_down.Draw('F')
gh_dr = rt.TGraphErrors()
for i in range(1, h_dr.GetNbinsX()+1):
    gh_dr.SetPoint(i-1, h_dr.GetBinCenter(i), h_dr.GetBinContent(i))
    gh_dr.SetPointError(i-1, h_dr.GetBinError(i))
gh_dr.SetLineColor(h_dr.GetLineColor())
gh_dr.Draw('P0')
ax = h_dr.GetYaxis()
ax.ChangeLabel(1, -1, -1, -1, -1, -1, ' ')
ax.ChangeLabel(4, -1, -1, -1, -1, -1, ' ')

l = rt.TLine()
l.SetLineColor(rt.kGray+1)
l.SetLineWidth(1)
l.SetLineStyle(9)
x_low = h_dr.GetBinCenter(1)-0.5*h_dr.GetBinWidth(1)
x_high = h_dr.GetBinCenter(h_dr.GetNbinsX())+0.5*h_dr.GetBinWidth(h_dr.GetNbinsX())
l.DrawLine(x_low, 1, x_high, 1)

c.Draw()
c.SaveAs(webFolder+'closure_' + cat.name +'.png')


makePlot(mB_var, binning=[75, 5.24, 5.34], wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'],
         axis_title=['mass(#mu#mu#piK)', 'Normalized entries'], saveFig='mass_mumupiK_postcal_'+cat.name+'.png')
makePlot('mass_mumu', binning=[70, 3.0, 3.2], wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'],
         axis_title=['mass(#mu#mu)', 'Normalized entries'], saveFig='mass_mumu_postcal_'+cat.name+'.png')
makePlot('mass_candKst', binning=[70, 0.82, 0.97], wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'],
         axis_title=['mass(#piK)', 'Normalized entries'], saveFig='mass_Kst_postcal_'+cat.name+'.png')


b = {'Low': array('d', list(np.arange(7, 9.4, 0.1)) + [9.4]),
     'Mid': array('d', list(np.arange(9, 12.2, 0.1)) +[12.2]),
     'High': array('d', list(10+np.logspace(np.log10(12.2-10), np.log10(50), 30)))
    }
makePlot('trgMu_pt', binning=b[cat.name],  wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'],
         axis_title=['trigger muon p_{T} [GeV]', 'Normalized entries'], saveFig='trgMu_pt_postcal_'+cat.name+'.png')

for n in ['otherMu', 'K', 'pi', 'Jpsi', 'Kst']:
    makePlot(n+'_pt', binning=3*[None], axis_title=[n+' p_{T}', 'Normalized entries'],
              wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'], saveFig=n+'_pt_postcal_'+cat.name+'.png')


for n in ['trgMu', 'otherMu', 'K', 'pi', 'Jpsi', 'Kst']:
    makePlot(n+'_eta', binning=[50, -2, 2], axis_title=[n+' #eta', 'Normalized entries'],
              wMC=[dfMC['w'], dfMC['w']*dfMC['wEta']*dfMC['wPt']], legMC=['MC precal', 'MC postcal'], saveFig=n+'_eta_postcal_'+cat.name+'.png')


################ Additional tracks ###################

selMC = np.array(dfMC['N_goodAddTks'] == 1).astype(np.int)
selRD = np.array(dfRD['N_goodAddTks'] == 1).astype(np.int)
# makePlot('tkPt_0', binning=[50, 0.5, 7.5],  wMC=[selMC*dfMC['w'], selMC*dfMC['w']*dfMC['wEta']*dfMC['wPt']], wRD=selRD,
#          legMC=['MC precal', 'MC cal (B p_{T} & #eta)'],
#          axis_title=['Additonal track p_{T} [GeV]', 'Normalized entries'], saveFig='addTk_pt_precal_'+cat.name+'.png')

makePlot('massVis_wTk_0', binning=[50, 4.5, 8.5],  wMC=[selMC*dfMC['w'], selMC*dfMC['w']*dfMC['wEta']*dfMC['wPt']], wRD=selRD,
         legMC=['MC precal', 'MC cal (B p_{T} & #eta)'],
         axis_title=['Visible mass with add. track [GeV]', 'Normalized entries'], saveFig='addTk_mVis_precal_'+cat.name+'.png')

# Calibration

b=[40,0.5,6.5]
binWdith = (b[2] - b[1])/float(b[0])
hRD = create_TH1D(dfRD['tkPt_0'], name='hRD', title='data', weights=selRD,
                  axis_title=['Additional track p_{T} (reco) [GeV]',
                              '1/#sigma d#sigma/dp_{T} / '+'({:.2f})'.format(binWdith)],
                  binning=b, scale_histo='norm', opt='overflow+underflow'
                 )
hRD.SetMarkerStyle(15)

hMC = create_TH1D(dfMC['tkPt_0'], name='hMC', weights=selMC*dfMC['w']*dfMC['wEta']*dfMC['wPt'],
                  axis_title=['Additional track p_{T} (reco) [GeV]',
                              '1/#sigma d#sigma/dp_{T} / '+'({:.2f})'.format(binWdith)],
                  title = 'B#rightarrow J/#psi K*',
                  scale_histo='norm', color=0,
                  h2clone=hRD, opt='overflow+underflow')

getPolyCorrection(hRD, hMC, 1, 'addTk_pt_polyCoeff', args.verbose)
cal_addTK_pt = calibrationReader(calibration_file=dataLoc+'calibration/kinematicCalibration_Bd/addTk_pt_polyCoeff_{}_{}.pkl'.format(cat.name, version))
dfMC['wTkPt'] = cal_addTK_pt.getWeights(dfMC['tkPt_0'])
print 'Add tk calibration: q={:.3f} m={:.3f}'.format(*(cal_addTK_pt.beta))

hMC_corr = create_TH1D(dfMC['tkPt_0'], name='hMC', weights=selMC*dfMC['w']*dfMC['wEta']*dfMC['wPt']*dfMC['wTkPt'],
                       axis_title=['Additional track p_{T} (reco) [GeV]',
                                   '1/#sigma d#sigma/dp_{T} / '+'({:.2f})'.format(binWdith)],
                       title = 'B#rightarrow J/#psi K* (corr)',
                       scale_histo='norm', color=1,
                       h2clone=hRD, opt='overflow+underflow')

CMS_lumi.extraText = '      Internal'
c = make_ratio_plot([hRD, hMC, hMC_corr], ratio_bounds='auto', draw_opt='E1')

CMS_lumi.CMS_lumi(c, -1, 0)
c.pad1.SetTopMargin(0.07)
c.pad1.SetRightMargin(0.035)
c.pad2.SetRightMargin(0.035)
# c.pad2.SetLogy()
c.leg.SetY1(0.3)
c.leg.SetY2(0.5)
c.leg.SetX1(0.35)
c.leg.SetX2(0.7)
c.Draw()

c.pad1.cd()
catText.SetTextSize(0.04)
catText.DrawLatexNDC(0.9, 0.8, 'Category: {}'.format(cat.name))
c.SaveAs(webFolder+'addTk_pt_calibration_' + cat.name + '.png')
