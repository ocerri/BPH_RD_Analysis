#!/usr/bin/env python
#############################################################################
####                              Imports                                ####
#############################################################################
import sys, os, pickle, time, re
from glob import glob
sys.path.append('../lib')
sys.path.append('../analysis')
import itertools
import json
from multiprocessing import Pool
import commands

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.interpolate import interp1d
from array import array

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

from analysis_utilities import drawOnCMSCanvas, getEff
from histo_utilities import create_TH1D, create_TH2D, std_color_list, SetMaxToMaxHist
from cebefo_style import Set_2D_colz_graphics
from gridVarQ2Plot import col_dic, plot_gridVarQ2
from progressBar import ProgressBar
from categoriesDef import categories
from B02DstMu_selection import candidate_selection, trigger_selection, candidateSelection_stringList, candidateSelection_nameList

import argparse
parser = argparse.ArgumentParser()
#Example: python B2DstMu_skimCAND_v1.py --maxEntries 80000 --applyCorr
parser.add_argument ('--function', type=str, default='main', help='Function to perform')
parser.add_argument ('-d', '--dataset', type=str, default=[], help='Dataset(s) to run on', nargs='+')
parser.add_argument ('-p', '--parallelType', choices=['pool', 'jobs'], default='jobs', help='Function to perform')
parser.add_argument ('--maxEntries', type=int, default=1e15, help='Max number of events to be processed')
parser.add_argument ('--cat', type=str, default=['low', 'mid', 'high', 'none'], choices=['low', 'mid', 'high'], help='Category(ies)', nargs='+')
parser.add_argument ('--applyCorr', default=False, action='store_true', help='Switch to apply crrections')
######## Arguments not for user #####################
parser.add_argument ('--tmpDir', type=str, default=None, help='Temporary directory')
parser.add_argument ('--jN', type=int, default=None, help='Job number')
args = parser.parse_args()

#############################################################################
####                          Datset declaration                         ####
#############################################################################
MCloc = '../data/cmsMC_private/BPH_Tag-'
RDloc = '../data/cmsRD/ParkingBPH*/'

filesLocMap = {
'mu_PU0' : MCloc+'B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU0_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'mu_PU20': MCloc+'B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU20_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'mu_PU35': MCloc+'B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU35_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'mu_HQETPU0': MCloc+'B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_HQET2_central_PU0_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'tau_PU0': MCloc+'B0_TauNuDmst-pD0bar-kp-t2mnn_pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU0_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'tau'   : MCloc+'B0_TauNuDmst-pD0bar-kp-t2mnn_pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU20_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'Hc'    : MCloc+'B0_DmstHc-pD0bar-kp-Hc2mu_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_PU20_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
'Dstst' : MCloc+'Bp_MuNuDstst_DmstPi_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_ISGW2_PU20_10-2-3/ntuples_B2DstMu/out_CAND_*.root',
#
#
#
'dataB2DstMu': RDloc+'*_RDntuplizer_B2DstMu_200327_CAND.root'
# 'dataB2DstMu': RDloc+'Run2018D-05May2019promptD-v1_RDntuplizer_B2DstMu_200320_CAND.root',
# 'dataCombDmstMum': RDloc + 'Run2018D-05May2019promptD-v1_RDntuplizer_combDmstMum_200320_CAND.root'
}

def getTLVfromField(ev, n, idx, mass):
    v = rt.TLorentzVector()
    v.SetPtEtaPhiM(getattr(ev, n+'_pt')[idx],
                   getattr(ev, n+'_eta')[idx],
                   getattr(ev, n+'_phi')[idx],
                   mass)
    return v

class Container(object):
    pass

fBfield = rt.TFile.Open('../data/calibration/bFieldMap_2Dover3D.root', 'r')
hBfieldMapsRatio = fBfield.Get('bfieldMap')
def get_bFieldCorr3D(phi, eta, verbose=False):
    idx = hBfieldMapsRatio.GetBin(hBfieldMapsRatio.GetXaxis().FindBin(phi), hBfieldMapsRatio.GetYaxis().FindBin(eta))
    return 1./hBfieldMapsRatio.GetBinContent(idx)

def correctPt(pt, eta, phi, corr=None, smear=0):
    if corr is None:
        return pt
    elif corr=='RD':
        return pt * get_bFieldCorr3D(phi, eta)
    elif corr=='MC':
        return pt * (smear*np.random.randn() + 1.)

def compMass(pt1, pt2, eta1, eta2, phi1, phi2, m1, m2):
    E1 = np.hypot(m1, pt1*np.cosh(eta1))
    E2 = np.hypot(m2, pt2*np.cosh(eta2))
    p1p2 = pt1 * pt2 * (np.cos(phi1 - phi2) + np.sinh(eta1) * np.sinh(eta2))
    return np.sqrt(m1**2 + m2**2 + 2*(E1*E2 - p1p2))

def compMass3(pt1, pt2, pt3, eta1, eta2, eta3, phi1, phi2, phi3, m1, m2, m3):
    E1 = np.hypot(m1, pt1*np.cosh(eta1))
    E2 = np.hypot(m2, pt2*np.cosh(eta2))
    E3 = np.hypot(m3, pt3*np.cosh(eta3))
    p1p2 = pt1 * pt2 * (np.cos(phi1 - phi2) + np.sinh(eta1) * np.sinh(eta2))
    p1p3 = pt1 * pt3 * (np.cos(phi1 - phi3) + np.sinh(eta1) * np.sinh(eta3))
    p2p3 = pt2 * pt3 * (np.cos(phi2 - phi3) + np.sinh(eta2) * np.sinh(eta3))
    return np.sqrt(m1**2 + m2**2 + m3**2 + 2*(E1*E2 - p1p2) + 2*(E1*E3 - p1p3) + 2*(E2*E3 - p2p3))

def SumPt(pt1, pt2, phi1, phi2):
    pSq = pt1**2 + pt2**2 + 2*pt1*pt2*np.cos(phi1-phi2)
    return np.sqrt(pSq)

def insertOrdered(list, el):
    if len(list) == 0:
        return 0, [el]
    else:
        for il in range(len(list)):
            if list[il] < el:
                break
        # list = list[:il] + [el] + list[il:]
        return il, list[:il] + [el] + list[il:]

def extractEventInfos(j, ev, corr=None):
    m_mu   = 0.105658
    m_pi   = 0.139570
    m_K    = 0.493677
    m_D0   = 1.86483
    m_Dst  = 2.01026
    m_jpsi = 3.096916
    m_B0   = 5.27963

    e = Container()
    # print ''
    e.K_eta = ev.K_refitpiK_eta[j]
    e.K_phi = ev.K_refitpiK_phi[j]
    e.K_pt = correctPt(ev.K_refitpiK_pt[j], e.K_eta, e.K_phi, corr, 3e-3)

    e.pi_eta = ev.pi_refitpiK_eta[j]
    e.pi_phi = ev.pi_refitpiK_phi[j]
    e.pi_pt = correctPt(ev.pi_refitpiK_pt[j], e.pi_eta, e.pi_phi, corr, 3e-3)
    # print 'pi pt: {:1.2e}'.format((e.pi_pt - ev.pi_refitpiK_pt[j])/ev.pi_refitpiK_pt[j])

    e.mass_piK = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_pi, m_K)
    # print 'mass piK: {:1.4f} {:1.4f}'.format(ev.mass_piK[j], e.mass_piK)

    ptPost = SumPt(e.pi_pt, e.K_pt, e.pi_phi, e.K_phi)
    e.D0_pt = (ptPost/ev.D0_pt[j])*ev.D0_refitD0pismu_pt[j]
    # print 'D0 pT: {:.4f} {:.4f} {:.4f} {:.4f}'.format(ev.D0_pt[j], ptPost, ev.D0_refitD0pismu_pt[j], e.D0_pt)
    e.D0_eta = ev.D0_refitD0pismu_eta[j]
    e.D0_phi = ev.D0_refitD0pismu_phi[j]

    e.pis_eta = ev.pis_refitD0pismu_eta[j]
    e.pis_phi = ev.pis_refitD0pismu_phi[j]
    e.pis_pt = correctPt(ev.pis_refitD0pismu_pt[j], e.pis_eta, e.pis_phi, corr, 6e-3)
    # print 'pis pt: {:1.2e}'.format((e.pis_pt - ev.pis_refitD0pismu_pt[j])/ev.pis_refitD0pismu_pt[j])

    e.mass_D0pis = compMass(e.pis_pt, e.D0_pt, e.pis_eta, e.D0_eta, e.pis_phi, e.D0_phi, m_pi, e.mass_piK)
    # e.mass_D0pis_v2 = compMass3(e.pi_pt, e.K_pt, e.pis_pt, e.pi_eta, e.K_eta, e.pis_eta, e.pi_phi, e.K_phi, e.pis_phi, m_pi, m_K, m_pi)
    # print 'mass D0pis: {:.4f} {:.4f} {:.4f}'.format(ev.mass_D0pis[j], e.mass_D0pis, e.mass_D0pis_v2)

    e.mu_eta = ev.mu_refitD0pismu_eta[j]
    e.mu_phi = ev.mu_refitD0pismu_phi[j]
    e.mu_pt = correctPt(ev.mu_refitD0pismu_pt[j], e.mu_eta, e.mu_phi, corr, 3e-3)

    p4_D0 = rt.TLorentzVector()
    p4_D0.SetPtEtaPhiM(e.D0_pt, e.D0_eta, e.D0_phi, e.mass_piK)
    # p4_D0.SetPtEtaPhiM(e.D0_pt, e.D0_eta, e.D0_phi, m_D0)
    p4_pis = rt.TLorentzVector()
    p4_pis.SetPtEtaPhiM(e.pis_pt, e.pis_eta, e.pis_phi, m_pi)
    p4_mu = rt.TLorentzVector()
    p4_mu.SetPtEtaPhiM(e.mu_pt, e.mu_eta, e.mu_phi, m_mu)
    p4_Dst = p4_D0 + p4_pis
    e.Dst_pt = p4_Dst.Pt()
    e.Dst_eta = p4_Dst.Eta()
    e.Dst_phi = p4_Dst.Phi()
    p4_vis = p4_Dst + p4_mu
    e.mass_D0pismu = p4_vis.M()
    # print 'm_vis: {:.4f} {:.4f}'.format(ev.mass_D0pismu[j], e.m_vis)

    e.B_pt = p4_vis.Pt() * m_B0/ p4_vis.M()
    e.B_eta = ev.B_D0pismu_eta[j]
    e.B_phi = ev.B_D0pismu_phi[j]
    # print 'B0 pt: {:.4f} {:.4f}'.format(ev.B_D0pismu_pt[j], pt_B)
    p4_B = rt.TLorentzVector()
    p4_B.SetPtEtaPhiM(e.B_pt, e.B_eta, e.B_phi, m_B0);

    e.M2_miss = (p4_B - p4_vis).M2()
    # print 'M2_miss: {:.4f} {:.4f}'.format(ev.M2_miss_D0pismu[j], e.M2_miss)
    e.q2 = (p4_B - p4_Dst).M2()
    # print 'q2: {:.4f} {:.4f}'.format(ev.q2_D0pismu[j], e.q2)

    p4st_mu = rt.TLorentzVector(p4_mu)
    p4st_mu.Boost(-1*p4_B.BoostVector())
    e.Est_mu = p4st_mu.E()
    # print 'Est: {:.4f} {:.4f}'.format(ev.Est_mu_D0pismu[j], e.Est_mu)

    #----------------- Additional tracks -------------------#
    idx_st = 0
    for jjj in range(j):
        idx_st += int(ev.nTksAdd[jjj])
    idx_stop = int(idx_st + ev.nTksAdd[j])

    e.N_goodAddTks = 0
    e.tkPt = [] #np.zeros(idx_stop - idx_st)
    e.massVis_wTk = [] #np.zeros(idx_stop - idx_st)
    e.massHad_wTk = [] #np.zeros(idx_stop - idx_st)
    e.massMuTk = [] #np.zeros(idx_stop - idx_st)

    for jj in range(idx_st, idx_stop):
        eta = ev.tksAdd_eta[jj]
        phi = ev.tksAdd_phi[jj]
        pt = correctPt(ev.tksAdd_pt[jj], eta, phi, corr, 6e-3)
        p4_tk = rt.TLorentzVector()
        p4_tk.SetPtEtaPhiM(pt, eta, phi, m_pi)

        mVis_wTk = (p4_vis + p4_tk).M()
        # print 'm_vis_wTk: {:.4f} {:.4f}'.format(ev.tksAdd_massVis[jj], mVis_wTk)

        if mVis_wTk < m_B0 and ev.tksAdd_cos_PV[jj]>0.95:
            e.N_goodAddTks += 1
            idx, e.tkPt = insertOrdered(e.tkPt, pt)
            e.massVis_wTk.insert(idx, mVis_wTk)
            e.massHad_wTk.insert(idx, (p4_Dst + p4_tk).M())
            e.massMuTk.insert(idx, (p4_mu + p4_tk).M())

    if e.N_goodAddTks < 2:
        for l in [e.tkPt, e.massVis_wTk, e.massHad_wTk, e.massMuTk]:
            l += [0, 0]

    return e

def makeSelection(inputs):
    n, tag, filepath, leafs_names, cat, idxInt, corr, skipCut, serial = inputs
    N_accepted_cand = []
    N_accepted_tot = 0

    tree = rt.TChain('outA/Tevts')
    lastIdxDisc = -1
    for fn in glob(filepath):
        tree.Add(fn)
        if tree.GetEntries() + lastIdxDisc < idxInt[0]:
            lastIdxDisc += tree.GetEntries()
            tree = rt.TChain('outA/Tevts')
        elif tree.GetEntries() + lastIdxDisc > idxInt[1]:
            break

    nDiscEvts = lastIdxDisc + 1

    if serial:
        pb = ProgressBar(maxEntry=idxInt[1]+1)
    else:
        perc = int((idxInt[1]-idxInt[0])*0.3)

    output = np.zeros((int(1.5*(idxInt[1]-idxInt[0]+1)), len(leafs_names)))

    for i_ev, ev in enumerate(tree):
        i_ev += nDiscEvts
        if i_ev < idxInt[0]:
            continue
        if i_ev > idxInt[1]:
            break

        if serial:
            pb.show(i_ev-idxInt[0])
        elif (i_ev-idxInt[0]) % perc == 0:
            print tag, ': {:.0f}%'.format(100*(i_ev+1-idxInt[0])/(idxInt[1]-idxInt[0]))
        N_acc = 0

        ev_output = []
        for j in range(ev.pval_piK.size()):
            idxTrg = int(ev.mu_trgMu_idx[j])
            if not cat is None:
                if not trigger_selection(idxTrg, ev, cat):
                    continue

            evEx = extractEventInfos(j, ev, corr)
            if not skipCut == 'all':
                if not candidate_selection(j, ev, evEx, skipCut):
                    continue

            N_acc += 1

            aux = (evEx.q2, evEx.Est_mu, evEx.M2_miss,
                   evEx.mu_pt, evEx.mu_eta, evEx.mu_phi, ev.trgMu_sigdxy[idxTrg],
                   evEx.B_pt, evEx.B_eta, evEx.B_phi,
                   evEx.Dst_pt, evEx.Dst_eta, evEx.Dst_phi,
                   evEx.D0_pt, evEx.D0_eta, evEx.D0_phi,
                   evEx.pi_pt, evEx.pi_eta, evEx.pi_phi,
                   evEx.K_pt, evEx.K_eta, evEx.K_phi,
                   ev.pval_piK[j], ev.sigdxy_vtxD0_PV[j],
                   evEx.pis_pt, evEx.pis_eta, evEx.pis_phi,
                   ev.pval_D0pis[j],
                   evEx.mass_piK, evEx.mass_D0pis, evEx.mass_D0pismu,
                   ev.pval_D0pismu[j], ev.cos_D0pismu_PV[j], ev.cosT_D0pismu_PV[j],
                   evEx.N_goodAddTks,
                   evEx.tkPt[0], evEx.tkPt[1],
                   evEx.massVis_wTk[0], evEx.massVis_wTk[1],
                   evEx.massHad_wTk[0], evEx.massHad_wTk[1],
                   evEx.massMuTk[0], evEx.massMuTk[1],
                   trigger_selection(idxTrg, ev, categories['low']),
                   trigger_selection(idxTrg, ev, categories['mid']),
                   trigger_selection(idxTrg, ev, categories['high']),
                   ev.N_vertexes
                  )
            if not 'data' in n:
                aux += (ev.MC_q2, ev.MC_Est_mu, ev.MC_M2_miss,
                        ev.MC_B_pt, ev.MC_B_eta, ev.MC_B_phi,
                        ev.MC_Dst_pt, ev.MC_Dst_eta, ev.MC_Dst_phi,
                        ev.MC_mu_pt, ev.MC_mu_eta, ev.MC_mu_phi, ev.MC_mu_IP,
                        ev.MC_pi_pt, ev.MC_pi_eta, ev.MC_pi_phi,
                        ev.MC_K_pt, ev.MC_K_eta, ev.MC_K_phi,
                        ev.MC_pis_pt, ev.MC_pis_eta, ev.MC_pis_phi,
                        ev.MC_idxCand == j
                       )
            if 'mu' in n or 'tau' in n:
                aux += (ev.wh_CLNCentral,
                        ev.wh_CLNR0Down, ev.wh_CLNR0Up,
                        ev.wh_CLNR1Down, ev.wh_CLNR1Up,
                        ev.wh_CLNR2Down, ev.wh_CLNR2Up,
                        ev.wh_CLNRhoSqDown, ev.wh_CLNRhoSqUp,
                       )
            ev_output.append(aux)

        N_acc = len(ev_output)
        idx = 0
        if N_acc > 1:
            if 'data' in n:
                idx = np.random.randint(len(ev_output))
            else:
                #Get matched can preferably
                varIdx = leafs_names.index('MC_idxMatch')
                goodIdx = np.nonzero([o[varIdx] for o in ev_output])[0]
                if goodIdx.shape[0] > 0:
                    auxIdx = np.random.randint(goodIdx.shape[0])
                    idx = goodIdx[auxIdx]
                else:
                    idx = np.random.randint(len(ev_output))
        if N_acc > 0:
            output[N_accepted_tot] = ev_output[idx]
            N_accepted_tot += 1
            N_accepted_cand.append(N_acc)

    output = output[:N_accepted_tot]
    if not serial:
        print tag, ': done'
    return [output, N_accepted_cand]

def create_dSet(n, filepath, cat, applyCorrections=False, skipCut=[], maxEntries=args.maxEntries):
    if cat is None:
        catName = 'None'
    else:
        catName = cat.name
    print n, catName
    if 'data' in n:
        loc = '../data/cmsRD/skimmed/'+ n.replace('data', '')
        out = re.search('20[01][1-9][0-3][0-9]', filepath)
        if out is None:
            print filepath
            raise
        fskimmed_name = loc + '_' + out.group(0) + '_' + catName
        N_evts_per_job = 300000
    else:
        d = os.path.dirname(filepath) + '/skimmed/'
        if not os.path.isdir(d):
            os.makedirs(d)
        fskimmed_name = d + catName
        N_evts_per_job = 30000
    if not skipCut == []:
        print 'Skipping cut(s)', skipCut
        if skipCut == 'all':
            fskimmed_name += '_skipall'
        else:
            fskimmed_name += '_skip'+'-'.join([str(i) for i in skipCut])
    if applyCorrections:
        print 'Appling corrections'
        fskimmed_name += '_corr'
    else:
        fskimmed_name += '_bare'
    fskimmed_name += '.root'
    logfile = fskimmed_name.replace('.root', '.log')
    if os.path.isfile(fskimmed_name) and not n in recreate:
        print 'Already present'
    else:
        tree = rt.TChain('outA/Tevts')
        for fn in glob(filepath):
            tree.Add(fn)
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEntries, tree.GetEntries())
        print n, ': Total number of candidate events =', N_cand_in

        leafs_names = ['q2', 'Est_mu', 'M2_miss',
                       'mu_pt', 'mu_eta', 'mu_phi', 'mu_sigdxy',
                       'B_pt', 'B_eta', 'B_phi',
                       'Dst_pt', 'Dst_eta', 'Dst_phi',
                       'D0_pt', 'D0_eta', 'D0_phi',
                       'pi_pt', 'pi_eta', 'pi_phi',
                       'K_pt', 'K_eta', 'K_phi',
                       'pval_piK', 'sigdxy_vtxD0_PV',
                       'pis_pt', 'pis_eta', 'pis_phi',
                       'pval_D0pis',
                       'mass_piK', 'mass_D0pis', 'mass_D0pismu',
                       'pval_D0pismu', 'cos_D0pismu_PV', 'cosT_D0pismu_PV',
                       'N_goodAddTks',
                       'tkPt_0', 'tkPt_1',
                       'tkMassVis_0', 'tkMassVis_1',
                       'tkMassHad_0', 'tkMassHad_1',
                       'tkMassMuTk_0', 'tkMassMuTk_1',
                       'cat_low', 'cat_mid', 'cat_high',
                       'N_vtx'
                      ]
        if not 'data' in n:
            leafs_names += ['MC_q2', 'MC_Est_mu', 'MC_M2_miss',
                            'MC_B_pt', 'MC_B_eta', 'MC_B_phi',
                            'MC_Dst_pt', 'MC_Dst_eta', 'MC_Dst_phi',
                            'MC_mu_pt', 'MC_mu_eta', 'MC_mu_phi', 'MC_mu_IP',
                            'MC_pi_pt', 'MC_pi_eta', 'MC_pi_phi',
                            'MC_K_pt', 'MC_K_eta', 'MC_K_phi',
                            'MC_pis_pt', 'MC_pis_eta', 'MC_pis_phi',
                            'MC_idxMatch'
                           ]
        if 'mu' in n or 'tau' in n:
            leafs_names += ['wh_CLNCentral',
                            'wh_CLNR0Down', 'wh_CLNR0Up',
                            'wh_CLNR1Down', 'wh_CLNR1Up',
                            'wh_CLNR2Down', 'wh_CLNR2Up',
                            'wh_CLNRhoSqDown', 'wh_CLNRhoSqUp']

        applyCorr = None
        if applyCorrections:
            applyCorr = 'MC'
            if 'data' in n:
                applyCorr = 'RD'

        if N_cand_in < 1.5*N_evts_per_job:
            output, N_accepted_cand = makeSelection([n, '', filepath, leafs_names, cat,
                                                     [0, N_cand_in-1], applyCorr, skipCut, True])
        else:
            pdiv = list(np.arange(0, N_cand_in, N_evts_per_job))
            if not pdiv[-1] == N_cand_in:
                pdiv.append(N_cand_in)
            print 'Will be divided into ' + str(len(pdiv)-1) + ' jobs'
            inputs = []
            for i in range(1, len(pdiv)):
                corr = 0
                if i == 1:
                    corr = -1
                inputs.append([n, str(i), filepath, leafs_names, cat, [pdiv[i-1]+1+corr, pdiv[i]], applyCorr, skipCut, False])
            print ' '

            start = time.time()
            if args.parallelType == 'pool':
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection, inputs)
            elif args.parallelType == 'jobs':
                tmpDir = 'tmp/B2DstMu_skimCAND_' + n
                os.system('rm -rf ' + tmpDir + '/out')
                os.system('rm -rf ' + tmpDir + '/*.p')
                os.makedirs(tmpDir + '/out')
                for ii, inAux in enumerate(inputs):
                    pickle.dump( inAux, open( tmpDir+'/input_{}.p'.format(ii), 'wb' ) )
                createSubmissionFile(tmpDir, len(inputs))
                print 'Submitting jobs'
                cmd = 'condor_submit {}/jobs.jdl'.format(tmpDir)
                cmd += ' -batch-name skim_' + n
                status, output = commands.getstatusoutput(cmd)
                if status !=0:
                    print 'Error in processing command:\n   ['+cmd+']'
                    print 'Output:\n   ['+output+'] \n'
                print 'Job submitted'
                print 'Waiting for jobs to be finished'
                time.sleep(20)
                proceed=False
                while not proceed:
                    status, output = commands.getstatusoutput('condor_q')
                    found = False
                    for line in output.split('\n'):
                        if 'skim_'+n in line:
                            print line
                            time.sleep(10)
                            found = True
                    proceed = not found

                outputs = []
                for ii in range(len(inputs)):
                    o = pickle.load( open( tmpDir+'/output_{}.p'.format(ii), 'rb' ) )
                    outputs.append(o)

            output = np.concatenate(tuple([o[0] for o in outputs]))
            N_accepted_cand = []
            for o in outputs: N_accepted_cand += o[1]
            print 'Total time: {:.1f} min'.format((time.time()-start)/60.)


        dset = pd.DataFrame(output, columns=leafs_names)
        if not os.path.isdir(os.path.dirname(fskimmed_name)):
            os.makedirs(os.path.dirname(fskimmed_name))
        rtnp.array2root(dset.to_records(), fskimmed_name, treename='Tevts', mode='RECREATE')

        with open(logfile, 'w') as f:
            ln = 'Number of candidates per events\n{'
            try:
                ln += ', '.join(['{}:{}'.format(i, N_accepted_cand.count(i)) for i in range(1, np.max(N_accepted_cand)+1)])
            except:
                ln += 'N/A'
            ln += '}\n'
            f.write(ln)
            f.write('N_analyzed: '+str(N_cand_in)+'\n')
            f.write('N_accepted: '+str(dset.shape[0])+'\n')
            e = getEff(dset.shape[0], N_cand_in)
            f.write('Eff: {:.3f} +/- {:.3f} %'.format(1e2*e[0], 1e2*e[1])+'\n')

    os.system('echo '+logfile+';cat '+logfile + ';echo ')

def createSubmissionFile(tmpDir, njobs):
    fjob = open(tmpDir+'/job.sh', 'w')
    fjob.write('#!/bin/bash\n')
    fjob.write('source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /storage/user/ocerri/CMSSW_10_2_3/; eval `scramv1 runtime -sh`\n')
    fjob.write('cd /storage/user/ocerri/BPhysics/scripts\n')
    fjob.write('python B2DstMu_skimCAND_v1.py --function makeSel --tmpDir $1 --jN $2\n')
    os.system('chmod +x {}/job.sh'.format(tmpDir))

    fsub = open(tmpDir+'/jobs.jdl', 'w')
    fsub.write('executable    = ' + tmpDir+'/job.sh')
    fsub.write('\n')
    fsub.write('arguments     = {} $(ProcId) '.format(tmpDir))
    fsub.write('\n')
    fsub.write('output        = {}/out/job_$(ProcId)_$(ClusterId).out'.format(tmpDir))
    fsub.write('\n')
    fsub.write('error         = {}/out/job_$(ProcId)_$(ClusterId).err'.format(tmpDir))
    fsub.write('\n')
    fsub.write('log           = {}/out/job_$(ProcId)_$(ClusterId).log'.format(tmpDir))
    fsub.write('\n')
    fsub.write('WHEN_TO_TRANSFER_OUTPUT = ON_EXIT_OR_EVICT')
    fsub.write('\n')
    fsub.write('+MaxRuntime   = 3600')
    fsub.write('\n')
    fsub.write('+RunAsOwner = True')
    fsub.write('\n')
    fsub.write('+InteractiveUser = True')
    fsub.write('\n')
    fsub.write('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7"')
    fsub.write('\n')
    fsub.write('+SingularityBindCVMFS = True')
    fsub.write('\n')
    fsub.write('run_as_owner = True')
    fsub.write('\n')
    fsub.write('RequestDisk = 2000000')
    fsub.write('\n')
    fsub.write('RequestMemory = 2500')
    fsub.write('\n')
    fsub.write('RequestCpus = 1')
    fsub.write('\n')
    fsub.write('x509userproxy = $ENV(X509_USER_PROXY)')
    fsub.write('\n')
    fsub.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)')
    fsub.write('\n')
    fsub.write('on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)')   # Send the job to Held state on failure.
    fsub.write('\n')
    fsub.write('periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*20))')   # Periodically retry the jobs for 3 times with an interval of 20 minutes.
    fsub.write('\n')
    fsub.write('+PeriodicRemove = ((JobStatus =?= 2) && ((MemoryUsage =!= UNDEFINED && MemoryUsage > 2.5*RequestMemory)))')
    fsub.write('\n')
    fsub.write('max_retries    = 3')
    fsub.write('\n')
    fsub.write('requirements   = Machine =!= LastRemoteHost')
    fsub.write('\n')
    fsub.write('universe = vanilla')
    fsub.write('\n')
    fsub.write('queue '+str(njobs))
    fsub.write('\n')
    fsub.close()

if __name__ == "__main__":
    if args.function == 'main':
        file_loc = {}
        for n in args.dataset:
            file_loc[n] = filesLocMap[n]
        if len(file_loc.keys()) == 0:
            print 'No dataset provided'
            exit()

        recreate = file_loc.keys()
        # recreate = []

        if 'none' in args.cat:
            print 'Running w/o category (ignoring other categories)'
            for n, fp in file_loc.iteritems():
                    create_dSet(n, fp, cat=None, applyCorrections=args.applyCorr)

            # for n, fp in file_loc.iteritems():
            #     create_dSet(n, fp, None, skipCut='all', applyCorrections=args.applyCorr)

        else:
            # for cn in args.cat:
            #     for n, fp in file_loc.iteritems():
            #         create_dSet(n, fp, categories[cn], applyCorrections=args.applyCorr)

            skip = []
            # skip.append([6, 11, 12]) #Mass D0, D* and D*-D0 (m piK)
            # skip.append([16]) #Visible mass (m D0pismu)
            skip.append([17]) #Additional tracks
            for idx in skip:
                for cn in args.cat:
                    for n, fp in file_loc.iteritems():
                        create_dSet(n, fp, categories[cn], skipCut=idx, applyCorrections=args.applyCorr)
    elif args.function == 'makeSel':
        tmpDir = args.tmpDir
        input = pickle.load( open( tmpDir+'/input_{}.p'.format(args.jN), 'rb' ) )
        output = makeSelection(input)
        pickle.dump(output, open( tmpDir+'/output_{}.p'.format(args.jN), 'wb' ) )

    else:
        print args.function, 'not recognized'
