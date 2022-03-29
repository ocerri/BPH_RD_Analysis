#!/usr/bin/env python
"""
Script for running the skimmer. After you have created the ntuples, you can run
this script by running:

    $ python B2DstMu_skimCAND_v1.py -d "Bd_MuNuDst\Z" --cat low

The argument to `-d` can be a regular expression, so if you want it to match
all the soft QCD MC events you can run:

    $ python B2DstMu_skimCAND_v1.py -d SoftQCDnonD

You can also pass a list of regular expressions:

    $ python B2DstMu_skimCAND_v1.py -d BdToDstarMuNu,BdToDstarTauNu

"""
import sys, os, pickle, time, re, yaml
from glob import glob
from multiprocessing import Pool
import commands
from os.path import join
import numpy as np
import pandas as pd

# try:
# except ImportError:
#     print >> sys.stderr, "Failed to import analysis_utilities."
#     print >> sys.stderr, "Did you remember to source the env.sh file in the repo?"
#     sys.exit(1)
sys.path.append('../lib')
sys.path.append('../analysis')
from analysis_utilities import getEff
from progressBar import ProgressBar
from categoriesDef import categories
from B2DstMu_selection import candidate_selection, trigger_selection

import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument ('--function', type=str, default='main', help='Function to perform')
parser.add_argument ('-d', '--dataset', type=str, default=[], help='Dataset(s) to run on or regular expression for them', nargs='+')
parser.add_argument ('--skimTag', type=str, default='', help='Tag to append at the name of the skimmed files directory')
parser.add_argument ('-p', '--parallelType', choices=['pool', 'jobs', 'serial'], default='jobs', help='Function to perform')
parser.add_argument ('--maxEvents', type=int, default=1e15, help='Max number of events to be processed')
parser.add_argument ('--recreate', default=False, action='store_true', help='Recreate even if file already present')
parser.add_argument ('--applyCorr', default=False, action='store_true', help='Switch to apply crrections')
parser.add_argument ('--region', type=str, default='all', choices=['signal', 'trkControl', 'all'], help='Region to skim: signal (0 tracks) or track control (1+)')
parser.add_argument ('--cat', type=str, default=['high', 'mid', 'low'], choices=['single', 'low', 'mid', 'high', 'none'], help='Category(ies)', nargs='+')
parser.add_argument ('--skipCut', type=str, default='', choices=['all', '11', '13', '14', '16', '17'], help='Cut to skip.\nAll: skip all the cuts\n16:Visible mass cut\n17: additional tracks cut')
######## Arguments not for user #####################
parser.add_argument ('--tmpDir', type=str, default=None, help='Temporary directory')
parser.add_argument ('--jN', type=int, default=None, help='Job number')
args = parser.parse_args()

# python B2DstMu_skimCAND_v1.py -d "Bd_MuNuDst\Z"
# python B2DstMu_skimCAND_v1.py -d "B[us]"
# python B2DstMu_skimCAND_v1.py -d "Bd_Dst" "Bd_DDs1" "Bd_Tau" "Bd_MuNuDstPi" "B_DstDXX"

#############################################################################
####                          Datset declaration                         ####
#############################################################################
filesLocMap = {}

root = '/storage/af/group/rdst_analysis/BPhysics/data'
MCloc = join(root,'cmsMC/')
MCend = 'ntuples_B2DstMu_220326/out_CAND_*.root'
MC_samples = ['Bd_MuNuDst',
              'Bd_TauNuDst',
              'Bu_MuNuDstPi',       'Bd_MuNuDstPi',
              'Bd_MuNuDstPiPi_v3',  'Bu_MuNuDstPiPi_v3',
              'Bu_TauNuDstPi',      'Bd_TauNuDstPi',
              'Bd_TauNuDstPiPi',    'Bu_TauNuDstPiPi',
              'Bs_MuNuDstK',        'Bs_TauNuDstK',
              'Bd_DstDu',           'Bd_DstDd',          'Bd_DstDs',
              'Bu_DstDu',           'Bu_DstDd',          'Bs_DstDs',
              'Bd_DDs1',            'Bu_DDs1',           'B_DstDXX',
              # Others
              # 'DstKu_KuToMu',
              # 'Mu_Enriched'
              ]

sampleFile = '/storage/af/user/ocerri/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/production/samples.yml'
samples = yaml.load(open(sampleFile))['samples']
for s in MC_samples:
    # if s.endswith('_v3'):
    #     filesLocMap[s] = join(MCloc, samples[s]['dataset'], MCend.replace('220201', '220201_ffD2S'))
    # else:
    filesLocMap[s] = join(MCloc, samples[s]['dataset'], MCend)

RDloc = join(root,'cmsRD/ParkingBPH*/')
filesLocMap['data'] = join(RDloc, '*_B2DstMu_220311_CAND.root')
filesLocMap['data_SS'] = join(RDloc, '*_SSDstMu_220311_CAND.root')

def getTLVfromField(ev, n, idx, mass):
    v = rt.TLorentzVector()
    v.SetPtEtaPhiM(getattr(ev, n+'_pt')[idx],
                   getattr(ev, n+'_eta')[idx],
                   getattr(ev, n+'_phi')[idx],
                   mass)
    return v

class Container(object):
    pass

# hBfieldMapsRatio = None
# if hBfieldMapsRatio is None:
fBfield = rt.TFile.Open('/storage/af/group/rdst_analysis/BPhysics/data/calibration/bFieldMap_2Dover3D.root', 'r')
hBfieldMapsRatio = fBfield.Get('bfieldMap')

def get_bFieldCorr3D(phi, eta, verbose=False):
    global hBfieldMapsRatio
    if np.abs(eta) > 2.4:
        eta = 2.39*np.sign(eta)
    if np.abs(phi) > np.pi:
        phi = phi - 2*np.pi*np.sign(phi)
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
        else:
            il += 1
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
    # print '------> <-----'
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
    e.D0pismu_eta = p4_vis.Eta()
    e.D0pismu_phi = p4_vis.Phi()
    e.mass_D0pismu = p4_vis.M()

    p4_mu_as_pi = rt.TLorentzVector()
    p4_mu_as_pi.SetPtEtaPhiM(e.mu_pt, e.mu_eta, e.mu_phi, m_pi)
    e.mass_D0pismu_muASpi = (p4_Dst + p4_mu_as_pi).M()

    p4_mu_as_K = rt.TLorentzVector()
    p4_mu_as_K.SetPtEtaPhiM(e.mu_pt, e.mu_eta, e.mu_phi, m_K)
    e.mass_D0pismu_muASK = (p4_Dst + p4_mu_as_K).M()
    # print 'm_vis: {:.4f} {:.4f}'.format(ev.mass_D0pismu[j], e.m_vis)

    e.B_pt = p4_vis.Pt() * m_B0/ p4_vis.M()

    # Using direction from vertex
    e.B_eta = ev.B_D0pismu_eta[j]
    e.B_phi = ev.B_D0pismu_phi[j]

    p4_B = rt.TLorentzVector()
    p4_B.SetPtEtaPhiM(e.B_pt, e.B_eta, e.B_phi, m_B0)

    e.M2_miss = (p4_B - p4_vis).M2()
    e.U_miss = (p4_B - p4_vis).E() - (p4_B - p4_vis).P()
    e.q2 = (p4_B - p4_Dst).M2()

    p4st_mu = rt.TLorentzVector(p4_mu)
    p4st_mu.Boost(-1*p4_B.BoostVector())
    e.Est_mu = p4st_mu.E()

    # Using collinearity approximation
    p4_B_coll = rt.TLorentzVector()
    p4_B_coll.SetPtEtaPhiM(e.B_pt, e.D0pismu_eta, e.D0pismu_phi, m_B0)

    e.M2_miss_coll = (p4_B_coll - p4_vis).M2()
    e.q2_coll = (p4_B_coll - p4_Dst).M2()

    p4st_mu = rt.TLorentzVector(p4_mu)
    p4st_mu.Boost(-1*p4_B_coll.BoostVector())
    e.Est_mu_coll = p4st_mu.E()

    #----------------- Additional tracks -------------------#
    idx_st = 0
    for jjj in range(j):
        idx_st += int(ev.nTksAdd[jjj])
    idx_stop = int(idx_st + ev.nTksAdd[j])

    e.N_goodAddTks = 0
    e.tkCharge = []
    e.tkPdgId = []
    e.tkPt = []
    e.tkPtError = []
    e.tkEta = []
    e.tkPhi = []
    e.tk_lostInnerHits = []
    e.tk_pval = []
    e.tk_sigIP3D_vtxB = []
    e.MC_tkFlag = []
    e.MC_tkFromMainB = []
    e.MC_tkPdgId = []
    e.MC_tk_dphi = []
    e.MC_tk_deta = []
    e.MC_tk_dpt = []
    e.MC_tkMotherPdgId = []
    e.MC_tkMotherMotherPdgId = []
    e.massVis_wTk = []
    e.massHad_wTk = []
    e.massMuTk = []
    e.massDTk = []
    e.mass2MissTk = []
    e.UmissTk = []

    p4_sumGoodTks = rt.TLorentzVector()
    p4_tks = []
    for jj in range(idx_st, idx_stop):
        pval = ev.tksAdd_pval[jj]
        if pval < 0.1:
            continue

        if ev.tksAdd_lostInnerHits[jj] > 0:
            continue

        eta = ev.tksAdd_eta[jj]
        if np.abs(eta) >= 2.4:
            continue
        phi = ev.tksAdd_phi[jj]
        pt = correctPt(ev.tksAdd_pt[jj], eta, phi, corr, 2e-3)
        # if pt < 1.0:
        if pt < 0.55:
            continue
        #Avoid tracks duplicates
        duplicate = False
        for n in ['mu', 'pi', 'K', 'pis']:
            dphi = phi - getattr(e, n+'_phi')
            if dphi > np.pi: dphi -= 2*np.pi
            if dphi < -np.pi: dphi += 2*np.pi
            dR = np.hypot(dphi, eta - getattr(e, n+'_eta'))
            dPt = np.abs(getattr(e, n+'_pt') - pt)/getattr(e, n+'_pt')
            if dPt < 0.03 and dR < 0.001:
                duplicate=True
                break
        if duplicate:
            continue
        #
        p4_tk = rt.TLorentzVector()
        p4_tk.SetPtEtaPhiM(pt, eta, phi, m_pi)

        mVis_wTk = (p4_vis + p4_tk).M()
        # print 'm_vis_wTk: {:.4f} {:.4f}'.format(ev.tksAdd_massVis[jj], mVis_wTk)

        if mVis_wTk < m_B0 and ev.tksAdd_cos_PV[jj]>0.95:
            e.N_goodAddTks += 1
            idx, e.tkPt = insertOrdered(e.tkPt, pt)
            e.tkPtError.insert(idx, ev.tksAdd_ptError[jj])
            e.tkEta.insert(idx, eta)
            e.tkPhi.insert(idx, phi)
            e.tk_lostInnerHits.insert(idx, ev.tksAdd_lostInnerHits[jj])
            e.tk_pval.insert(idx, ev.tksAdd_pval[jj])
            e.tk_sigIP3D_vtxB.insert(idx, ev.tksAdd_sigIP3D_vtxB[jj])
            e.tkCharge.insert(idx, ev.tksAdd_charge[jj]*ev.mu_charge[j])
            e.tkPdgId.insert(idx, ev.tksAdd_pdgId[jj])
            e.massVis_wTk.insert(idx, mVis_wTk)
            e.massHad_wTk.insert(idx, (p4_Dst + p4_tk).M())
            e.massMuTk.insert(idx, (p4_mu + p4_tk).M())
            e.massDTk.insert(idx, (p4_D0 + p4_tk).M())
            p_miss = p4_B - p4_vis - p4_tk
            e.mass2MissTk.insert(idx, p_miss.M2())
            e.UmissTk.insert(idx, p_miss.E() - p_miss.P())
            if hasattr(ev, 'MC_addTkFlag'):
                e.MC_tkFlag.insert(idx, ev.MC_addTkFlag[jj])
                e.MC_tkFromMainB.insert(idx, ev.MC_addTk_fromMainB[jj])
                e.MC_tk_dphi.insert(idx, ev.MC_addTk_dPhi[jj])
                e.MC_tk_deta.insert(idx, ev.MC_addTk_dEta[jj])
                e.MC_tk_dpt.insert(idx, ev.MC_addTk_dPt[jj])
                e.MC_tkPdgId.insert(idx, ev.MC_addTk_pdgId[jj])
                e.MC_tkMotherPdgId.insert(idx, ev.MC_addTk_pdgIdMother[jj])
                e.MC_tkMotherMotherPdgId.insert(idx, ev.MC_addTk_pdgIdMotherMother[jj])

            p4_tks.insert(idx,p4_tk)
            p4_sumGoodTks += p4_tk

    while len(p4_tks) < 3:
        p4_tks.append(rt.TLorentzVector())

    p4_vis_wTks = p4_vis + p4_sumGoodTks
    e.massVisTks1 = (p4_vis + p4_tks[0]).M()
    e.massVisTks2 = (p4_vis + p4_tks[0] + p4_tks[1]).M()
    e.massVisTks3 = (p4_vis + p4_tks[0] + p4_tks[1] + p4_tks[2]).M()
    e.massVisTks = p4_vis_wTks.M()
    e.massHadTks = (p4_Dst + p4_sumGoodTks).M()
    e.massHadTks1 = (p4_Dst + p4_tks[0]).M()
    e.massHadTks2 = (p4_Dst + p4_tks[0] + p4_tks[1]).M()
    e.massHadTks3 = (p4_Dst + p4_tks[0] + p4_tks[1] + p4_tks[2]).M()
    e.massHadTks_DstMassConstraint = (p4_Dst + p4_sumGoodTks).M() - p4_Dst.M() + m_Dst

    if e.N_goodAddTks == 2:
        e.massTks_pipi = p4_sumGoodTks.M()

        tk0_pi = rt.TLorentzVector()
        tk0_pi.SetPtEtaPhiM(e.tkPt[0], e.tkEta[0], e.tkPhi[0], m_pi)

        tk0_K = rt.TLorentzVector()
        tk0_K.SetPtEtaPhiM(e.tkPt[0], e.tkEta[0], e.tkPhi[0], m_K)

        tk1_pi = rt.TLorentzVector()
        tk1_pi.SetPtEtaPhiM(e.tkPt[1], e.tkEta[1], e.tkPhi[1], m_pi)

        tk1_K = rt.TLorentzVector()
        tk1_K.SetPtEtaPhiM(e.tkPt[1], e.tkEta[1], e.tkPhi[1], m_K)

        e.massTks_KK = (tk0_K + tk1_K).M()

        if e.tkCharge[0] > 0 and e.tkCharge[1] < 0:
            e.massTks_piK = (tk0_pi + tk1_K).M()
            e.massTks_Kpi = (tk0_K + tk1_pi).M()
        elif e.tkCharge[0] < 0 and e.tkCharge[1] > 0:
            e.massTks_Kpi = (tk0_pi + tk1_K).M()
            e.massTks_piK = (tk0_K + tk1_pi).M()
        else:
            e.massTks_piK = (tk0_pi + tk1_K).M()
            e.massTks_Kpi = (tk0_K + tk1_pi).M()
    else:
        e.massTks_pipi = 0
        e.massTks_KK = 0
        e.massTks_piK = 0
        e.massTks_Kpi = 0



    e.BwTks_pt = p4_vis_wTks.Pt() * m_B0/ p4_vis_wTks.M()
    p4_BwTks = rt.TLorentzVector()
    p4_BwTks.SetPtEtaPhiM(e.BwTks_pt, e.B_eta, e.B_phi, m_B0);
    p_miss_wTks = p4_BwTks - p4_vis_wTks
    e.EmissTks = p_miss_wTks.E()
    e.PmissTks = p_miss_wTks.P()
    e.UmissTks = p_miss_wTks.E() - p_miss_wTks.P()

    e.M2missTks = (p4_BwTks - p4_vis_wTks).M2()
    e.q2Tks = (p4_BwTks - p4_Dst - p4_sumGoodTks).M2()


    if e.N_goodAddTks < 3:
        auxList = [e.tkCharge, e.tkPdgId, e.tkPt, e.tkPtError, e.tkEta, e.tkPhi, e.tk_lostInnerHits, e.tk_pval, e.tk_sigIP3D_vtxB, e.massVis_wTk, e.massHad_wTk, e.massMuTk, e.massDTk, e.mass2MissTk, e.UmissTk]
        auxList += [e.MC_tkFlag, e.MC_tkFromMainB, e.MC_tkPdgId, e.MC_tkMotherPdgId, e.MC_tkMotherMotherPdgId, e.MC_tk_dphi, e.MC_tk_deta, e.MC_tk_dpt]
        for l in auxList:
            l += [0, 0, 0]


    #----------------- Additional neutrals -------------------#
    idx_st = 0
    for jjj in range(j):
        idx_st += int(ev.nNeuAdd[jjj])
    idx_stop = int(idx_st + ev.nNeuAdd[j])

    e.N_goodAddNeu = 0
    e.neuPdgId = []
    e.neuPt = []
    e.neuEta = []
    e.neuPhi = []
    e.neuEnergy = []
    e.neuEt2 = []
    e.massVis_wNeu = []
    e.massHad_wNeu = []
    e.massMuNeu = []

    p4_sumGoodNeu = rt.TLorentzVector()
    for jj in range(idx_st, idx_stop):

        eta = ev.neuAdd_eta[jj]
        if np.abs(eta) >= 2.4:
            continue
        phi = ev.neuAdd_phi[jj]
        pt = ev.neuAdd_pt[jj]

        pdgId = ev.neuAdd_pdgId[jj]

        #Avoid duplicates
        duplicate = False
        for n in ['mu', 'pi', 'K', 'pis']:
            dphi = phi - getattr(e, n+'_phi')
            if dphi > np.pi: dphi -= 2*np.pi
            if dphi < -np.pi: dphi += 2*np.pi
            dR = np.hypot(dphi, eta - getattr(e, n+'_eta'))
            dPt = np.abs(getattr(e, n+'_pt') - pt)/getattr(e, n+'_pt')
            if dPt < 0.03 and dR < 0.001:
                duplicate=True
                break
        if duplicate:
            continue

        if pdgId == 22: # photon
            m_neu = 0
        elif pdgId == 130: # Neutral hadron, a K_L
            m_neu = 0.497611
        else:
            print 'Pdg not recognized in add neutral:', pdgId
            raise
        p4_neu = rt.TLorentzVector()
        p4_neu.SetPtEtaPhiM(pt, eta, phi, m_pi)

        mVis_wNeu = (p4_vis + p4_neu).M()

        e.N_goodAddNeu += 1
        idx, e.neuPt = insertOrdered(e.neuPt, pt)
        e.neuEta.insert(idx, eta)
        e.neuPhi.insert(idx, phi)
        e.neuPdgId.insert(idx, pdgId)
        e.neuEnergy.insert(idx, ev.neuAdd_energy[jj])
        e.neuEt2.insert(idx, ev.neuAdd_et2[jj])
        e.massVis_wNeu.insert(idx, mVis_wNeu)
        e.massHad_wNeu.insert(idx, (p4_Dst + p4_neu).M())
        e.massMuNeu.insert(idx, (p4_mu + p4_neu).M())

        p4_sumGoodNeu += p4_neu


    p4_vis_wNeu = p4_vis + p4_sumGoodNeu
    e.massVisNeu = p4_vis_wNeu.M()
    e.massHadNeu = (p4_Dst + p4_sumGoodNeu).M()

    if e.N_goodAddNeu < 3:
        auxList = [e.neuPdgId, e.neuPt, e.neuEta, e.neuPhi, e.neuEnergy, e.neuEt2, e.massVis_wNeu, e.massHad_wNeu, e.massMuNeu]
        for l in auxList:
            l += [0, 0, 0]

    return e

def makeSelection(inputs):
    n, tag, filepath, leafs_names, cat, idxInt, corr, skipCut, trkControlRegion, serial = inputs
    N_accepted_cand = []
    N_accepted_tot = 0

    start, stop = idxInt

    tree = rt.TChain('outA/Tevts')
    for fn in glob(filepath):
        tree.Add(fn)

    if serial:
        pb = ProgressBar(maxEntry=stop)
    else:
        perc = int((stop-start)*0.35)

    output = np.zeros((int(1.5*(stop-start+1)), len(leafs_names)))

    for i_ev in range(start,stop):
        bytesRead = tree.GetEntry(i_ev, 1)
        if bytesRead == 0:
            print 'Error reading entry', i_ev
        ev = tree

        if serial:
            pb.show(i_ev-start)
        elif (i_ev-start) % perc == 0:
            print tag, ': {:.0f}%'.format(100*(i_ev+1-start)/(stop-start))
        N_acc = 0

        ev_output = []
        for j in range(ev.pval_piK.size()):
            idxTrg = int(ev.mu_trgMu_idx[j]) if hasattr(ev, 'mu_trgMu_idx') else int(ev.mu_trgCand_idx[j])
            evEx = extractEventInfos(j, ev, corr)

            if not cat is None:
                if not trigger_selection(idxTrg, ev, evEx, cat):
                    continue

            if not skipCut == 'all':
                if not candidate_selection(j, ev, evEx, skipCut, trkControlRegion):
                    continue

            N_acc += 1

            aux = (ev.runNum, ev.lumiNum, ev.eventNum, ev.LumiBlock,
                   evEx.q2, evEx.Est_mu, evEx.M2_miss, evEx.U_miss,
                   evEx.q2_coll, evEx.Est_mu_coll, evEx.M2_miss_coll,
                   ev.mu_charge[j], evEx.mu_pt, evEx.mu_eta, evEx.mu_phi, ev.trgMu_sigdxy_BS[idxTrg],
                   ev.mu_dca_vtxDst[j], ev.mu_sigdca_vtxDst[j],
                   ev.mu_dcaT_vtxDst[j], ev.mu_sigdcaT_vtxDst[j],
                   ev.mu_dca_vtxDstMu[j], ev.mu_sigdca_vtxDstMu[j],
                   ev.mu_dcaT_vtxDstMu[j], ev.mu_sigdcaT_vtxDstMu[j],
                   ev.mu_lostInnerHits[j],
                   ev.mu_kickFinder[j], ev.mu_segmentCompatibility[j],
                   ev.mu_trackerStandalonePosLocalChi2[j], ev.mu_tightId[j],
                   ev.mu_IP3D_vtxDst[j], ev.mu_sigIP3D_vtxDst[j],
                   ev.mu_db_iso04[j], ev.mu_db_corr_iso04[j], ev.mu_db_iso03[j],
                   ev.mu_db_iso04[j]/evEx.mu_pt, ev.mu_db_corr_iso04[j]/evEx.mu_pt, ev.mu_db_iso03[j]/evEx.mu_pt,
                   evEx.B_pt, evEx.B_eta, evEx.B_phi,
                   evEx.Dst_pt, evEx.Dst_eta, evEx.Dst_phi,
                   evEx.D0_pt, evEx.D0_eta, evEx.D0_phi,
                   evEx.pi_pt, evEx.pi_eta, evEx.pi_phi,
                   ev.pi_lostInnerHits[j],
                   evEx.K_pt, evEx.K_eta, evEx.K_phi,
                   ev.K_lostInnerHits[j],
                   ev.pval_piK[j], ev.sigdxy_vtxD0_PV[j],
                   evEx.pis_pt, evEx.pis_eta, evEx.pis_phi,
                   ev.pis_lostInnerHits[j],
                   ev.pval_D0pis[j],
                   evEx.mass_piK, evEx.mass_D0pis, evEx.mass_D0pismu,
                   evEx.mass_D0pis - evEx.mass_piK,
                   evEx.mass_D0pismu_muASpi, evEx.mass_D0pismu_muASK,
                   evEx.D0pismu_eta, evEx.D0pismu_phi,
                   ev.pval_D0pismu[j], ev.chi2_D0pismu[j],
                   ev.d_vtxD0pismu_PV[j], ev.dxy_vtxD0pismu_PV[j],
                   ev.cos_D0pismu_PV[j], ev.cosT_D0pismu_PV[j],
                   evEx.N_goodAddTks,
                   evEx.tkCharge[0], evEx.tkCharge[1], evEx.tkCharge[2],
                   evEx.tkPdgId[0], evEx.tkPdgId[1], evEx.tkPdgId[2],
                   evEx.tkPt[0], evEx.tkPt[1], evEx.tkPt[2],
                   evEx.tkPtError[0], evEx.tkPtError[1], evEx.tkPtError[2],
                   evEx.tkEta[0], evEx.tkEta[1], evEx.tkEta[2],
                   evEx.tkPhi[0], evEx.tkPhi[1], evEx.tkPhi[2],
                   evEx.tk_lostInnerHits[0], evEx.tk_lostInnerHits[1], evEx.tk_lostInnerHits[2],
                   evEx.tk_pval[0], evEx.tk_pval[1], evEx.tk_pval[2],
                   evEx.tk_sigIP3D_vtxB[0], evEx.tk_sigIP3D_vtxB[1], evEx.tk_sigIP3D_vtxB[2],
                   evEx.massVis_wTk[0], evEx.massVis_wTk[1], evEx.massVis_wTk[2],
                   evEx.massHad_wTk[0], evEx.massHad_wTk[1], evEx.massHad_wTk[2],
                   evEx.massMuTk[0], evEx.massMuTk[1], evEx.massMuTk[2],
                   evEx.massDTk[0], evEx.massDTk[1], evEx.massDTk[2],
                   evEx.mass2MissTk[0], evEx.mass2MissTk[1], evEx.mass2MissTk[2],
                   evEx.UmissTk[0], evEx.UmissTk[1], evEx.UmissTk[2],
                   evEx.massTks_pipi, evEx.massTks_KK, evEx.massTks_piK, evEx.massTks_Kpi,
                   evEx.massVisTks, evEx.massHadTks,
                   evEx.massVisTks1, evEx.massVisTks2, evEx.massVisTks3,
                   evEx.massHadTks1, evEx.massHadTks2, evEx.massHadTks3,
                   evEx.massHadTks_DstMassConstraint,
                   evEx.M2missTks, evEx.q2Tks,
                   evEx.EmissTks, evEx.PmissTks, evEx.UmissTks,
                   evEx.N_goodAddNeu,
                   evEx.neuPdgId[0], evEx.neuPdgId[1], evEx.neuPdgId[2],
                   evEx.neuPt[0], evEx.neuPt[1], evEx.neuPt[2],
                   evEx.neuEta[0], evEx.neuEta[1], evEx.neuEta[2],
                   evEx.neuPhi[0], evEx.neuPhi[1], evEx.neuPhi[2],
                   evEx.neuEnergy[0], evEx.neuEnergy[1], evEx.neuEnergy[2],
                   evEx.neuEt2[0], evEx.neuEt2[1], evEx.neuEt2[2],
                   evEx.massVis_wNeu[0], evEx.massVis_wNeu[1], evEx.massVis_wNeu[2],
                   evEx.massHad_wNeu[0], evEx.massHad_wNeu[1], evEx.massHad_wNeu[2],
                   evEx.massMuNeu[0], evEx.massMuNeu[1], evEx.massMuNeu[2],
                   evEx.massVisNeu, evEx.massHadNeu,
                   trigger_selection(idxTrg, ev, evEx, categories['low']),
                   trigger_selection(idxTrg, ev, evEx, categories['mid']),
                   trigger_selection(idxTrg, ev, evEx, categories['high']),
                   ev.trgMu_HLT_Mu12_IP6[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu12_IP6') else ev.trgObj_HLT_Mu12_IP6[idxTrg],
                   ev.trgMu_HLT_Mu9_IP6[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu9_IP6') else ev.trgObj_HLT_Mu9_IP6[idxTrg],
                   ev.trgMu_HLT_Mu7_IP4[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu7_IP4') else ev.trgObj_HLT_Mu7_IP4[idxTrg],
                   ev.N_vertexes, ev.N_goodVtx[j], ev.PV_chi2[j], ev.PV_ndof[j], ev.localVertexDensity[j],
                   ev.localVertexDensity_10mm[j], ev.localVertexDensity_5mm[j], ev.localVertexDensity_1mm[j],
                   ev.localVertexDensity_cos800[j], ev.localVertexDensity_cos990[j], ev.localVertexDensity_cos999[j],
                   ev.vtx_PV_x[j], ev.vtx_PV_y[j], ev.vtx_PV_z[j],
                   ev.vtx_B_decay_x[j], ev.vtx_B_decay_y[j], ev.vtx_B_decay_z[j],
                   ev.beamSpot_x, ev.beamSpot_y, ev.beamSpot_z,
                  )
            if not 'data' in n:
                muSisPdgId = []
                for id in ev.MC_muSistersPdgId:
                    if np.abs(id) in [14, 16]: continue #neutrino
                    muSisPdgId.append(id)
                while len(muSisPdgId) < 2:
                    muSisPdgId.append(0)
                if abs(muSisPdgId[0]) < abs(muSisPdgId[1]):
                    auxSwap = muSisPdgId[0]
                    muSisPdgId[0] = muSisPdgId[1]
                    muSisPdgId[1] = auxSwap

                aux += (
                        ev.MC_nAddOgB, ev.MC_bestBB_dR, ev.MC_bestBB_dphi, ev.MC_bestBB_mass,
                        ev.MC_q2, ev.MC_Est_mu, ev.MC_M2_miss,
                        ev.MC_B_pt, ev.MC_B_eta, ev.MC_B_phi, ev.MC_B_ctau,
                        ev.MC_Dst_pt, ev.MC_Dst_eta, ev.MC_Dst_phi,
                        ev.MC_mu_pt, ev.MC_mu_eta, ev.MC_mu_phi,
                        ev.MC_mu_TransvIP_PV, ev.MC_mu_TransvIP_vtxDst, ev.MC_mu_IP_vtxDst,
                        ev.MC_pi_pt, ev.MC_pi_eta, ev.MC_pi_phi,
                        ev.MC_K_pt, ev.MC_K_eta, ev.MC_K_phi,
                        ev.MC_pis_pt, ev.MC_pis_eta, ev.MC_pis_phi,
                        ev.MC_idxCand == j,
                        ev.MC_muMotherPdgId,
                        muSisPdgId[0], muSisPdgId[1],
                        ev.MC_MassCharmedBDaugther,
                        ev.MC_DstMotherPdgId, ev.MC_CharmedDstSisPdgId, ev.MC_StrangeDstSisPdgId,
                        ev.MC_nAddCharged, ev.MC_addCharged_SumQ, ev.MC_nAddNeutral,
                        evEx.MC_tkFlag[0], evEx.MC_tkFlag[1],
                        evEx.MC_tkFromMainB[0], evEx.MC_tkFromMainB[1],
                        evEx.MC_tk_dpt[0], evEx.MC_tk_dpt[1],
                        evEx.MC_tk_deta[0], evEx.MC_tk_deta[1],
                        evEx.MC_tk_dphi[0], evEx.MC_tk_dphi[1],
                        evEx.MC_tkPdgId[0], evEx.MC_tkPdgId[1],
                        evEx.MC_tkMotherPdgId[0], evEx.MC_tkMotherPdgId[1],
                        evEx.MC_tkMotherMotherPdgId[0], evEx.MC_tkMotherMotherPdgId[1],
                        ev.nTrueIntMC
                       )
            if n in ['Bd_MuNuDst', 'Bd_TauNuDst']:
                id_0, id_1 = np.abs(muSisPdgId).astype(np.int)[:2]
                m_Dst = np.abs(int(ev.MC_DstMotherPdgId))
                m_mu = np.abs(int(ev.MC_muMotherPdgId))
                decayRawList = list(np.sort(np.abs(ev.MC_decay)[4:])[::-1].astype(np.int))
                if (m_mu == 15) and (decayRawList == [413, 16, 15]) and (m_Dst == 511):
                    process = 0
                elif id_0 == 413 and (id_1 in [0, 22]) and (m_mu == m_Dst) and (m_mu == 511):
                    process = 1
                elif len(ev.MC_decay) == 0:
                    process = -1
                else:
                    process = -2

                aux += (ev.wh_CLNCentral,
                        ev.wh_CLNR0Down, ev.wh_CLNR0Up,
                        ev.wh_CLNeig1Down, ev.wh_CLNeig1Up,
                        ev.wh_CLNeig2Down, ev.wh_CLNeig2Up,
                        ev.wh_CLNeig3Down, ev.wh_CLNeig3Up,
                        ev.wh_BLPRCentral,
                        ev.wh_BLPReig1Down, ev.wh_BLPReig1Up,
                        ev.wh_BLPReig2Down, ev.wh_BLPReig2Up,
                        ev.wh_BLPReig3Down, ev.wh_BLPReig3Up,
                        ev.wh_BLPReig4Down, ev.wh_BLPReig4Up,
                        ev.wh_BLPReig5Down, ev.wh_BLPReig5Up,
                        ev.wh_BLPReig6Down, ev.wh_BLPReig6Up,
                        process
                       )
            if 'MuNuDstPi' in n:
                aux += (ev.wh_Dstst_BLRCentral,
                        ev.wh_DststN_BLReig1Down, ev.wh_DststN_BLReig1Up,
                        ev.wh_DststN_BLReig2Down, ev.wh_DststN_BLReig2Up,
                        ev.wh_DststN_BLReig3Down, ev.wh_DststN_BLReig3Up,
                        ev.wh_DststN_BLReig4Down, ev.wh_DststN_BLReig4Up,
                        ev.wh_DststW_BLReig1Down, ev.wh_DststW_BLReig1Up,
                        ev.wh_DststW_BLReig2Down, ev.wh_DststW_BLReig2Up,
                        ev.wh_DststW_BLReig3Down, ev.wh_DststW_BLReig3Up,
                        ev.wh_D2S_BLOPCentral,
                        ev.wh_D2S_BLOPRhoSqDown, ev.wh_D2S_BLOPRhoSqUp,
                        ev.wh_D2S_BLOPchi11Down, ev.wh_D2S_BLOPchi11Up,
                        ev.wh_D2S_BLOPchi21Down, ev.wh_D2S_BLOPchi21Up,
                        ev.wh_D2S_BLOPchi31Down, ev.wh_D2S_BLOPchi31Up,
                        ev.wh_D2S_BLOPeta1Down, ev.wh_D2S_BLOPeta1Up,
                       )

                id_0, id_1 = np.abs(muSisPdgId).astype(np.int)[:2]
                id_m = np.abs(int(ev.MC_DstMotherPdgId))
                id_muM = np.abs(int(ev.MC_muMotherPdgId))
                # B -> D1 mu nu
                if (id_0 in [10423, 10413]) and (id_1 in [0, 22]) and id_0 == id_m:
                    process = 1
                # B -> D1* mu nu
                elif (id_0 in [20423, 20413]) and (id_1 in [0, 22]) and id_0 == id_m:
                    process = 2
                # B -> D2* mu nu
                elif (id_0 in [425, 415]) and (id_1 in [0, 22]) and id_0 == id_m:
                    process = 3
                # B -> D*(2S) mu nu
                elif (id_0 in [100423, 100413]) and (id_1 in [0, 22]) and id_0 == id_m:
                    process = 4
                # B -> D*(2S) mu nu with D*(2S) -> D1 pi
                elif (id_0 in [100423, 100413]) and (id_1 in [0, 22]) and id_m in [10423, 10413]:
                    process = 11
                # B -> D*(2S) mu nu with D*(2S) -> D1* pi
                elif (id_0 in [100423, 100413]) and (id_1 in [0, 22]) and id_m in [20423, 20413]:
                    process = 12
                # B -> D*(2S) mu nu with D*(2S) -> D*2 pi
                elif (id_0 in [100423, 100413]) and (id_1 in [0, 22]) and id_m in [425, 415]:
                    process = 13
                # B -> D(2S) mu nu with D(2S) -> D*2 pi
                elif (id_0 in [100421, 100411]) and (id_1 in [0, 22]) and id_m in [425, 415]:
                    process = 23
                # B -> D1 pi mu nu
                elif (id_0 in [10423, 10413]) and (id_1 in [111, 211]) and id_0 == id_m:
                    process = 31
                # B -> D1* pi mu nu
                elif (id_0 in [20423, 20413]) and (id_1 in [111, 211]) and id_0 == id_m:
                    process = 32
                # B -> D2* pi mu nu
                elif (id_0 in [425, 415]) and (id_1 in [111, 211]) and id_0 == id_m:
                    process = 33
                # B -> D* (1+)pi mu nu
                elif id_0 == 413 and (id_1 in [111, 211]) and (id_m in [511, 521]):
                    process = 0
                # No MC decay matching
                elif len(ev.MC_decay) == 0:
                    process = -1
                # B -> D* mu nu (slipped in)
                elif id_0 == 413 and (id_1 in [0, 22]) and (id_m in [511, 521]):
                    process = -2
                # B -> D* tau nu (slipped in)
                elif list(np.sort(np.abs(ev.MC_decay)[2:])[::-1].astype(np.int)) == [511, 413, 16, 15, 0]:
                    process = -3
                # Muon from Tau
                elif id_muM == 15:
                    process = -4
                # Muon from D meson
                elif (id_muM %1000)/100 == 4:
                    process = -5
                # Muon from K meson
                elif (id_muM %1000)/100 == 3:
                    process = -6
                else:
                    print '\n\n'
                    print id_0, id_1
                    print id_m
                    print id_muM
                    print [x for x in ev.MC_decay]
                    raise

                aux += (process,)
            if re.match('B[usd]_DstD[usd]', n):
                id_0, id_1 = np.abs(muSisPdgId).astype(np.int)[:2]
                m_Dst = np.abs(int(ev.MC_DstMotherPdgId))
                m_mu = np.abs(int(ev.MC_muMotherPdgId))
                sis_c = np.abs(ev.MC_CharmedDstSisPdgId)
                sis_s = np.abs(ev.MC_StrangeDstSisPdgId)

                decayRawList = list(np.sort(np.abs(ev.MC_decay)[4:])[::-1].astype(np.int))

                if (m_Dst == 511 and sis_c == 421 and sis_s == 321):
                    process = 101
                elif (m_Dst == 511 and sis_c == 421 and sis_s == 323):
                    process = 102
                elif (m_Dst == 511 and sis_c == 423 and sis_s == 321):
                    process = 103
                elif (m_Dst == 511 and sis_c == 423 and sis_s == 323):
                    process = 104
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 311):
                    process = 105 # 1.5, 1.8
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 313):
                    process = 106 # 1.6, 1.9
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 0):
                    process = 107 # 1.7, 1.10

                elif (m_Dst == 511 and sis_c == 411 and sis_s == 311):
                    process = 201 # 2.1, 2.2
                elif (m_Dst == 511 and sis_c == 411 and sis_s == 313):
                    process = 203 # 2.3, 2.4
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 311):
                    process = 205 # 2.5, 2.8
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 313):
                    process = 206 # 2.6, 2.9
                elif (m_Dst == 511 and sis_c == 413 and sis_s == 0):
                    process = 207 # 2.7, 2.10
                elif (m_Dst == 511 and sis_c == 411 and sis_s == 0):
                    process = 211 # 2.11, 2.12

                elif (m_Dst == 511 and sis_c == 431 and sis_s == 0):
                    process = 301
                elif (m_Dst == 511 and sis_c == 433 and sis_s == 0):
                    process = 302
                elif (m_Dst == 511 and sis_c == 10431 and sis_s == 0):
                    process = 303

                elif (m_Dst == 521 and sis_c == 421 and sis_s == 311):
                    process = 401
                elif (m_Dst == 521 and sis_c == 421 and sis_s == 313):
                    process = 402
                elif (m_Dst == 521 and sis_c == 423 and sis_s == 311):
                    process = 403
                elif (m_Dst == 521 and sis_c == 423 and sis_s == 313):
                    process = 404
                elif (m_Dst == 521 and sis_c == 413 and sis_s == 321):
                    process = 405 # 4.5, 4.7
                elif (m_Dst == 521 and sis_c == 413 and sis_s == 323):
                    process = 406 # 4.6, 4.8
                elif (m_Dst == 521 and sis_c == 421 and sis_s == 0):
                    process = 409
                elif (m_Dst == 521 and sis_c == 423 and sis_s == 0):
                    process = 410

                # Not present in out MC
                # B+ -> D* D_s (Cabibbo soppressed)
                elif (m_Dst == 521 and sis_c == 431 and sis_s == 0):
                    process = 450

                elif (m_Dst == 521 and sis_c == 411 and sis_s == 321):
                    process = 501 # 5.1, 5.3
                elif (m_Dst == 521 and sis_c == 411 and sis_s == 323):
                    process = 502 # 5.2, 5.4
                elif (m_Dst == 521 and sis_c == 413 and sis_s == 321):
                    process = 505 # 5.5, 5.7
                elif (m_Dst == 521 and sis_c == 413 and sis_s == 323):
                    process = 506 # 5.6, 5.8
                # Not present in our MC
                # B+ -> D* D_s* (Cabibbo soppressed)
                elif (m_Dst == 521 and sis_c == 433 and sis_s == 0):
                    process = 520

                elif (m_Dst == 531 and sis_c == 433 and sis_s == 311):
                    process = 601
                elif (m_Dst == 531 and sis_c == 433 and sis_s == 313):
                    process = 602
                elif (m_Dst == 531 and sis_c == 431 and sis_s == 311):
                    process = 603
                elif (m_Dst == 531 and sis_c == 431 and sis_s == 313):
                    process = 604
                elif (m_Dst == 531 and sis_c == 431 and sis_s == 0):
                    process = 605
                elif (m_Dst == 531 and sis_c == 433 and sis_s == 0):
                    process = 606
                # Not present in out MC
                elif (m_Dst == 531 and sis_c == 421 and sis_s == 321):
                    process = 620
                elif (m_Dst == 531 and sis_c == 411 and sis_s == 311):
                    process = 621

                elif (m_Dst == 511 and sis_c == 20433 and sis_s == 0):
                    process = 701
                elif (m_Dst == 511 and sis_c == 10433 and sis_s == 0):
                    process = 702

                # No MC decay matching
                elif len(ev.MC_decay) == 0:
                    process = -1
                # B -> D* mu nu (slipped in)
                elif id_0 == 413 and (id_1 in [0, 22, 111, 211]) and (m_Dst in [511, 521]) and (m_mu == m_Dst):
                    process = -2
                # B -> D* tau nu (slipped in)
                elif decayRawList == [413, 16, 15]:
                    process = -3
                # B -> D** mu nu (slipped in)
                elif (m_mu in [511, 521]) and (id_0 in [10413, 10423, 20413, 20423, 415, 425]) and (id_1 in [0, 22]) and (m_Dst == id_0):
                    process = -4
                # Bs -> Ds** mu nu (slipped in)
                elif (m_mu == 531) and (id_0 in [10433, 435]) and (id_1 in [0, 22]) and (m_Dst == id_0):
                    process = -5
                # D** tau nu splipped in
                elif (m_mu == 15) and (m_Dst in [10413, 10423, 20413, 20423, 415, 425, 10433]) and decayRawList == [m_Dst, 16, 15]:
                    process = -6
                # 4 or more bodies decay
                elif len(decayRawList) > 3:
                    process = -100
                # D* from D**
                elif m_Dst in [10413, 10423, 20413, 20423, 415, 425, 10433]:
                    process = -101
                else:
                    print '\n\nUnrecognized decay'
                    print 'Mu sisters:', id_0, id_1
                    print 'Mu mother:', m_mu
                    print 'Dst sisters:', sis_c, sis_s
                    print 'Dst mother:', m_Dst
                    print 'Raw decay list:', decayRawList
                    print ' '.join([str(int(x)) for x in ev.MC_decay])
                    raise

                aux += (process,)
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

def create_dSet(n, filepath, cat, applyCorrections=False, skipCut=[], trkControlRegion=False, maxEvents=1e15):
    if cat is None:
        catName = 'NoCat'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print n, catName
    if 'data' in n:
        loc = join(root,'cmsRD/skimmed'+args.skimTag+'/B2DstMu'+ n.replace('data', ''))
        out = re.search('2[12][01][0-9][0-3][0-9]', filepath)
        if out is None:
            print filepath
            raise
        fskimmed_name = loc + '_' + out.group(0) + '_' + catName
        N_evts_per_job = 75000
    else:
        d = join(os.path.dirname(filepath),'skimmed'+args.skimTag+'/')
        if not os.path.isdir(d):
            os.makedirs(d)
        fskimmed_name = d + catName
        N_evts_per_job = 20000 # looser tracks
        if n == 'Bd_MuNuDst':
            N_evts_per_job = 40000
    if not skipCut == []:
        print 'Skipping cut(s)', skipCut
        if skipCut == 'all':
            fskimmed_name += '_skipall'
        else:
            fskimmed_name += '_skip'+'-'.join([str(i) for i in skipCut])
    if trkControlRegion:
        print 'Track control region'
        fskimmed_name += '_trkCtrl'
        N_evts_per_job *= 2
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
        hAllNvtx = rt.TH1D('hAllNvtx', 'hAllNvtx', 101, -0.5, 100.5)
        hAllVtxZ = rt.TH1D('hAllVtxZ', 'hAllVtxZ', 100, -25, 25)
        hAllNTrueIntMC = rt.TH1D('hAllNTrueIntMC', 'hAllNTrueIntMC', 101, -0.5, 100.5)

        filenames = glob(filepath)
        print "Analyzing %i files matching '%s'" % (len(filenames),filepath)
        pb = ProgressBar(maxEntry=len(filenames))
        nTryMax = 10 if 'test' in args.skimTag else 500
        print 'Merging vertexes histograms from the first {} files.'.format(nTryMax)
        for i, fn in enumerate(filenames):
            pb.show(i)
            try:
                tree.Add(fn)
                if i < nTryMax:
                    fAux = rt.TFile.Open(fn, 'READ')
                    hAux = fAux.Get('trgF/hAllNvts')
                    hAllNvtx.Add(hAux)
                    hAux = fAux.Get('trgF/hAllVtxZ')
                    hAllVtxZ.Add(hAux)
                    if not 'data' in n:
                        hAux = fAux.Get('trgF/hAllNTrueIntMC')
                        hAllNTrueIntMC.Add(hAux)
                    fAux.Close()
            except:
                print >> sys.stderr, '[ERROR] Problem with vertexes histograms in %s' % fn
                raise
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEvents, tree.GetEntries())
        print n, ': Total number of candidate events =', N_cand_in

        leafs_names = ['runNum', 'lumiNum', 'eventNum', 'lumiBlock',
                       'q2', 'Est_mu', 'M2_miss', 'U_miss',
                       'q2_coll', 'Est_mu_coll', 'M2_miss_coll',
                       'mu_charge', 'mu_pt', 'mu_eta', 'mu_phi', 'mu_sigdxy',
                       'mu_dca_vtxDst', 'mu_sigdca_vtxDst',
                       'mu_dcaT_vtxDst', 'mu_sigdcaT_vtxDst',
                       'mu_dca_vtxDstMu', 'mu_sigdca_vtxDstMu',
                       'mu_dcaT_vtxDstMu', 'mu_sigdcaT_vtxDstMu',
                       'mu_lostInnerHits',
                       'mu_kickFinder', 'mu_segmentCompatibility',
                       'mu_trackerStandalonePosLocalChi2', 'mu_tightId',
                       'mu_IP3D_vtxDst', 'mu_sigIP3D_vtxDst',
                       'mu_db_iso04', 'mu_db_corr_iso04', 'mu_db_iso03',
                       'mu_db_iso04_rel', 'mu_db_corr_iso04_rel', 'mu_db_iso03_rel',
                       'B_pt', 'B_eta', 'B_phi',
                       'Dst_pt', 'Dst_eta', 'Dst_phi',
                       'D0_pt', 'D0_eta', 'D0_phi',
                       'pi_pt', 'pi_eta', 'pi_phi',
                       'pi_lostInnerHits',
                       'K_pt', 'K_eta', 'K_phi',
                       'K_lostInnerHits',
                       'pval_piK', 'sigdxy_vtxD0_PV',
                       'pis_pt', 'pis_eta', 'pis_phi',
                       'pis_lostInnerHits',
                       'pval_D0pis',
                       'mass_piK', 'mass_D0pis', 'mass_D0pismu',
                       'deltaM_DstD',
                       'mass_D0pismu_muASpi', 'mass_D0pismu_muASK',
                       'D0pismu_eta', 'D0pismu_phi',
                       'pval_D0pismu', 'chi2_D0pismu',
                       'd_vtxD0pismu_PV', 'dxy_vtxD0pismu_PV',
                       'cos_D0pismu_PV', 'cosT_D0pismu_PV',
                       'N_goodAddTks',
                       'tkCharge_0', 'tkCharge_1', 'tkCharge_2',
                       'tkPdgId_0', 'tkPdgId_1', 'tkPdgId_2',
                       'tkPt_0', 'tkPt_1', 'tkPt_2',
                       'tkPtError_0', 'tkPtError_1', 'tkPtError_2',
                       'tkEta_0', 'tkEta_1', 'tkEta_2',
                       'tkPhi_0', 'tkPhi_1', 'tkPhi_2',
                       'tk_lostInnerHits_0', 'tk_lostInnerHits_1', 'tk_lostInnerHits_2',
                       'tk_pval_0', 'tk_pval_1', 'tk_pval_2',
                       'tk_sigIP3D_vtxB_0', 'tk_sigIP3D_vtxB_1', 'tk_sigIP3D_vtxB_2',
                       'tkMassVis_0', 'tkMassVis_1', 'tkMassVis_2',
                       'tkMassHad_0', 'tkMassHad_1', 'tkMassHad_2',
                       'tkMassMuTk_0', 'tkMassMuTk_1', 'tkMassMuTk_2',
                       'tkMassDTk_0', 'tkMassDTk_1', 'tkMassDTk_2',
                       'tkMassMiss2_0', 'tkMassMiss2_1', 'tkMassMiss2_2',
                       'tkUmiss_0', 'tkUmiss_1', 'tkUmiss_2',
                       'massTks_pipi', 'massTks_KK', 'massTks_piK', 'massTks_Kpi',
                       'massVisTks', 'massHadTks',
                       'massVisTks1',
                       'massVisTks2',
                       'massVisTks3',
                       'massHadTks1',
                       'massHadTks2',
                       'massHadTks3',
                       'massHadTks_DstMassConstraint',
                       'M2missTks', 'q2Tks',
                       'EmissTks', 'PmissTks', 'UmissTks',
                       'N_goodAddNeu',
                       'neuPdgId_0', 'neuPdgId_1', 'neuPdgId_2',
                       'neuPt_0', 'neuPt_1', 'neuPt_2',
                       'neuEta_0', 'neuEta_1', 'neuEta_2',
                       'neuPhi_0', 'neuPhi_1', 'neuPhi_2',
                       'neuEnergy_0', 'neuEnergy_1', 'neuEnergy_2',
                       'neuEt2_0', 'neuEt2_1', 'neuEt2_2',
                       'massVis_wNeu_0', 'massVis_wNeu_1', 'massVis_wNeu_2',
                       'massHad_wNeu_0', 'massHad_wNeu_1', 'massHad_wNeu_2',
                       'massMuNeu_0', 'massMuNeu_1', 'massMuNeu_2',
                       'massVisNeu', 'massHadNeu',
                       'cat_low', 'cat_mid', 'cat_high',
                       'muPass_Mu12_IP6', 'muPass_Mu9_IP6', 'muPass_Mu7_IP4',
                       'N_vtx', 'N_goodVtx', 'PV_chi2', 'PV_ndof', 'localVertexDensity',
                       'localVertexDensity_10mm', 'localVertexDensity_5mm', 'localVertexDensity_1mm',
                       'localVertexDensity_cos800', 'localVertexDensity_cos990', 'localVertexDensity_cos999',
                       'vtx_PV_x', 'vtx_PV_y', 'vtx_PV_z',
                       'vtx_B_decay_x', 'vtx_B_decay_y', 'vtx_B_decay_z',
                       'beamSpot_x', 'beamSpot_y', 'beamSpot_z',
                      ]
        if not 'data' in n:
            leafs_names += [
                            'MC_nAddOgB', 'MC_bestBB_dR', 'MC_bestBB_dphi', 'MC_bestBB_mass',
                            'MC_q2', 'MC_Est_mu', 'MC_M2_miss',
                            'MC_B_pt', 'MC_B_eta', 'MC_B_phi', 'MC_B_ctau',
                            'MC_Dst_pt', 'MC_Dst_eta', 'MC_Dst_phi',
                            'MC_mu_pt', 'MC_mu_eta', 'MC_mu_phi',
                            'MC_mu_TransvIP_PV', 'MC_mu_TransvIP_vtxDst', 'MC_mu_IP_vtxDst',
                            'MC_pi_pt', 'MC_pi_eta', 'MC_pi_phi',
                            'MC_K_pt', 'MC_K_eta', 'MC_K_phi',
                            'MC_pis_pt', 'MC_pis_eta', 'MC_pis_phi',
                            'MC_idxMatch',
                            'MC_muMotherPdgId',
                            'MC_munuSisterPdgId_0', 'MC_munuSisterPdgId_1',
                            'MC_MassCharmedBDaughter',
                            'MC_DstMotherPdgId', 'MC_CharmedDstSisPdgId', 'MC_StrangeDstSisPdgId',
                            'MC_nAddCharged', 'MC_addCharged_SumQ', 'MC_nAddNeutral',
                            'MC_tkFlag_0', 'MC_tkFlag_1',
                            'MC_tkFromMainB_0', 'MC_tkFromMainB_1',
                            'MC_tk_dpt_0', 'MC_tk_dpt_1',
                            'MC_tk_deta_0', 'MC_tk_deta_1',
                            'MC_tk_dphi_0', 'MC_tk_dphi_1',
                            'MC_tkPdgId_0', 'MC_tkPdgId_1',
                            'MC_tkMotherPdgId_0', 'MC_tkMotherPdgId_1',
                            'MC_tkMotherMotherPdgId_0', 'MC_tkMotherMotherPdgId_1',
                            'MC_nInteractions'
                           ]
        if n in ['Bd_MuNuDst', 'Bd_TauNuDst']:
            leafs_names += ['wh_CLNCentral',
                            'wh_CLNR0Down', 'wh_CLNR0Up',
                            'wh_CLNeig1Down', 'wh_CLNeig1Up',
                            'wh_CLNeig2Down', 'wh_CLNeig2Up',
                            'wh_CLNeig3Down', 'wh_CLNeig3Up',
                            'wh_BLPRCentral',
                            'wh_BLPReig1Down', 'wh_BLPReig1Up',
                            'wh_BLPReig2Down', 'wh_BLPReig2Up',
                            'wh_BLPReig3Down', 'wh_BLPReig3Up',
                            'wh_BLPReig4Down', 'wh_BLPReig4Up',
                            'wh_BLPReig5Down', 'wh_BLPReig5Up',
                            'wh_BLPReig6Down', 'wh_BLPReig6Up',
                            'procId_signal'
                            ]
        elif 'MuNuDstPi' in n:
            leafs_names += ['wh_Dstst_BLRCentral',
                            'wh_DststN_BLReig1Down', 'wh_DststN_BLReig1Up',
                            'wh_DststN_BLReig2Down', 'wh_DststN_BLReig2Up',
                            'wh_DststN_BLReig3Down', 'wh_DststN_BLReig3Up',
                            'wh_DststN_BLReig4Down', 'wh_DststN_BLReig4Up',
                            'wh_DststW_BLReig1Down', 'wh_DststW_BLReig1Up',
                            'wh_DststW_BLReig2Down', 'wh_DststW_BLReig2Up',
                            'wh_DststW_BLReig3Down', 'wh_DststW_BLReig3Up',
                            'wh_D2S_BLOPCentral',
                            'wh_D2S_BLOPRhoSqDown', 'wh_D2S_BLOPRhoSqUp',
                            'wh_D2S_BLOPchi11Down', 'wh_D2S_BLOPchi11Up',
                            'wh_D2S_BLOPchi21Down', 'wh_D2S_BLOPchi21Up',
                            'wh_D2S_BLOPchi31Down', 'wh_D2S_BLOPchi31Up',
                            'wh_D2S_BLOPeta1Down', 'wh_D2S_BLOPeta1Up',
                            'procId_Dstst'
                            ]
        elif re.match('B[usd]_DstD[usd]', n):
            leafs_names += ['procId_DstHc']

        applyCorr = None
        if applyCorrections:
            applyCorr = 'MC'
            if 'data' in n:
                applyCorr = 'RD'

        if N_cand_in < 1.5*N_evts_per_job or args.parallelType == 'serial':
            output, N_accepted_cand = makeSelection([n, '', filepath, leafs_names, cat,
                                                     [0, N_cand_in-1], applyCorr, skipCut, trkControlRegion, True])
        else:
            pdiv = list(range(0, N_cand_in, N_evts_per_job))
            if not pdiv[-1] == N_cand_in:
                pdiv.append(N_cand_in)
            print 'Will be divided into ' + str(len(pdiv)-1) + ' jobs'
            inputs = []
            for i, (start, stop) in enumerate(zip(pdiv[:-1],pdiv[1:])):
                inputs.append([n, str(i), filepath, leafs_names, cat, [start, stop], applyCorr, skipCut, trkControlRegion, False])
            print ' '

            start = time.time()
            if args.parallelType == 'pool' or len(inputs) < 10:
                p = Pool(min(15,len(inputs)))
                outputs = p.map(makeSelection, inputs)
            elif args.parallelType == 'jobs':
                tmpDir = 'tmp/B2DstMu_skimCAND_%s_%s' % (n,catName)
                if trkControlRegion:
                    tmpDir += '_trkControl'
                os.system('rm -rf ' + tmpDir + '/out')
                os.system('rm -rf ' + tmpDir + '/*.p')
                os.makedirs(tmpDir + '/out')
                tmpDir = os.path.abspath(tmpDir)
                for ii, inAux in enumerate(inputs):
                    with open(join(tmpDir,'input_%i.p' % ii), 'wb') as f:
                        pickle.dump(inAux, f)
                createSubmissionFile(tmpDir, len(inputs))
                print 'Submitting jobs'
                cmd = 'condor_submit %s' % join(tmpDir,'jobs.jdl')
                batch_name = 'skim_%s_%s' % (n,catName)
                if trkControlRegion:
                    batch_name += '_trkControl'
                cmd += ' -batch-name %s' % batch_name
                status, output = commands.getstatusoutput(cmd)
                if status != 0:
                    print >> sys.stderr, "Error in processing command: '%s'" % cmd
                    print >> sys.stderr, "Output: %s" % output
                    sys.exit(1)
                print 'Job submitted'
                print 'Waiting for jobs to be finished'
                time.sleep(20)
                proceed = False
                while not proceed:
                    status, output = commands.getstatusoutput('condor_q')
                    found = False
                    for line in output.split('\n'):
                        if batch_name in line:
                            print line
                            time.sleep(10)
                            found = True
                    proceed = not found
                print 'All jobs finished'

                outputs = []
                print 'Fetching jobs output'
                pb = ProgressBar(len(inputs))
                for ii in range(len(inputs)):
                    pb.show(ii)
                    with open(join(tmpDir,'output_%i.p' % ii), 'rb') as f:
                        o = pickle.load(f)
                        outputs.append(o)

            print 'Concatenating the outputs'
            output = np.concatenate(tuple([o[0] for o in outputs]))
            N_accepted_cand = []
            for o in outputs: N_accepted_cand += o[1]
            print 'Total time: {:.1f} min'.format((time.time()-start)/60.)


        dset = pd.DataFrame(output, columns=leafs_names)
        if not os.path.isdir(os.path.dirname(fskimmed_name)):
            os.makedirs(os.path.dirname(fskimmed_name))
        rtnp.array2root(dset.to_records(), fskimmed_name, treename='Tevts', mode='RECREATE')
        fAux = rt.TFile.Open(fskimmed_name, 'UPDATE')
        hAllNvtx.Write()
        hAllVtxZ.Write()
        if not 'data' in n:
            hAllNTrueIntMC.Write()
        fAux.Close()

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
    job_file = join(tmpDir,'job.sh')
    with open(job_file, 'w') as fjob:
        fjob.write('#!/bin/bash\n')
        fjob.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        if os.environ['USER'] == 'ocerri':
            fjob.write('cd /storage/af/user/ocerri/CMSSW_10_2_3/; eval `scramv1 runtime -sh`\n')
            fjob.write('cd '+os.path.dirname(os.path.abspath(__file__))+'\n')
            fjob.write('python B2DstMu_skimCAND_v1.py --function makeSel --tmpDir $1 --jN $2\n')
        else:
            fjob.write('cd %s/RDstAnalysis/CMSSW_10_2_3/\n' % os.environ['HOME'])
            fjob.write('eval `scramv1 runtime -sh`\n')
            fjob.write('cd %s/RDstAnalysis/BPH_RD_Analysis/\n' % os.environ['HOME'])
            fjob.write('export PYTHONPATH=%s/RDstAnalysis/BPH_RD_Analysis/lib:$PYTHONPATH\n' % os.environ['HOME'])
            fjob.write('export PYTHONPATH=%s/RDstAnalysis/BPH_RD_Analysis/analysis:$PYTHONPATH\n' % os.environ['HOME'])
            fjob.write('python ./scripts/B2DstMu_skimCAND_v1.py --function makeSel --tmpDir $1 --jN $2\n')
    os.system('chmod +x {}/job.sh'.format(tmpDir))

    sub_file = join(tmpDir,'jobs.jdl')
    with open(sub_file, 'w') as fsub:
        fsub.write('executable    = %s\n' % join(tmpDir,'job.sh'))
        fsub.write('arguments     = {} $(ProcId)\n'.format(tmpDir))
        fsub.write('output        = {}/out/job_$(ProcId)_$(ClusterId).out\n'.format(tmpDir))
        fsub.write('error         = {}/out/job_$(ProcId)_$(ClusterId).err\n'.format(tmpDir))
        fsub.write('log           = {}/out/job_$(ProcId)_$(ClusterId).log\n'.format(tmpDir))
        fsub.write('WHEN_TO_TRANSFER_OUTPUT = ON_EXIT_OR_EVICT\n')
        fsub.write('+JobQueue="Normal"\n')
        # fsub.write('+RequestWalltime   = 7000\n')
        fsub.write('+MaxRuntime   = 7000\n')
        fsub.write('+RunAsOwner = True\n')
        fsub.write('+InteractiveUser = True\n')
        fsub.write('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"\n')
        fsub.write('+SingularityBindCVMFS = True\n')
        fsub.write('run_as_owner = True\n')
        fsub.write('RequestDisk = 2000000\n')
        fsub.write('RequestMemory = 2500\n')
        fsub.write('RequestCpus = 1\n')
        fsub.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        fsub.write('on_exit_remove = ((ExitBySignal == False) && (ExitCode == 0)) || (JobStatus=?=3)\n')
        # Send the job to Held state on failure.
        fsub.write('on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')
        # Periodically retry the jobs for 2 times with an interval of 20 minutes.
        fsub.write('periodic_release =  (NumJobStarts < 2) && ((CurrentTime - EnteredCurrentStatus) > (60*20))\n')
        fsub.write('+PeriodicRemove = ((JobStatus =?= 2) && ((MemoryUsage =!= UNDEFINED && MemoryUsage > 2.5*RequestMemory)))\n')
        fsub.write('max_retries    = 3\n')
        fsub.write('requirements   = Machine =!= LastRemoteHost\n')
        fsub.write('universe = vanilla\n')
        fsub.write('queue %i\n' % njobs)
        fsub.close()

if __name__ == "__main__":
    if args.function == 'main':
        file_loc = {}
        nCount = 0
        for n in args.dataset:
            for kn in filesLocMap:
                if re.match(n, kn):
                    print 'Adding %s' % kn
                    file_loc[kn] = filesLocMap[kn]
                    nCount += 1
        print 'Running over {} datasets'.format(nCount)

        if len(args.dataset) == 0:
            print >> sys.stderr, 'No dataset provided, rerun with -d'
            sys.exit(1)

        if len(file_loc) == 0:
            print >> sys.stderr, "No datasets found matching '%s'" % str(args.dataset)
            sys.exit(1)

        recreate = []
        if args.recreate:
            recreate = file_loc.keys()
        print '-'*50 + '\n'

        skip = []
        if args.skipCut == 'all':
            skip.append('all')
        elif args.skipCut:
            skip.append([int(args.skipCut)])
        else:
            skip.append([])

        trackControlFlag = []
        if args.region == 'all' or args.region=='signal':
            trackControlFlag.append(False)
        if args.region == 'all' or args.region=='trkControl':
            trackControlFlag.append(True)

        for idx in skip:
            for cn in args.cat:
                for iFile, (n, fp) in enumerate(file_loc.iteritems()):
                    print '>>>> Sample {}/{}'.format(iFile+1, len(file_loc.keys()))
                    for trkSelectionFlag in trackControlFlag:
                        create_dSet(n, fp, categories[cn], skipCut=idx, applyCorrections=args.applyCorr, trkControlRegion=trkSelectionFlag, maxEvents=args.maxEvents)

    elif args.function == 'makeSel':
        tmpDir = args.tmpDir
        with open(join(tmpDir,'input_%i.p' % args.jN), 'rb') as f:
            input = pickle.load(f)
        output = makeSelection(input)
        with open(join(tmpDir,'output_%i.p' % args.jN), 'wb') as f:
            pickle.dump(output, f)

    else:
        print args.function, 'not recognized'
