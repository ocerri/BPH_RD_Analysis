#!/usr/bin/env python
"""
Script for running the skimmer. After you have created the ntuples, you can run
this script by running:

    $ python skimmer.py ~/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples/B2DstMu/out*.root -o ~/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/skimmed/B2DstMu/ --applyCorr --selection B2DstMu

The argument to `-d` can be a regular expression, so if you want it to match
all the soft QCD MC events you can run:

    $ python skimmer.py -d SoftQCDnonD --applyCorr --selection B2DstMu

You can also pass a list of regular expressions:

    $ python skimmer.py -d BdToDstarMuNu,BdToDstarTauNu --applyCorr --selection B2DstMu
"""
import sys, os, pickle, time, re
from multiprocessing import Pool
import commands
from os.path import join
import numpy as np
import pandas as pd
from functools import partial
import uuid

try:
    from analysis_utilities import getEff
except ImportError:
    print >> sys.stderr, "Failed to import analysis_utilities."
    print >> sys.stderr, "Did you remember to source the env.sh file in the repo?"
    sys.exit(1)
from progressBar import ProgressBar
from categoriesDef import categories
import B02DstMu_selection
import B02JpsiKst_selection
import B02JpsiK_selection

import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)
import root_numpy as rtnp

root = join(os.environ['HOME'],'BPhysics/data')

class Container(object):
    pass

hBfieldMapsRatio = None

def get_bFieldCorr3D(phi, eta, verbose=False):
    global hBfieldMapsRatio
    if np.abs(eta) > 2.4:
        eta = 2.39*np.sign(eta)
    if np.abs(phi) > np.pi:
        phi = phi - 2*np.pi*np.sign(phi)
    if hBfieldMapsRatio is None:
        fBfield = rt.TFile.Open(join(root,'calibration/bFieldMap_2Dover3D.root'), 'r')
        hBfieldMapsRatio = fBfield.Get('bfieldMap')
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

def extractEventInfos_JpsiKst(j, ev, corr=None):
    m_mu   = 0.105658
    m_pi   = 0.139570
    m_K    = 0.493677
    m_Kst  = 0.8955
    m_jpsi = 3.096916
    m_B0   = 5.27963

    e = Container()
    # print '------> <-----'
    e.mup_eta = ev.mupRefit_eta[j]
    e.mup_phi = ev.mupRefit_phi[j]
    e.mup_pt = correctPt(ev.mupRefit_pt[j], e.mup_eta, e.mup_phi, corr, 3e-3)
    p4_mup = rt.TLorentzVector()
    p4_mup.SetPtEtaPhiM(e.mup_pt, e.mup_eta, e.mup_phi, m_mu)

    e.mum_eta = ev.mumRefit_eta[j]
    e.mum_phi = ev.mumRefit_phi[j]
    e.mum_pt = correctPt(ev.mumRefit_pt[j], e.mum_eta, e.mum_phi, corr, 3e-3)
    p4_mum = rt.TLorentzVector()
    p4_mum.SetPtEtaPhiM(e.mum_pt, e.mum_eta, e.mum_phi, m_mu)

    e.mass_mumu = compMass(e.mup_pt, e.mum_pt, e.mup_eta, e.mum_eta, e.mup_phi, e.mum_phi, m_mu, m_mu)
    e.Jpsi_pt = SumPt(e.mup_pt, e.mum_pt, e.mup_phi, e.mum_phi)

    e.K_eta = ev.KRefit_eta[j]
    e.K_phi = ev.KRefit_phi[j]
    e.K_pt = correctPt(ev.KRefit_pt[j], e.K_eta, e.K_phi, corr, 3e-3)
    p4_K = rt.TLorentzVector()
    p4_K.SetPtEtaPhiM(e.K_pt, e.K_eta, e.K_phi, m_K)

    e.pi_eta = ev.piRefit_eta[j]
    e.pi_phi = ev.piRefit_phi[j]
    e.pi_pt = correctPt(ev.piRefit_pt[j], e.pi_eta, e.pi_phi, corr, 3e-3)
    p4_pi = rt.TLorentzVector()
    p4_pi.SetPtEtaPhiM(e.pi_pt, e.pi_eta, e.pi_phi, m_pi)

    e.mass_piK = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_pi, m_K)
    # mass_piK_recomp = compMass(ev.piRefit_pt[j], ev.KRefit_pt[j], ev.piRefit_eta[j], ev.KRefit_eta[j], ev.piRefit_phi[j], ev.KRefit_phi[j], m_pi, m_K)
    e.Kst_pt = SumPt(e.pi_pt, e.K_pt, e.pi_phi, e.K_phi)
    e.Kst_eta = ev.KstRefit_eta[j]
    e.Kst_phi = ev.KstRefit_phi[j]

    p4_B = p4_mup + p4_mum + p4_K + p4_pi
    e.mass_mumupiK = p4_B.M()
    e.mass_mumupiK_cJpsi = p4_B.M() - e.mass_mumu + m_jpsi
    e.mass_mumupiK_cJpsi_cKst = p4_B.M() - e.mass_mumu + m_jpsi - e.mass_piK + m_Kst
    e.B_pt = p4_B.Pt()
    e.B_eta = p4_B.Eta()
    e.B_phi = p4_B.Phi()


    e.mass_Kpi = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_K, m_pi)
    # mass_Kpi_recomp = compMass(ev.piRefit_pt[j], ev.KRefit_pt[j], ev.piRefit_eta[j], ev.KRefit_eta[j], ev.piRefit_phi[j], ev.KRefit_phi[j], m_K, m_pi)
    p4_K.SetPtEtaPhiM(e.pi_pt, e.pi_eta, e.pi_phi, m_K)
    p4_pi.SetPtEtaPhiM(e.K_pt, e.K_eta, e.K_phi, m_pi)
    p4_B = p4_mup + p4_mum + p4_K + p4_pi
    e.mass_mumuKpi = p4_B.M()
    e.mass_mumuKpi_cJpsi = p4_B.M() - e.mass_mumu + m_jpsi
    e.mass_mumuKpi_cJpsi_cKst = p4_B.M() - e.mass_mumu + m_jpsi - e.mass_Kpi + m_Kst

    e.mass_KK = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_K, m_K)

    if ev.K_charge[j] > 0 and ev.pi_charge[j] < 0:
        e.isAntiB = 0
    elif ev.K_charge[j] < 0 and ev.pi_charge[j] > 0:
        e.isAntiB = 1
    else:
        print 'Charges not matching: K:{}  pi:{}'.format(ev.K_charge[j], ev.pi_charge[j])
        raise
    e.mass_candKst = e.mass_piK
    e.mass_candB = e.mass_mumupiK_cJpsi

    return e

def extractEventInfos_JpsiK(j, ev, corr=None):
    m_mu   = 0.105658
    m_pi   = 0.139570
    m_K    = 0.493677
    m_Kst  = 0.8955
    m_jpsi = 3.096916
    m_B0   = 5.27963

    e = Container()
    # print '------> <-----'
    e.mup_eta = ev.mupRefit_eta[j]
    e.mup_phi = ev.mupRefit_phi[j]
    e.mup_pt = correctPt(ev.mupRefit_pt[j], e.mup_eta, e.mup_phi, corr, 3e-3)
    p4_mup = rt.TLorentzVector()
    p4_mup.SetPtEtaPhiM(e.mup_pt, e.mup_eta, e.mup_phi, m_mu)

    e.mum_eta = ev.mumRefit_eta[j]
    e.mum_phi = ev.mumRefit_phi[j]
    e.mum_pt = correctPt(ev.mumRefit_pt[j], e.mum_eta, e.mum_phi, corr, 3e-3)
    p4_mum = rt.TLorentzVector()
    p4_mum.SetPtEtaPhiM(e.mum_pt, e.mum_eta, e.mum_phi, m_mu)

    e.mass_mumu = compMass(e.mup_pt, e.mum_pt, e.mup_eta, e.mum_eta, e.mup_phi, e.mum_phi, m_mu, m_mu)
    e.Jpsi_pt = SumPt(e.mup_pt, e.mum_pt, e.mup_phi, e.mum_phi)

    e.K_eta = ev.KRefit_eta[j]
    e.K_phi = ev.KRefit_phi[j]
    e.K_pt = correctPt(ev.KRefit_pt[j], e.K_eta, e.K_phi, corr, 3e-3)
    p4_K = rt.TLorentzVector()
    p4_K.SetPtEtaPhiM(e.K_pt, e.K_eta, e.K_phi, m_K)

    p4_B = p4_mup + p4_mum + p4_K
    e.mass_mumuK = p4_B.M()
    e.mass_mumuK_cJpsi = p4_B.M() - e.mass_mumu + m_jpsi

    e.B_pt = p4_B.Pt()
    e.B_eta = p4_B.Eta()
    e.B_phi = p4_B.Phi()

    #----------------- Additional tracks -------------------#
    idx_st = 0
    for jjj in range(j):
        idx_st += int(ev.nTksAdd[jjj])
    idx_stop = int(idx_st + ev.nTksAdd[j])

    e.N_goodAddTks = 0
    e.tkCharge = []
    e.tkPt = []
    e.tkEta = []
    e.tkPhi = []
    e.massVis_wTk = []
    e.massHad_wTk = []

    p4_sumGoodTks = rt.TLorentzVector()
    for jj in range(idx_st, idx_stop):
        eta = ev.tksAdd_eta[jj]
        if np.abs(eta) >= 2.4:
            continue
        phi = ev.tksAdd_phi[jj]
        pt = correctPt(ev.tksAdd_pt[jj], eta, phi, corr, 2e-3)
        if pt < 0.5:
            continue
        #Avoid tracks duplicates
        duplicate = False
        for n in ['mup', 'mum', 'K']:
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

        mVis_wTk = (p4_B + p4_tk).M()

        if mVis_wTk < m_B0 and ev.tksAdd_cos_PV[jj]>0.95:
            e.N_goodAddTks += 1
            idx, e.tkPt = insertOrdered(e.tkPt, pt)
            e.tkEta.insert(idx, eta)
            e.tkPhi.insert(idx, phi)
            e.tkCharge.insert(idx, ev.tksAdd_charge[jj])
            e.massVis_wTk.insert(idx, mVis_wTk)
            e.massHad_wTk.insert(idx, (p4_K + p4_tk).M())

            p4_sumGoodTks += p4_tk


    p4_vis_wTks = p4_B + p4_sumGoodTks
    e.massVisTks = p4_vis_wTks.M()
    e.massHadTks = (p4_K + p4_sumGoodTks).M()

    e.BwTks_pt = p4_vis_wTks.Pt()

    if e.N_goodAddTks < 3:
        auxList = [e.tkCharge, e.tkPt, e.tkEta, e.tkPhi, e.massVis_wTk, e.massHad_wTk]
        for l in auxList:
            l += [0, 0, 0]

    return e

def extractEventInfos_DstMu(j, ev, corr=None):
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
    e.tkPt = []
    e.tkPtError = []
    e.tkEta = []
    e.tkPhi = []
    e.MC_tkFlag = []
    e.MC_tkPdgId = []
    e.MC_tk_dphi = []
    e.MC_tk_deta = []
    e.MC_tk_dpt = []
    e.MC_tkMotherPdgId = []
    e.MC_tkMotherMotherPdgId = []
    e.massVis_wTk = []
    e.massHad_wTk = []
    e.massMuTk = []
    e.mass2MissTk = []
    e.UmissTk = []

    p4_sumGoodTks = rt.TLorentzVector()
    for jj in range(idx_st, idx_stop):
        eta = ev.tksAdd_eta[jj]
        if np.abs(eta) >= 2.4:
            continue
        phi = ev.tksAdd_phi[jj]
        pt = correctPt(ev.tksAdd_pt[jj], eta, phi, corr, 2e-3)
        if pt < 0.5:
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
            e.tkCharge.insert(idx, ev.tksAdd_charge[jj]*ev.mu_charge[j])
            e.massVis_wTk.insert(idx, mVis_wTk)
            e.massHad_wTk.insert(idx, (p4_Dst + p4_tk).M())
            e.massMuTk.insert(idx, (p4_mu + p4_tk).M())
            p_miss = p4_B - p4_vis - p4_tk
            e.mass2MissTk.insert(idx, p_miss.M2())
            e.UmissTk.insert(idx, p_miss.E() - p_miss.P())
            if hasattr(ev, 'MC_addTkFlag'):
                e.MC_tkFlag.insert(idx, ev.MC_addTkFlag[jj])
                e.MC_tk_dphi.insert(idx, ev.MC_addTk_dPhi[jj])
                e.MC_tk_deta.insert(idx, ev.MC_addTk_dEta[jj])
                e.MC_tk_dpt.insert(idx, ev.MC_addTk_dPt[jj])
                e.MC_tkPdgId.insert(idx, ev.MC_addTk_pdgId[jj])
                e.MC_tkMotherPdgId.insert(idx, ev.MC_addTk_pdgIdMother[jj])
                e.MC_tkMotherMotherPdgId.insert(idx, ev.MC_addTk_pdgIdMotherMother[jj])

            p4_sumGoodTks += p4_tk


        # if e.N_goodAddTks == 1:
        #     p4_tk1 = p4_tk.Clone('p4_tk1')
        # elif e.N_goodAddTks == 2:
        #     e.massVisTk12 = (p4_vis + p4_tk1 + p4_tk).M()
        #     e.massHadTk12 = (p4_Dst + p4_tk1 + p4_tk).M()
        #     p_miss = p4_B - p4_vis - p4_tk - p4_tk1
        #     e.UmissTk12 = p_miss.E() - p_miss.P()
    p4_vis_wTks = p4_vis + p4_sumGoodTks
    e.massVisTks = p4_vis_wTks.M()
    e.massHadTks = (p4_Dst + p4_sumGoodTks).M()
    e.massHadTks_DstMassConstraint = (p4_Dst + p4_sumGoodTks).M() - p4_Dst.M() + m_Dst

    e.BwTks_pt = p4_vis_wTks.Pt() * m_B0/ p4_vis_wTks.M()
    p4_BwTks = rt.TLorentzVector()
    p4_BwTks.SetPtEtaPhiM(e.BwTks_pt, e.B_eta, e.B_phi, m_B0);
    p_miss_wTks = p4_BwTks - p4_vis_wTks
    e.UmissTks = p_miss_wTks.E() - p_miss_wTks.P()

    if e.N_goodAddTks < 3:
        auxList = [e.tkCharge, e.tkPt, e.tkPtError, e.tkEta, e.tkPhi, e.massVis_wTk, e.massHad_wTk, e.massMuTk, e.mass2MissTk, e.UmissTk]
        auxList += [e.MC_tkFlag, e.MC_tkPdgId, e.MC_tkMotherPdgId, e.MC_tkMotherMotherPdgId, e.MC_tk_dphi, e.MC_tk_deta, e.MC_tk_dpt]
        for l in auxList:
            l += [0, 0, 0]

    return e

def makeSelection_JpsiKst_wrapper(kwargs):
    """
    Wrapper for the makeSelection_JpsiKst function so that we can use
    pool.map() which only allows you to pass a single argument to each
    function.
    """
    return makeSelection_JpsiKst(**kwargs)

def makeSelection_JpsiKst(tag, filenames, leafs_names, cat, start, stop, applyCorr, skipCut, serial, data):
    N_accepted_cand = []
    N_accepted_tot = 0

    tree = rt.TChain('outA/Tevts')
    for fn in filenames:
        tree.Add(fn)

    if serial:
        pb = ProgressBar(maxEntry=stop+1)
    else:
        perc = int((stop-start)*0.35)

    output = np.zeros((int(1.5*(stop-start+1)), len(leafs_names)))

    for i_ev in range(start,stop):
        tree.GetEntry(i_ev)
        ev = tree

        if serial:
            pb.show(i_ev-start)
        elif (i_ev-start) % perc == 0:
            print tag, ': {:.0f}%'.format(100*(i_ev+1-start)/(stop-start))
        N_acc = 0

        ev_output = []
        for j in range(ev.pval_piK.size()):
            evEx = extractEventInfos_JpsiKst(j, ev, applyCorr)

            if not cat is None:
                if not B02JpsiKst_selection.category_selection(j, ev, evEx, cat, True):
                    continue

            if not skipCut == 'all':
                if not B02JpsiKst_selection.candidate_selection(j, ev, evEx, skipCut):
                    continue

            N_acc += 1

            aux = (evEx.trgMu_pt, evEx.trgMu_eta, evEx.trgMu_sigdxy,
                   evEx.otherMu_pt,
                   evEx.mum_pt, evEx.mum_eta, evEx.mum_phi, ev.mum_dxy_PV[j],
                   evEx.mup_pt, evEx.mup_eta, evEx.mup_phi, ev.mup_dxy_PV[j],
                   ev.pval_mumu[j], evEx.mass_mumu,
                   evEx.Jpsi_pt, ev.cosT_Jpsi_PV[j],
                   evEx.K_pt, evEx.K_eta, evEx.K_phi, ev.K_sigdxy_PV[j],
                   evEx.pi_pt, evEx.pi_eta, evEx.pi_phi, ev.pi_sigdxy_PV[j],
                   ev.pval_piK[j], evEx.mass_piK, evEx.mass_Kpi, evEx.mass_KK,
                   evEx.Kst_pt, evEx.Kst_eta, evEx.Kst_pt, ev.sigdxy_vtxKst_PV[j],
                   ev.pval_mumupiK[j],
                   evEx.mass_mumupiK, evEx.mass_mumupiK_cJpsi, evEx.mass_mumupiK_cJpsi_cKst,
                   evEx.B_pt, evEx.B_eta, evEx.B_phi,
                   ev.cos_B_PV_mumupiK[j], ev.sigd_vtxB_PV_mumupiK[j], ev.sigdxy_vtxB_PV[j],
                   evEx.mass_mumuKpi, evEx.mass_mumuKpi_cJpsi, evEx.mass_mumuKpi_cJpsi_cKst,
                   evEx.isAntiB, evEx.mass_candKst, evEx.mass_candB,
                   B02JpsiKst_selection.category_selection(j, ev, evEx, categories['low']),
                   B02JpsiKst_selection.category_selection(j, ev, evEx, categories['mid']),
                   B02JpsiKst_selection.category_selection(j, ev, evEx, categories['high']),
                   ev.N_vertexes
                  )
            if not data:
                aux += (ev.MC_B_pt, ev.MC_B_eta,
                        ev.MC_idxCand == j,
                        ev.MC_mup_pt, ev.MC_mup_eta,
                        ev.MC_mum_pt, ev.MC_mum_eta,
                        ev.MC_d_vtxB, ev.MC_dxy_vtxB,
                        ev.d_vtxB_PV_mumupiK[j], ev.dxy_vtxB_PV[j],
                        ev.nTrueIntMC
                        )
            ev_output.append(aux)

        N_acc = len(ev_output)
        idx = 0
        if N_acc > 1:
            if data:
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

def makeSelection_JpsiK_wrapper(kwargs):
    """
    Wrapper for the makeSelection_JpsiK function so that we can use pool.map()
    which only allows you to pass a single argument to each function.
    """
    return makeSelection_JpsiK(**kwargs)

def makeSelection_JpsiK(tag, filenames, leafs_names, cat, start, stop, applyCorr, skipCut, serial, data):
    N_accepted_cand = []
    N_accepted_tot = 0

    tree = rt.TChain('outA/Tevts')
    for fn in filenames:
        tree.Add(fn)

    if serial:
        pb = ProgressBar(maxEntry=stop+1)
    else:
        perc = int((stop-start)*0.35)

    output = np.zeros((int(1.5*(stop-start+1)), len(leafs_names)))

    for i_ev in range(start,stop):
        tree.GetEntry(i_ev)
        ev = tree

        if serial:
            pb.show(i_ev-start)
        elif (i_ev-start) % perc == 0:
            print tag, ': {:.0f}%'.format(100*(i_ev+1-start)/(stop-start))
        N_acc = 0

        ev_output = []
        for j in range(ev.pval_mumuK.size()):
            evEx = extractEventInfos_JpsiK(j, ev, applyCorr)

            if not cat is None:
                if not B02JpsiK_selection.category_selection(j, ev, evEx, cat, True):
                    continue

            if not skipCut == 'all':
                if not B02JpsiK_selection.candidate_selection(j, ev, evEx, skipCut):
                    continue

            N_acc += 1

            aux = (evEx.trgMu_pt, evEx.trgMu_eta, evEx.trgMu_sigdxy,
                   evEx.mum_pt, evEx.mum_eta, evEx.mum_phi, ev.mum_dxy[j],
                   evEx.mup_pt, evEx.mup_eta, evEx.mup_phi, ev.mup_dxy[j],
                   ev.pval_mumu[j], evEx.mass_mumu,
                   evEx.Jpsi_pt, ev.cosT_Jpsi_PV[j],
                   ev.K_charge[j],
                   evEx.K_pt, evEx.K_eta, evEx.K_phi, ev.K_sigdxy_PV[j],
                   ev.pval_mumuK[j], evEx.mass_mumuK,
                   evEx.mass_mumuK_cJpsi,
                   evEx.B_pt, evEx.B_eta,
                   ev.cos_B_PV_mumuK[j], np.abs(ev.sigd_vtxB_PV_mumuK[j]),
                   ev.cosT_B_PV_mumuK[j], np.abs(ev.sigdxy_vtxB_PV_mumuK[j]),
                   B02JpsiK_selection.category_selection(j, ev, evEx, categories['low']),
                   B02JpsiK_selection.category_selection(j, ev, evEx, categories['mid']),
                   B02JpsiK_selection.category_selection(j, ev, evEx, categories['high']),
                   ev.N_vertexes
                  )
            if not data:
                aux += (ev.MC_B_pt, ev.MC_B_eta,
                        ev.MC_idxCand == j,
                        ev.MC_mup_pt, ev.MC_mup_eta,
                        ev.MC_mum_pt, ev.MC_mum_eta
                        )
            ev_output.append(aux)

        N_acc = len(ev_output)
        idx = 0
        if N_acc > 1:
            if data:
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

def makeSelection_DstMu_wrapper(kwargs):
    """
    Wrapper for the makeSelection_DstMu function so that we can use pool.map()
    which only allows you to pass a single argument to each function.
    """
    return makeSelection_DstMu(**kwargs)

def makeSelection_DstMu(tag, filenames, leafs_names, cat, start, stop, applyCorr, skipCut, trkControlRegion, serial, data):
    N_accepted_cand = []
    N_accepted_tot = 0

    tree = rt.TChain('outA/Tevts')
    for fn in filenames:
        tree.Add(fn)

    if serial:
        pb = ProgressBar(maxEntry=stop+1)
    else:
        perc = int((stop-start)*0.35)

    output = np.zeros((int(1.5*(stop-start+1)), len(leafs_names)))

    for i_ev in range(start,stop):
        tree.GetEntry(i_ev)
        ev = tree

        if serial:
            pb.show(i_ev-start)
        elif (i_ev-start) % perc == 0:
            print tag, ': {:.0f}%'.format(100*(i_ev+1-start)/(stop-start))
        N_acc = 0

        ev_output = []
        for j in range(ev.pval_piK.size()):
            idxTrg = int(ev.mu_trgMu_idx[j]) if hasattr(ev, 'mu_trgMu_idx') else int(ev.mu_trgCand_idx[j])
            evEx = extractEventInfos_DstMu(j, ev, applyCorr)

            if not cat is None:
                if not B02DstMu_selection.trigger_selection(idxTrg, ev, evEx, cat):
                    continue

            if not skipCut == 'all':
                if not B02DstMu_selection.candidate_selection(j, ev, evEx, skipCut, trkControlRegion):
                    continue

            N_acc += 1

            aux = (evEx.q2, evEx.Est_mu, evEx.M2_miss, evEx.U_miss,
                   evEx.q2_coll, evEx.Est_mu_coll, evEx.M2_miss_coll,
                   ev.mu_charge[j], evEx.mu_pt, evEx.mu_eta, evEx.mu_phi, ev.trgMu_sigdxy_BS[idxTrg],
                   ev.mu_dca_vtxDst[j], ev.mu_sigdca_vtxDst[j],
                   ev.mu_dcaT_vtxDst[j], ev.mu_sigdcaT_vtxDst[j],
                   ev.mu_dca_vtxDstMu[j], ev.mu_sigdca_vtxDstMu[j],
                   ev.mu_dcaT_vtxDstMu[j], ev.mu_sigdcaT_vtxDstMu[j],
                   evEx.B_pt, evEx.B_eta, evEx.B_phi,
                   evEx.Dst_pt, evEx.Dst_eta, evEx.Dst_phi,
                   evEx.D0_pt, evEx.D0_eta, evEx.D0_phi,
                   evEx.pi_pt, evEx.pi_eta, evEx.pi_phi,
                   evEx.K_pt, evEx.K_eta, evEx.K_phi,
                   ev.pval_piK[j], ev.sigdxy_vtxD0_PV[j],
                   evEx.pis_pt, evEx.pis_eta, evEx.pis_phi,
                   ev.pval_D0pis[j],
                   evEx.mass_piK, evEx.mass_D0pis, evEx.mass_D0pismu,
                   evEx.mass_D0pismu_muASpi, evEx.mass_D0pismu_muASK,
                   evEx.D0pismu_eta, evEx.D0pismu_phi,
                   ev.pval_D0pismu[j], ev.chi2_D0pismu[j],
                   ev.d_vtxD0pismu_PV[j], ev.dxy_vtxD0pismu_PV[j],
                   ev.cos_D0pismu_PV[j], ev.cosT_D0pismu_PV[j],
                   evEx.N_goodAddTks,
                   evEx.tkCharge[0], evEx.tkCharge[1], evEx.tkCharge[2],
                   evEx.tkPt[0], evEx.tkPt[1], evEx.tkPt[2],
                   evEx.tkPtError[0], evEx.tkPtError[1], evEx.tkPtError[2],
                   evEx.tkEta[0], evEx.tkEta[1], evEx.tkEta[2],
                   evEx.tkPhi[0], evEx.tkPhi[1], evEx.tkPhi[2],
                   evEx.massVis_wTk[0], evEx.massVis_wTk[1], evEx.massVis_wTk[2],
                   evEx.massHad_wTk[0], evEx.massHad_wTk[1], evEx.massHad_wTk[2],
                   evEx.massMuTk[0], evEx.massMuTk[1], evEx.massMuTk[2],
                   evEx.mass2MissTk[0], evEx.mass2MissTk[1], evEx.mass2MissTk[2],
                   evEx.UmissTk[0], evEx.UmissTk[1], evEx.UmissTk[2],
                   evEx.massVisTks,
                   evEx.massHadTks,
                   evEx.massHadTks_DstMassConstraint,
                   evEx.UmissTks,
                   B02DstMu_selection.trigger_selection(idxTrg, ev, evEx, categories['low']),
                   B02DstMu_selection.trigger_selection(idxTrg, ev, evEx, categories['mid']),
                   B02DstMu_selection.trigger_selection(idxTrg, ev, evEx, categories['high']),
                   ev.trgMu_HLT_Mu12_IP6[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu12_IP6') else ev.trgObj_HLT_Mu12_IP6[idxTrg],
                   ev.trgMu_HLT_Mu9_IP6[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu9_IP6') else ev.trgObj_HLT_Mu9_IP6[idxTrg],
                   ev.trgMu_HLT_Mu7_IP4[idxTrg] if hasattr(ev, 'trgMu_HLT_Mu7_IP4') else ev.trgObj_HLT_Mu7_IP4[idxTrg],
                   ev.N_vertexes, ev.localVertexDensity[j]
                  )
            if not data:
                muSisPdgId = []
                for id in ev.MC_muSistersPdgId:
                    if np.abs(id) == 14: continue #neutrino
                    muSisPdgId.append(id)
                while len(muSisPdgId) < 2:
                    muSisPdgId.append(0)

                aux += (ev.MC_q2, ev.MC_Est_mu, ev.MC_M2_miss,
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
                        evEx.MC_tkFlag[0], evEx.MC_tkFlag[1],
                        evEx.MC_tk_dpt[0], evEx.MC_tk_dpt[1],
                        evEx.MC_tk_deta[0], evEx.MC_tk_deta[1],
                        evEx.MC_tk_dphi[0], evEx.MC_tk_dphi[1],
                        evEx.MC_tkPdgId[0], evEx.MC_tkPdgId[1],
                        evEx.MC_tkMotherPdgId[0], evEx.MC_tkMotherPdgId[1],
                        evEx.MC_tkMotherMotherPdgId[0], evEx.MC_tkMotherMotherPdgId[1],
                        ev.nTrueIntMC
                       )
            # FIXME: Should we use something better here?
            basename = os.path.basename(filenames[0])
            if 'mu' in basename.lower() or 'tau' in basename.lower():
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
                       )
            ev_output.append(aux)

        N_acc = len(ev_output)
        idx = 0
        if N_acc > 1:
            if data:
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

def create_dSet_JpsiKst(filenames, outdir, cat, applyCorrections=False, skipCut=None, maxEvents=1e15, recreate=False, data=False):
    if cat is None:
        catName = 'None'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print catName
    if data:
        N_evts_per_job = 100000
    else:
        N_evts_per_job = 30000

    fskimmed_name = join(outdir,catName)

    if skipCut is not None:
        print 'Skipping cut(s)', skipCut
        if skipCut == 'all':
            fskimmed_name += '_skipall'
        else:
            fskimmed_name += '_skip'+'-'.join([str(i) for i in skipCut])
    if applyCorrections:
        print 'Applying corrections'
        fskimmed_name += '_corr'
    else:
        fskimmed_name += '_bare'
    fskimmed_name += '.root'
    logfile = fskimmed_name.replace('.root', '.log')
    if os.path.isfile(fskimmed_name) and not recreate:
        print 'Already present'
    else:
        tree = rt.TChain('outA/Tevts')
        hAllNvtx = rt.TH1D('hAllNvtx', 'hAllNvtx', 101, -0.5, 100.5)
        hAllVtxZ = rt.TH1D('hAllVtxZ', 'hAllVtxZ', 100, -25, 25)
        hAllNTrueIntMC = rt.TH1D('hAllNTrueIntMC', 'hAllNTrueIntMC', 101, -0.5, 100.5)

        print "Analyzing %i files" % len(filenames)
        for i, fn in enumerate(filenames):
            try:
                tree.Add(fn)
                fAux = rt.TFile.Open(fn, 'READ')
                hAux = fAux.Get('trgF/hAllNvts')
                hAllNvtx.Add(hAux)
                hAux = fAux.Get('trgF/hAllVtxZ')
                hAllVtxZ.Add(hAux)
                if not data:
                    hAux = fAux.Get('trgF/hAllNTrueIntMC')
                    hAllNTrueIntMC.Add(hAux)
                fAux.Close()
            except:
                print >> sys.stderr, '[ERROR] Problem with vertexes histograms in %s' % fn
                raise
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEvents, tree.GetEntries())
        print 'Total number of candidate events = %i' % N_cand_in

        leafs_names = [ 'trgMu_pt', 'trgMu_eta', 'trgMu_sigdxy',
                        'otherMu_pt',
                        'mum_pt', 'mum_eta', 'mum_phi', 'mum_dxy',
                        'mup_pt', 'mup_eta', 'mup_phi', 'mup_dxy',
                        'pval_mumu', 'mass_mumu',
                        'Jpsi_pt', 'cosT_Jpsi_PV',
                        'K_pt','K_eta','K_phi', 'K_sigdxy_PV',
                        'pi_pt','pi_eta','pi_phi', 'pi_sigdxy_PV',
                        'pval_piK', 'mass_piK', 'mass_Kpi', 'mass_KK',
                        'Kst_pt', 'Kst_eta', 'Kst_phi', 'sigdxy_vtxKst_PV',
                        'pval_mumupiK', 'mass_mumupiK',
                        'mass_mumupiK_cJpsi', 'mass_mumupiK_cJpsi_cKst',
                        'B_pt', 'B_eta', 'B_phi',
                        'cos_B_PV', 'sigd_vtxB_PV', 'sigdxy_vtxB_PV',
                        'mass_mumuKpi', 'mass_mumuKpi_cJpsi', 'mass_mumuKpi_cJpsi_cKst',
                        'isAntiB', 'mass_candKst', 'mass_candB',
                        'cat_low', 'cat_mid', 'cat_high',
                        'N_vtx'
                      ]
        if not data:
            leafs_names += ['MC_B_pt', 'MC_B_eta',
                            'MC_idxMatch',
                            'MC_mup_pt', 'MC_mup_eta',
                            'MC_mum_pt', 'MC_mum_eta',
                            'MC_d_vtxB', 'MC_dxy_vtxB',
                            'd_vtxB_PV', 'dxy_vtxB_PV',
                            'MC_nInteractions'
                           ]

        # FIXME: Should this be here?
        applyCorr = None
        if applyCorrections:
            applyCorr = 'MC'
            if data:
                applyCorr = 'RD'

        if N_cand_in < 1.5*N_evts_per_job:
            # If there are very few events per job, just run the selection
            # directly without using pool.map() or condor
            output, N_accepted_cand = makeSelection_JpsiKst(tag='', filenames=filenames, leafs_names=leafs_names, cat=cat, start=0, stop=N_cand_in, applyCorr=applyCorr, skipCut=skipCut, serial=True)
        else:
            pdiv = list(range(0, N_cand_in, N_evts_per_job))
            if not pdiv[-1] == N_cand_in:
                pdiv.append(N_cand_in)
            print 'Will be divided into ' + str(len(pdiv)-1) + ' jobs'
            inputs = []
            for i, (start, stop) in enumerate(zip(pdiv[:-1],pdiv[1:])):
                inputs.append({'tag':str(i), 'filenames':filenames, 'leafs_names':leafs_names, 'cat':cat, 'start': start, 'stop': stop, 'applyCorr': applyCorr, 'skipCut': skipCut, 'serial':False, 'data':data})
            print ' '

            start = time.time()
            if args.parallelType == 'pool' or len(inputs) < 15:
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection_JspiKst_wrapper, inputs)
            elif args.parallelType == 'jobs':
                id = uuid.uuid1().hex
                tmpDir = 'tmp/B2JpsiKst_skimCAND_' + id
                os.makedirs(tmpDir + '/out')
                tmpDir = os.path.abspath(tmpDir)
                for ii, inAux in enumerate(inputs):
                    with open(join(tmpDir,'input_%i.p' % ii), 'wb') as f:
                        pickle.dump(inAux, f)
                createSubmissionFile(tmpDir, len(inputs), 'B2JpsiKst')
                print 'Submitting jobs'
                cmd = 'condor_submit {}/jobs.jdl'.format(tmpDir)
                cmd += ' -batch-name skim_' + id
                status, output = commands.getstatusoutput(cmd)
                if status != 0:
                    print >> sys.stderr, "Error in processing command: '%s'" % cmd
                    print >> sys.stderr, "Output: %s" % output
                print 'Job submitted'
                print 'Waiting for jobs to be finished'
                time.sleep(20)
                proceed=False
                while not proceed:
                    status, output = commands.getstatusoutput('condor_q')
                    found = False
                    for line in output.split('\n'):
                        if 'skim_' + id in line:
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
            for o in outputs:
                N_accepted_cand += o[1]
            print 'Total time: {:.1f} min'.format((time.time()-start)/60.)


        dset = pd.DataFrame(output, columns=leafs_names)
        if not os.path.isdir(os.path.dirname(fskimmed_name)):
            os.makedirs(os.path.dirname(fskimmed_name))
        rtnp.array2root(dset.to_records(), fskimmed_name, treename='Tevts', mode='RECREATE')
        fAux = rt.TFile.Open(fskimmed_name, 'UPDATE')
        hAllNvtx.Write()
        hAllVtxZ.Write()
        if not data:
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

def create_dSet_JpsiK(filenames, outdir, cat, applyCorrections=False, skipCut=None, maxEvents=1e15, recreate=False, data=False):
    if cat is None:
        catName = 'None'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print catName
    if data:
        N_evts_per_job = 150000
    else:
        N_evts_per_job = 20000

    fskimmed_name = join(outdir,catName)

    if skipCut is not None:
        print 'Skipping cut(s)', skipCut
        if skipCut == 'all':
            fskimmed_name += '_skipall'
        else:
            fskimmed_name += '_skip'+'-'.join([str(i) for i in skipCut])
    if applyCorrections:
        print 'Applying corrections'
        fskimmed_name += '_corr'
    else:
        fskimmed_name += '_bare'
    fskimmed_name += '.root'
    logfile = fskimmed_name.replace('.root', '.log')
    if os.path.isfile(fskimmed_name) and not recreate:
        print 'Already present'
    else:
        tree = rt.TChain('outA/Tevts')
        hAllNvtx = rt.TH1D('hAllNvtx', 'hAllNvtx', 101, -0.5, 100.5)
        hAllVtxZ = rt.TH1D('hAllVtxZ', 'hAllVtxZ', 100, -25, 25)

        print "Analyzing %i files" % len(filenames)
        for i, fn in enumerate(filenames):
            try:
                tree.Add(fn)
                fAux = rt.TFile.Open(fn, 'READ')
                hAux = fAux.Get('trgF/hAllNvts')
                hAllNvtx.Add(hAux)
                hAux = fAux.Get('trgF/hAllVtxZ')
                hAllVtxZ.Add(hAux)
                fAux.Close()
            except:
                print >> sys.stderr, '[ERROR] Problem with vertexes histograms in %s' % fn
                raise
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEvents, tree.GetEntries())
        print 'Total number of candidate events = %i' % N_cand_in

        leafs_names = [ 'trgMu_pt', 'trgMu_eta', 'trgMu_sigdxy',
                        'mum_pt', 'mum_eta', 'mum_phi', 'mum_dxy',
                        'mup_pt', 'mup_eta', 'mup_phi', 'mup_dxy',
                        'pval_mumu', 'mass_mumu',
                        'Jpsi_pt', 'cosT_Jpsi_PV',
                        'K_charge',
                        'K_pt','K_eta','K_phi', 'K_sigdxy_PV',
                        'pval_mumuK', 'mass_mumuK',
                        'mass_mumuK_cJpsi',
                        'B_pt', 'B_eta',
                        'cos_B_PV', 'sigd_vtxB_PV',
                        'cosT_B_PV', 'sigdxy_vtxB_PV',
                        'cat_low', 'cat_mid', 'cat_high',
                        'N_vtx'
                      ]
        if not data:
            leafs_names += ['MC_B_pt', 'MC_B_eta',
                            'MC_idxMatch',
                            'MC_mup_pt', 'MC_mup_eta',
                            'MC_mum_pt', 'MC_mum_eta'
                           ]

        applyCorr = None
        if applyCorrections:
            applyCorr = 'MC'
            if data:
                applyCorr = 'RD'

        if N_cand_in < 1.5*N_evts_per_job:
            output, N_accepted_cand = makeSelection_JpsiK(tag='', filenames=filenames, leafs_names=leafs_names, cat=cat, start=0, stop=N_cand_in, applyCorr=applyCorr, skipCut=skipCut, serial=True)
        else:
            pdiv = list(range(0, N_cand_in, N_evts_per_job))
            if not pdiv[-1] == N_cand_in:
                pdiv.append(N_cand_in)
            print 'Will be divided into ' + str(len(pdiv)-1) + ' jobs'
            inputs = []
            for i, (start, stop) in enumerate(zip(pdiv[:-1],pdiv[1:])):
                inputs.append({'tag':str(i), 'filenames':filenames, 'leafs_names':leafs_names, 'cat':cat, 'start': start, 'stop': stop, 'applyCorr': applyCorr, 'skipCut': skipCut, 'serial':False, 'data':data})
            print ' '

            start = time.time()
            if args.parallelType == 'pool' or len(inputs) < 15:
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection_JpsiK_wrapper, inputs)
            elif args.parallelType == 'jobs':
                id = uuid.uuid1().hex
                tmpDir = 'tmp/B2JpsiK_skimCAND_' + id
                os.makedirs(tmpDir + '/out')
                tmpDir = os.path.abspath(tmpDir)
                for ii, inAux in enumerate(inputs):
                    with open(join(tmpDir,'input_%i.p' % ii), 'wb') as f:
                        pickle.dump(inAux, f)
                createSubmissionFile(tmpDir, len(inputs), 'B2JpsiK')
                print 'Submitting jobs'
                cmd = 'condor_submit {}/jobs.jdl'.format(tmpDir)
                cmd += ' -batch-name skim_' + id
                status, output = commands.getstatusoutput(cmd)
                if status != 0:
                    print >> sys.stderr, "Error in processing command: '%s'" % cmd
                    print >> sys.stderr, "Output: %s" % output
                print 'Job submitted'
                print 'Waiting for jobs to be finished'
                time.sleep(20)
                proceed=False
                while not proceed:
                    status, output = commands.getstatusoutput('condor_q')
                    found = False
                    for line in output.split('\n'):
                        if 'skim_' + id in line:
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
            for o in outputs:
                N_accepted_cand += o[1]
            print 'Total time: {:.1f} min'.format((time.time()-start)/60.)


        dset = pd.DataFrame(output, columns=leafs_names)
        if not os.path.isdir(os.path.dirname(fskimmed_name)):
            os.makedirs(os.path.dirname(fskimmed_name))
        rtnp.array2root(dset.to_records(), fskimmed_name, treename='Tevts', mode='RECREATE')
        fAux = rt.TFile.Open(fskimmed_name, 'UPDATE')
        hAllNvtx.Write()
        hAllVtxZ.Write()
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

def create_dSet_DstMu(filenames, outdir, cat, applyCorrections=False, skipCut=None, trkControlRegion=False, maxEvents=1e15, recreate=False, data=False):
    if cat is None:
        catName = 'NoCat'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print catName
    if data:
        N_evts_per_job = 100000
    else:
        N_evts_per_job = 30000

    fskimmed_name = join(outdir,catName)

    if skipCut is not None:
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
        print 'Applying corrections'
        fskimmed_name += '_corr'
    else:
        fskimmed_name += '_bare'
    fskimmed_name += '.root'
    logfile = fskimmed_name.replace('.root', '.log')
    if os.path.isfile(fskimmed_name) and not recreate:
        print 'Already present'
    else:
        tree = rt.TChain('outA/Tevts')
        hAllNvtx = rt.TH1D('hAllNvtx', 'hAllNvtx', 101, -0.5, 100.5)
        hAllVtxZ = rt.TH1D('hAllVtxZ', 'hAllVtxZ', 100, -25, 25)
        hAllNTrueIntMC = rt.TH1D('hAllNTrueIntMC', 'hAllNTrueIntMC', 101, -0.5, 100.5)

        print "Analyzing %i files" % len(filenames)
        for i, fn in enumerate(filenames):
            try:
                tree.Add(fn)
                fAux = rt.TFile.Open(fn, 'READ')
                hAux = fAux.Get('trgF/hAllNvts')
                hAllNvtx.Add(hAux)
                hAux = fAux.Get('trgF/hAllVtxZ')
                hAllVtxZ.Add(hAux)
                if not data:
                    hAux = fAux.Get('trgF/hAllNTrueIntMC')
                    hAllNTrueIntMC.Add(hAux)
                fAux.Close()
            except:
                print >> sys.stderr, '[ERROR] Problem with vertexes histograms in %s' % fn
                raise
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEvents, tree.GetEntries())
        print 'Total number of candidate events = %i' % N_cand_in

        leafs_names = ['q2', 'Est_mu', 'M2_miss', 'U_miss',
                       'q2_coll', 'Est_mu_coll', 'M2_miss_coll',
                       'mu_charge', 'mu_pt', 'mu_eta', 'mu_phi', 'mu_sigdxy',
                       'mu_dca_vtxDst', 'mu_sigdca_vtxDst',
                       'mu_dcaT_vtxDst', 'mu_sigdcaT_vtxDst',
                       'mu_dca_vtxDstMu', 'mu_sigdca_vtxDstMu',
                       'mu_dcaT_vtxDstMu', 'mu_sigdcaT_vtxDstMu',
                       'B_pt', 'B_eta', 'B_phi',
                       'Dst_pt', 'Dst_eta', 'Dst_phi',
                       'D0_pt', 'D0_eta', 'D0_phi',
                       'pi_pt', 'pi_eta', 'pi_phi',
                       'K_pt', 'K_eta', 'K_phi',
                       'pval_piK', 'sigdxy_vtxD0_PV',
                       'pis_pt', 'pis_eta', 'pis_phi',
                       'pval_D0pis',
                       'mass_piK', 'mass_D0pis', 'mass_D0pismu',
                       'mass_D0pismu_muASpi', 'mass_D0pismu_muASK',
                       'D0pismu_eta', 'D0pismu_phi',
                       'pval_D0pismu', 'chi2_D0pismu',
                       'd_vtxD0pismu_PV', 'dxy_vtxD0pismu_PV',
                       'cos_D0pismu_PV', 'cosT_D0pismu_PV',
                       'N_goodAddTks',
                       'tkCharge_0', 'tkCharge_1', 'tkCharge_2',
                       'tkPt_0', 'tkPt_1', 'tkPt_2',
                       'tkPtError_0', 'tkPtError_1', 'tkPtError_2',
                       'tkEta_0', 'tkEta_1', 'tkEta_2',
                       'tkPhi_0', 'tkPhi_1', 'tkPhi_2',
                       'tkMassVis_0', 'tkMassVis_1', 'tkMassVis_2',
                       'tkMassHad_0', 'tkMassHad_1', 'tkMassHad_2',
                       'tkMassMuTk_0', 'tkMassMuTk_1', 'tkMassMuTk_2',
                       'tkMassMiss2_0', 'tkMassMiss2_1', 'tkMassMiss2_2',
                       'tkUmiss_0', 'tkUmiss_1', 'tkUmiss_2',
                       'massVisTks',
                       'massHadTks',
                       'massHadTks_DstMassConstraint',
                       'UmissTks',
                       'cat_low', 'cat_mid', 'cat_high',
                       'muPass_Mu12_IP6', 'muPass_Mu9_IP6', 'muPass_Mu7_IP4',
                       'N_vtx', 'localVertexDensity'
                      ]
        if not data:
            leafs_names += ['MC_q2', 'MC_Est_mu', 'MC_M2_miss',
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
                            'MC_tkFlag_0', 'MC_tkFlag_1',
                            'MC_tk_dpt_0', 'MC_tk_dpt_1',
                            'MC_tk_deta_0', 'MC_tk_deta_1',
                            'MC_tk_dphi_0', 'MC_tk_dphi_1',
                            'MC_tkPdgId_0', 'MC_tkPdgId_1',
                            'MC_tkMotherPdgId_0', 'MC_tkMotherPdgId_1',
                            'MC_tkMotherMotherPdgId_0', 'MC_tkMotherMotherPdgId_1',
                            'MC_nInteractions'
                           ]
        # FIXME: Should we use something better here?
        basename = os.path.basename(filenames[0])
        if 'mu' in basename.lower() or 'tau' in basename.lower():
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
                            ]

        applyCorr = None
        if applyCorrections:
            applyCorr = 'MC'
            if data:
                applyCorr = 'RD'

        if N_cand_in < 1.5*N_evts_per_job:
            output, N_accepted_cand = makeSelection_DstMu(tag='', filenames=filenames, leafs_names=leafs_names, cat=cat, start=0, stop=N_cand_in, applyCorr=applyCorr, skipCut=skipCut, serial=True, trkControlRegion=trkControlRegion, data=data)
        else:
            pdiv = list(range(0, N_cand_in, N_evts_per_job))
            if not pdiv[-1] == N_cand_in:
                pdiv.append(N_cand_in)
            print 'Will be divided into ' + str(len(pdiv)-1) + ' jobs'
            inputs = []
            for i, (start, stop) in enumerate(zip(pdiv[:-1],pdiv[1:])):
                inputs.append({'tag':str(i), 'filenames':filenames, 'leafs_names':leafs_names, 'cat':cat, 'start': start, 'stop': stop, 'applyCorr': applyCorr, 'skipCut': skipCut, 'serial':False, 'trkControlRegion': trkControlRegion, 'data':data})
            print ' '

            start = time.time()
            if args.parallelType == 'pool' or len(inputs) < 15:
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection_DstMu_wrapper, inputs)
            elif args.parallelType == 'jobs':
                id = uuid.uuid1().hex
                tmpDir = 'tmp/B2DstMu_skimCAND_' + id
                if args.trkControlRegion:
                    tmpDir += '_trkControl'
                os.makedirs(tmpDir + '/out')
                tmpDir = os.path.abspath(tmpDir)
                for ii, inAux in enumerate(inputs):
                    with open(join(tmpDir,'input_%i.p' % ii), 'wb') as f:
                        pickle.dump(inAux, f)
                createSubmissionFile(tmpDir, len(inputs), 'B2DstMu')
                print 'Submitting jobs'
                cmd = 'condor_submit %s' % join(tmpDir,'jobs.jdl')
                cmd += ' -batch-name skim_' + id
                if args.trkControlRegion:
                    cmd += '_trkControl'
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
                        aux = 'skim_' + id + ('_trkControl' if args.trkControlRegion else '')
                        if aux in line:
                            print line
                            time.sleep(10)
                            found = True
                    proceed = not found

                outputs = []
                for ii in range(len(inputs)):
                    with open(join(tmpDir,'output_%i.p' % ii), 'rb') as f:
                        o = pickle.load(f)
                        outputs.append(o)

            output = np.concatenate(tuple([o[0] for o in outputs]))
            N_accepted_cand = []
            for o in outputs:
                N_accepted_cand += o[1]
            print 'Total time: {:.1f} min'.format((time.time()-start)/60.)


        dset = pd.DataFrame(output, columns=leafs_names)
        if not os.path.isdir(os.path.dirname(fskimmed_name)):
            os.makedirs(os.path.dirname(fskimmed_name))
        rtnp.array2root(dset.to_records(), fskimmed_name, treename='Tevts', mode='RECREATE')
        fAux = rt.TFile.Open(fskimmed_name, 'UPDATE')
        hAllNvtx.Write()
        hAllVtxZ.Write()
        if not data:
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

def createSubmissionFile(tmpDir, njobs, selection):
    job_file = join(tmpDir,'job.sh')
    with open(job_file, 'w') as fjob:
        fjob.write('#!/bin/bash\n')
        fjob.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        fjob.write('cd %s/RDstAnalysis/CMSSW_10_2_3/\n' % os.environ['HOME'])
        fjob.write('eval `scramv1 runtime -sh`\n')
        fjob.write('cd %s/RDstAnalysis/BPH_RD_Analysis/\n' % os.environ['HOME'])
        fjob.write('export PYTHONPATH=%s/RDstAnalysis/BPH_RD_Analysis/lib:$PYTHONPATH\n' % os.environ['HOME'])
        fjob.write('export PYTHONPATH=%s/RDstAnalysis/BPH_RD_Analysis/analysis:$PYTHONPATH\n' % os.environ['HOME'])
        fjob.write('python ./scripts/skimmer.py --selection %s --function makeSel --tmpDir $1 --jN $2\n' % selection)
        os.system('chmod +x {}/job.sh'.format(tmpDir))

    sub_file = join(tmpDir,'jobs.jdl')
    with open(sub_file, 'w') as fsub:
        fsub.write('executable    = %s\n' % join(tmpDir,'job.sh'))
        fsub.write('arguments     = {} $(ProcId)\n'.format(tmpDir))
        fsub.write('output        = {}/out/job_$(ProcId)_$(ClusterId).out\n'.format(tmpDir))
        fsub.write('error         = {}/out/job_$(ProcId)_$(ClusterId).err\n'.format(tmpDir))
        fsub.write('log           = {}/out/job_$(ProcId)_$(ClusterId).log\n'.format(tmpDir))
        fsub.write('WHEN_TO_TRANSFER_OUTPUT = ON_EXIT_OR_EVICT\n')
        fsub.write('+MaxRuntime   = 3600\n')
        fsub.write('+RunAsOwner = True\n')
        fsub.write('+InteractiveUser = True\n')
        fsub.write('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7-m20200605"\n')
        fsub.write('+SingularityBindCVMFS = True\n')
        fsub.write('run_as_owner = True\n')
        fsub.write('RequestDisk = 2000000\n')
        fsub.write('RequestMemory = 2500\n')
        fsub.write('RequestCpus = 1\n')
        fsub.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        fsub.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        # Send the job to Held state on failure.
        fsub.write('on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')
        # Periodically retry the jobs for 2 times with an interval of 20 minutes.
        fsub.write('periodic_release =  (NumJobStarts < 2) && ((CurrentTime - EnteredCurrentStatus) > (60*20))\n')
        fsub.write('+PeriodicRemove = ((JobStatus =?= 2) && ((MemoryUsage =!= UNDEFINED && MemoryUsage > 2.5*RequestMemory)))\n')
        fsub.write('max_retries    = 3\n')
        fsub.write('requirements   = Machine =!= LastRemoteHost && regexp("blade-.*", TARGET.Machine)\n')
        fsub.write('universe = vanilla\n')
        fsub.write('queue %i\n' % njobs)
        fsub.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument ('filenames', nargs='*', help='Files to skim')
    parser.add_argument ('-o', '--outdir', default=None, required=False, help='Output file directory')
    parser.add_argument ('--data', default=False, action='store_true', help='Files are real data, not MC')
    parser.add_argument ('--function', type=str, default='main', help='Function to perform')
    parser.add_argument ('-p', '--parallelType', choices=['pool', 'jobs'], default='jobs', help='Function to perform')
    parser.add_argument ('--maxEvents', type=int, default=1e15, help='Max number of events to be processed')
    parser.add_argument ('--recreate', default=False, action='store_true', help='Recreate even if file already present')
    parser.add_argument ('--applyCorr', default=False, action='store_true', help='Switch to apply crrections')
    parser.add_argument ('--trkControlRegion', default=False, action='store_true', help='Track control region selection')
    parser.add_argument ('--cat', type=str, default=['low', 'mid', 'high'], choices=['single', 'low', 'mid', 'high', 'none'], help='Category(ies)', nargs='+')
    parser.add_argument ('--skipCut', type=str, default='', choices=['all', '11', '13', '14', '16', '17'], help='Cut to skip.\nAll: skip all the cuts\n16:Visible mass cut\n17: additional tracks cut')
    parser.add_argument ('--selection', type=str, default=None, choices=['B2DstMu', 'B2JpsiK', 'B2JpsiKst'], help='Selection ("B2DstMu", "B2JpsiK", or "B2JpsiKst")', required=True)
    ######## Arguments not for user #####################
    parser.add_argument ('--tmpDir', type=str, default=None, help='Temporary directory')
    parser.add_argument ('--jN', type=int, default=None, help='Job number')
    args = parser.parse_args()

    if args.function == 'main':
        if not args.outdir:
            parser.print_help()
            print >> sys.stderr, "skimmer.py: error: argument -o/--outdir is required"
            sys.exit(1)

        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

        print '-'*50 + '\n'

        skip = []
        if args.skipCut == 'all':
            skip.append('all')
        elif args.skipCut:
            skip.append([int(args.skipCut)])
        else:
            skip.append([])

        for idx in skip:
            for cn in args.cat:
                if args.selection == 'B2DstMu':
                    create_dSet_DstMu(args.filenames, args.outdir, categories[cn], skipCut=idx, applyCorrections=args.applyCorr, trkControlRegion=args.trkControlRegion, maxEvents=args.maxEvents, recreate=args.recreate, data=args.data)
                elif args.selection == 'B2JpsiK':
                    create_dSet_JpsiK(args.filenames, args.outdir, categories[cn], skipCut=idx, applyCorrections=args.applyCorr, maxEvents=args.maxEvents, recreate=args.recreate, data=args.data)
                elif args.selection == 'B2JpsiKst':
                    create_dSet_JpsiKst(args.filenames, args.outdir, categories[cn], skipCut=idx, applyCorrections=args.applyCorr, maxEvents=args.maxEvents, recreate=args.recreate, data=args.data)
    elif args.function == 'makeSel':
        with open(join(args.tmpDir,'input_%i.p' % args.jN), 'rb') as f:
            input = pickle.load(f)
            if args.selection == 'B2DstMu':
                output = makeSelection_DstMu(**input)
            elif args.selection == 'B2JpsiK':
                output = makeSelection_JpsiK(**input)
            elif args.selection == 'B2JpsiKst':
                output = makeSelection_JpsiKst(**input)
            else:
                print >> sys.stderr, "unknown selectoin '%s'" % args.selection

        with open(join(args.tmpDir,'output_%i.p' % args.jN), 'wb') as f:
            pickle.dump(output, f)
    else:
        print >> sys.stderr, "function '%s' not recognized" % args.function
