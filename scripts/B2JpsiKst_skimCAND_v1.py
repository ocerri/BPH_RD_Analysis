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
from B02JpsiKst_selection import candidate_selection, category_selection

import argparse
parser = argparse.ArgumentParser()
#Example: python B2JpsiKst_skimCAND_v1.py -d n_PU20 --maxEvents 80000 --applyCorr
parser.add_argument ('--function', type=str, default='main', help='Function to perform')
parser.add_argument ('-d', '--dataset', type=str, default=[], help='Dataset(s) to run on or regular expression for them', nargs='+')
parser.add_argument ('-p', '--parallelType', choices=['pool', 'jobs'], default='jobs', help='Function to perform')
parser.add_argument ('--maxEvents', type=int, default=1e15, help='Max number of events to be processed')
parser.add_argument ('--recreate', default=False, action='store_true', help='Recreate even if file already present')
parser.add_argument ('--applyCorr', default=False, action='store_true', help='Switch to apply crrections')
parser.add_argument ('--cat', type=str, default=['low', 'mid', 'high'], choices=['single', 'low', 'mid', 'high', 'probe', 'none'], help='Category(ies)', nargs='+')
parser.add_argument ('--skipCut', type=str, default='', choices=['all', '16', '17'], help='Cut to skip')
######## Arguments not for user #####################
parser.add_argument ('--tmpDir', type=str, default=None, help='Temporary directory')
parser.add_argument ('--jN', type=int, default=None, help='Job number')
args = parser.parse_args()

#############################################################################
####                          Datset declaration                         ####
#############################################################################
MCloc = '../data/cmsMC_private/'
MCend = '/ntuples_B2JpsiKst/out_CAND_*.root'
RDloc = '../data/cmsRD/ParkingBPH*/'

filesLocMap = {
'n_PU0'         : MCloc+'BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV_PU0_10-2-3'+MCend,
'n_PU20'        : MCloc+'BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV_PU20_10-2-3'+MCend,
'n_PUc0'        : MCloc+'BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3'+MCend,
'n_PU35'        : MCloc+'BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV_PU35_10-2-3'+MCend,
#
'fsr_PU20'      : MCloc+'BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgenFSR_SVV_PU20_10-2-3'+MCend,
#
's_PU0'         : MCloc+'BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_SVV_PU0_10-2-3'+MCend,
#
#

# 'data' : RDloc+'*2018*B2JpsiKst_200521_CAND.root'
'data' : RDloc+'*2018*B2JpsiKst_200622_CAND.root'
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
        # list = list[:il] + [el] + list[il:]
        return il, list[:il] + [el] + list[il:]

def extractEventInfos(j, ev, corr=None):
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
    e.mass_Kpi = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_K, m_pi)
    # mass_Kpi_recomp = compMass(ev.piRefit_pt[j], ev.KRefit_pt[j], ev.piRefit_eta[j], ev.KRefit_eta[j], ev.piRefit_phi[j], ev.KRefit_phi[j], m_K, m_pi)
    e.mass_KK = compMass(e.pi_pt, e.K_pt, e.pi_eta, e.K_eta, e.pi_phi, e.K_phi, m_K, m_K)

    p4_B = p4_mup + p4_mum + p4_K + p4_pi
    e.mass_mumupiK = p4_B.M()
    e.B_pt = p4_B.Pt()
    e.B_eta = p4_B.Eta()
    e.B_phi = p4_B.Phi()

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
            evEx = extractEventInfos(j, ev, corr)

            if not cat is None:
                if not category_selection(j, ev, evEx, cat, True):
                    continue

            if not skipCut == 'all':
                if not candidate_selection(j, ev, evEx, skipCut):
                    continue

            N_acc += 1

            aux = (evEx.trgMu_pt, evEx.trgMu_eta, evEx.trgMu_sigdxy,
                   evEx.mum_pt, evEx.mum_eta, evEx.mum_phi, ev.mum_dxy[j],
                   evEx.mup_pt, evEx.mup_eta, evEx.mup_phi, ev.mup_dxy[j],
                   ev.pval_mumu[j], evEx.mass_mumu,
                   evEx.Jpsi_pt, ev.cosT_Jpsi_PV[j],
                   evEx.K_pt, evEx.K_eta, evEx.K_phi, ev.K_sigdxy_PV[j],
                   evEx.pi_pt, evEx.pi_eta, evEx.pi_phi, ev.pi_sigdxy_PV[j],
                   ev.pval_piK[j], evEx.mass_piK, evEx.mass_Kpi, evEx.mass_KK,
                   ev.sigdxy_vtxKst_PV[j],
                   ev.pval_mumupiK[j], evEx.mass_mumupiK,
                   evEx.B_pt, evEx.B_eta,
                   ev.cos_B_PV_mumupiK[j], ev.sigd_vtxB_PV_mumupiK[j],
                   # ev.cosT_B_PV_mumupiK[j], ev.sigdxy_vtxB_PV[j],
                   category_selection(j, ev, evEx, categories['low']),
                   category_selection(j, ev, evEx, categories['mid']),
                   category_selection(j, ev, evEx, categories['high']),
                   ev.N_vertexes
                  )
            if not 'data' in n:
                aux += (ev.MC_B_pt, ev.MC_B_eta,
                        ev.MC_idxCand == j,
                        ev.MC_mup_pt, ev.MC_mup_eta,
                        ev.MC_mum_pt, ev.MC_mum_eta
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

def create_dSet(n, filepath, cat, applyCorrections=False, skipCut=[], maxEvents=args.maxEvents):
    if cat is None:
        catName = 'None'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print n, catName
    if 'data' in n:
        loc = '../data/cmsRD/skimmed/B2JpsiKst'+ n.replace('data', '')
        out = re.search('20[01][1-9][0-3][0-9]', filepath)
        if out is None:
            print filepath
            raise
        fskimmed_name = loc + '_' + out.group(0) + '_' + catName
        N_evts_per_job = 150000
    else:
        d = os.path.dirname(filepath) + '/skimmed/'
        if not os.path.isdir(d):
            os.makedirs(d)
        fskimmed_name = d + catName
        N_evts_per_job = 20000
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
        hAllNvtx = rt.TH1D('hAllNvtx', 'hAllNvtx', 101, -0.5, 100.5)
        hAllVtxZ = rt.TH1D('hAllVtxZ', 'hAllVtxZ', 100, -25, 25)
        for fn in glob(filepath):
            tree.Add(fn)
            fAux = rt.TFile.Open(fn, 'READ')
            hAux = fAux.Get('trgF/hAllNvts')
            hAllNvtx.Add(hAux)
            hAux = fAux.Get('trgF/hAllVtxZ')
            hAllVtxZ.Add(hAux)
            fAux.Close()
        print 'Computing events from {} files'.format(tree.GetNtrees())
        N_cand_in = min(maxEvents, tree.GetEntries())
        print n, ': Total number of candidate events =', N_cand_in

        leafs_names = [ 'trgMu_pt', 'trgMu_eta', 'trgMu_sigdxy',
                        'mum_pt', 'mum_eta', 'mum_phi', 'mum_dxy',
                        'mup_pt', 'mup_eta', 'mup_phi', 'mup_dxy',
                        'pval_mumu', 'mass_mumu',
                        'Jpsi_pt', 'cosT_Jpsi_PV',
                        'K_pt','K_eta','K_phi', 'K_sigdxy_PV',
                        'pi_pt','pi_eta','pi_phi', 'pi_sigdxy_PV',
                        'pval_piK', 'mass_piK', 'mass_piK_CPconj', 'mass_KK',
                        'sigdxy_vtxKst_PV',
                        'pval_mumupiK', 'mass_mumupiK',
                        'B_pt', 'B_eta',
                        'cos_B_PV', 'sigd_vtxB_PV',
                        # 'cosT_B_PV', 'sigdxy_vtxB_PV',
                        'cat_low', 'cat_mid', 'cat_high',
                        'N_vtx'
                      ]
        if not 'data' in n:
            leafs_names += ['MC_B_pt', 'MC_B_eta',
                            'MC_idxMatch',
                            'MC_mup_pt', 'MC_mup_eta',
                            'MC_mum_pt', 'MC_mum_eta'
                           ]

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
            if args.parallelType == 'pool' or len(inputs) < 15:
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection, inputs)
            elif args.parallelType == 'jobs':
                tmpDir = 'tmp/B2JspiKst_skimCAND_' + n
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

def createSubmissionFile(tmpDir, njobs):
    fjob = open(tmpDir+'/job.sh', 'w')
    fjob.write('#!/bin/bash\n')
    fjob.write('source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /storage/user/ocerri/CMSSW_10_2_3/; eval `scramv1 runtime -sh`\n')
    fjob.write('cd /storage/user/ocerri/BPhysics/scripts\n')
    fjob.write('python B2JpsiKst_skimCAND_v1.py --function makeSel --tmpDir $1 --jN $2\n')
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
            for kn in filesLocMap.keys():
                if not re.match(n, kn) is None:
                    print 'Adding', kn
                    file_loc[kn] = filesLocMap[kn]
        if len(file_loc.keys()) == 0:
            print 'No dataset provided'
            exit()

        recreate = []
        if args.recreate:
            recreate = file_loc.keys()
        print '-'*50 + '\n'

        if 'none' in args.cat:
            print 'Running w/o category (ignoring other categories)'
            sc = []
            if args.skipCut == 'all':
                sc = 'all'
            for n, fp in file_loc.iteritems():
                    create_dSet(n, fp, cat=None, skipCut=sc, applyCorrections=args.applyCorr)

        else:
            skip = []
            if args.skipCut == 'all':
                skip.append('all')
            # skip.append([6, 11, 12]) #Mass D0, D* and D*-D0 (m piK)
            # skip.append([16]) #Visible mass (m D0pismu)
            # skip.append([17]) #Additional tracks
            elif args.skipCut:
                skip.append([int(args.skipCut)])
            else:
                skip.append([])

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
