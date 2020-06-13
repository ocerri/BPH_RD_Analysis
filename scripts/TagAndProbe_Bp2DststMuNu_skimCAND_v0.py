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

import argparse
parser = argparse.ArgumentParser()
#Example: python B2JpsiKst_skimCAND_v1.py -d n_PU20 --maxEvents 80000
parser.add_argument ('--function', type=str, default='main', help='Function to perform')
parser.add_argument ('-d', '--dataset', type=str, default=[], help='Dataset(s) to run on or regular expression for them', nargs='+')
parser.add_argument ('-p', '--parallelType', choices=['pool', 'jobs'], default='jobs', help='Function to perform')
parser.add_argument ('--maxEvents', type=int, default=1e15, help='Max number of events to be processed')
parser.add_argument ('--recreate', default=False, action='store_true', help='Recreate even if file already present')
# parser.add_argument ('--cat', type=str, default=['low', 'mid', 'high'], choices=['single', 'low', 'mid', 'high', 'probe', 'none'], help='Category(ies)', nargs='+')
######## Arguments not for user #####################
parser.add_argument ('--tmpDir', type=str, default=None, help='Temporary directory')
parser.add_argument ('--jN', type=int, default=None, help='Job number')
args = parser.parse_args()

#############################################################################
####                          Datset declaration                         ####
#############################################################################
MCloc = '../data/cmsMC_private/'
MCend = '/ntuples_TagAndProbe_Bp_MuNuDstst/out_CAND_*.root'
RDloc = '../data/cmsRD/ParkingBPH*/'

filesLocMap = {
'n_PUc0'        : MCloc+'BP_Tag_Bp_MuNuDstst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3'+MCend,
#
#
'data' : RDloc+'Run2018D-05May2019promptD-v1_RDntuplizer_TagAndProbe_Bp2MuNuDstst_Pip_200522_CAND.root'
}

def makeSelection(inputs):
    n, tag, filepath, leafs_names, cat, idxInt, serial = inputs
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
            if not ev.K_pt[j] > 0.8: continue
            if not ev.pi_pt[j] > 0.8: continue
            if not ev.pval_piK[j] > 0.1: continue
            if not np.abs(ev.K_eta[j]) < 2.4: continue
            if not np.abs(ev.pi_eta[j]) < 2.4: continue
            if not np.abs(ev.mass_piK[j] - 1.86483) < 0.035: continue
            if not ev.sigdxy_vtxD0_PV[j] > 7: continue
            if not ev.cosT_D0_PV[j] > 0.9: continue
            if not ev.pip_pt[j] > 0.5: continue
            if not np.abs(ev.pip_eta[j]) < 2.4: continue
            if not ev.mass_D0pipmu[j] < 5.27963: continue
            if not ev.pval_D0pipmu[j] > 0.1: continue
            if not ev.cos_D0pipmu_PV[j] > 0.99: continue

            N_goodPi = 0
            jPis_good = -1
            for jPis in range(ev.pis_pt.size()):
                if not 'data' in n:
                    if not ev.pis_TagCandIdx[jPis] == j: continue
                if not ev.pis_pt[jPis] > 0.4: continue
                if not np.abs(ev.pis_eta[jPis]) < 2.4: continue
                if not ev.pval_D0pis[jPis] > 0.1: continue
                if not np.abs(ev.dm_D0pis_piK[jPis] - 0.14543) < 6e-3: continue
                N_goodPi += 1
                jPis_good = jPis

            N_acc += 1

            aux = (
                   ev.D0_refitD0pipmu_pt[j], ev.D0_refitD0pipmu_eta[j], ev.D0_refitD0pipmu_phi[j],
                   ev.pip_refitD0pipmu_pt[j], ev.pip_refitD0pipmu_eta[j], ev.pip_refitD0pipmu_phi[j],
                   ev.mu_refitD0pipmu_pt[j], ev.mu_refitD0pipmu_eta[j], ev.mu_refitD0pipmu_phi[j],
                   ev.Dstst_expD0pip_pt[j], ev.Dstst_expD0pip_eta[j], ev.Dstst_expD0pip_phi[j],
                   ev.mass_piK[j], ev.mass_D0pip[j], ev.mass_expDstpip[j], ev.dm_expDstpip_pik[j], ev.mass_D0pipmu[j],
                   ev.pis_expD0pipmu_pt[j], ev.pis_expD0pipmu_eta[j], ev.pis_expD0pipmu_phi[j], 0.3*3.8/ev.pis_expD0pipmu_pt[j],
                   ev.N_vertexes
                  )
            if N_goodPi == 0:
                aux += (-1, -999, -999, -1, -999)
            else:
                aux += (ev.pis_pt[jPis_good], ev.pis_eta[jPis_good],ev.pis_phi[jPis_good], ev.sigdxy_pis_PV[jPis_good],
                        ev.dm_D0pis_piK[jPis_good]
                )
            ev_output.append(aux)

        N_acc = len(ev_output)
        idx = 0
        if N_acc > 1:
            idx = np.random.randint(len(ev_output))
        if N_acc > 0:
            output[N_accepted_tot] = ev_output[idx]
            N_accepted_tot += 1
            N_accepted_cand.append(N_acc)

    output = output[:N_accepted_tot]
    if not serial:
        print tag, ': done'
    return [output, N_accepted_cand]

def create_dSet(n, filepath, cat, maxEvents=args.maxEvents):
    if cat is None:
        catName = 'All'
    else:
        catName = cat.name
    print '\n' + 50*'-'
    print n, catName
    if 'data' in n:
        loc = '../data/cmsRD/skimmed/TnP_Bp2DststMuNu'+ n.replace('data', '')
        out = re.search('20[01][1-9][0-3][0-9]', filepath)
        if out is None:
            print filepath
            raise
        fskimmed_name = loc + '_' + out.group(0) + '_' + catName
        N_evts_per_job = 100000
    else:
        d = os.path.dirname(filepath) + '/skimmed/'
        if not os.path.isdir(d):
            os.makedirs(d)
        fskimmed_name = d + catName
        N_evts_per_job = 20000
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

        leafs_names = [ 'D0_pt', 'D0_eta', 'D0_phi',
                        'pip_pt', 'pip_eta', 'pip_phi',
                        'mu_pt', 'mu_eta', 'mu_phi',
                        'expDstst_pt', 'expDstst_eta', 'expDstst_phi',
                        'mass_piK', 'mass_D0pip', 'mass_expDstpip', 'dm_expDstpip_pik', 'mass_D0pipmu',
                        'expPis_pt', 'expPis_eta', 'expPis_phi', 'expPis_k',
                        'N_vtx',
                        'pis_pt', 'pis_eta','pis_phi', 'sigdxy_pis_PV',
                        'dm_D0pis_piK'
                      ]

        if N_cand_in < 1.5*N_evts_per_job:
            output, N_accepted_cand = makeSelection([n, '', filepath, leafs_names, cat,
                                                     [0, N_cand_in-1], True])
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
                inputs.append([n, str(i), filepath, leafs_names, cat, [pdiv[i-1]+1+corr, pdiv[i]], False])
            print ' '

            start = time.time()
            if args.parallelType == 'pool' or len(inputs) < 15:
                p = Pool(min(20,len(inputs)))
                outputs = p.map(makeSelection, inputs)
            elif args.parallelType == 'jobs':
                tmpDir = 'tmp/TnP_Bp2DststMuNu_skimCAND_' + n
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
    fjob.write('python TagAndProbe_Bp2DststMuNu_skimCAND_v0.py --function makeSel --tmpDir $1 --jN $2\n')
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

        # if 'none' in args.cat:
        print 'Running w/o category (ignoring other categories)'
        for n, fp in file_loc.iteritems():
                create_dSet(n, fp, cat=None)

    elif args.function == 'makeSel':
        tmpDir = args.tmpDir
        input = pickle.load( open( tmpDir+'/input_{}.p'.format(args.jN), 'rb' ) )
        output = makeSelection(input)
        pickle.dump(output, open( tmpDir+'/output_{}.p'.format(args.jN), 'wb' ) )

    else:
        print args.function, 'not recognized'
