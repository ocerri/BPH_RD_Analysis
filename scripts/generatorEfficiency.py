#!/usr/bin/env python
import sys, os, re, yaml, pickle
import commands
from glob import glob
from prettytable import PrettyTable
sys.path.append('../lib')

import time, datetime

import signal

import numpy as np
from scipy.stats import mode
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from progressBar import ProgressBar

import uproot as ur
import ROOT as rt
rt.gErrorIgnoreLevel = rt.kError
rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

# load FWLite C++ libraries
rt.gSystem.Load("libFWCoreFWLite.so");
rt.gSystem.Load("libDataFormatsFWLite.so");
rt.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Lumis
from DataFormats.FWLite import Handle

from analysis_utilities import DSetLoader, NTUPLE_TAG

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

order = ['mu_c0', 'tau_c0', 'DstPip_c0', 'DstPi0_c0', 'DststPipPi0_c0', 'DststPipPim_c0', 'DststPi0Pi0_c0', 'Bp_TauNuDstst_Pip_PUc0', 'B0_TauNuDstst_Pi0_PUc0', 'DstmDsp','DstmD0','DstmDp','BpHc','BmHc','antiB0Hc']

class TimeoutError(Exception):
    pass

class timeout:
    def __init__(self, seconds=1, error_message='Timeout'):
        self.seconds = seconds
        self.error_message = error_message
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

inDic = {}

candDir = 'ntuples_B2DstMu_%s' % NTUPLE_TAG

######## Signals
# inDic['Bd_MuNuDst'] = DSetLoader('Bd_MuNuDst', candDir=candDir)
# inDic['Bd_TauNuDst'] = DSetLoader('Bd_TauNuDst', candDir=candDir)
# ######## D** background
# inDic['Bu_MuDstPi'] = DSetLoader('Bu_MuNuDstPi', candDir=candDir)
# inDic['Bd_MuDstPi'] = DSetLoader('Bd_MuNuDstPi', candDir=candDir)
# inDic['Bd_MuDstPiPi'] = DSetLoader('Bd_MuNuDstPiPi_v3', candDir=candDir)
# inDic['Bu_MuDstPiPi'] = DSetLoader('Bu_MuNuDstPiPi_v3', candDir=candDir)
# inDic['Bu_TauDstPi'] = DSetLoader('Bu_TauNuDstPi', candDir=candDir)
# inDic['Bd_TauDstPi'] = DSetLoader('Bd_TauNuDstPi', candDir=candDir)
# inDic['Bd_TauDstPiPi'] = DSetLoader('Bd_TauNuDstPiPi', candDir=candDir)
# inDic['Bu_TauDstPiPi'] = DSetLoader('Bu_TauNuDstPiPi', candDir=candDir)
# inDic['Bs_MuDstK'] = DSetLoader('Bs_MuNuDstK', candDir=candDir)
# inDic['Bs_TauDstK'] = DSetLoader('Bs_TauNuDstK', candDir=candDir)
#
# ######## D*Hc background
# inDic['Bd_DstDu'] = DSetLoader('Bd_DstDu', candDir=candDir)
# inDic['Bd_DstDd'] = DSetLoader('Bd_DstDd', candDir=candDir)
# inDic['Bd_DstDs'] = DSetLoader('Bd_DstDs', candDir=candDir)
# inDic['Bu_DstDu'] = DSetLoader('Bu_DstDu', candDir=candDir)
# inDic['Bu_DstDd'] = DSetLoader('Bu_DstDd', candDir=candDir)
# inDic['Bs_DstDs'] = DSetLoader('Bs_DstDs', candDir=candDir)
# inDic['Bd_DDs1'] = DSetLoader('Bd_DDs1', candDir=candDir)
# inDic['Bu_DDs1'] = DSetLoader('Bu_DDs1', candDir=candDir)
# inDic['B_DstDXX'] = DSetLoader('B_DstDXX', candDir=candDir)

######## Others background
# inDic['DstKu_KuToMu'] = DSetLoader('DstKu_KuToMu', candDir='ntuples_B2DstMu_211118')

# inDic['JPsiKst'] = DSetLoader('Bd_JpsiKst_General', candDir='ntuples_Bd2JpsiKst_%s' % NTUPLE_TAG)
inDic['JPsiKst'] = DSetLoader('Bd_JpsiKst_General', candDir='ntuples_Bd2JpsiKst_220531')

def getEff(k,N):
    e = k/float(N)
    de = np.sqrt(e*(1-e)/N)
    return [e, de]

# Generator Efficiency

handle = {}
handle['genFilter'] = [Handle('GenFilterInfo'), ('genFilterEfficiencyProducer', '', 'SIM')]
handle['genProduct'] = [Handle('GenLumiInfoProduct'), ('generator', '', 'SIM')]

def analyzeMINIAODs(fileList):
    print 'Analizing', len(fileList), 'MINIAOD'
    N_gen = 0
    N_cuts = 0
    xsec = []
    xsec_err = []
    pb = ProgressBar(maxEntry=len(fileList))
    skippedFiles = []
    for i_j, fileName in enumerate(fileList):
        if not os.path.exists(fileName):
            fileName = 'root://cmsxrootd.fnal.gov/' + fileName
        pb.show(i_j)
        with timeout(seconds=30):
            try:
#                 cmd = 'python generatorEfficiency_MINIAODSIM.py ' + fileName
#                 status, output = commands.getstatusoutput(cmd)
#                 aux = output.split(' ')
#                 N_gen += float(aux[0])
#                 N_cuts += float(aux[1])
#                 xsec.append(float(aux[2]))
#                 xsec_err.append(float(aux[4]))
                for lumi in Lumis(fileName):
                    prods = {}
                    for k,v in handle.iteritems():
                        lumi.getByLabel(v[1], v[0])
                        prods[k] = v[0].product()
                    N_cuts += prods['genFilter'].numEventsPassed()
                    N_gen += prods['genFilter'].numEventsTotal()
                    xs = prods['genProduct'].getProcessInfos()[0].lheXSec()
                    xsec.append(xs.value())
                    xsec_err.append(xs.error())
            except TimeoutError:
                skippedFiles.append(fileName)
    print 'Skipped {} files'.format(len(skippedFiles))
    print 'Total events in analyzed MINIAODs', N_cuts
    xsec = np.array(xsec)
    xsec_err = np.array(xsec_err)
    return N_gen, N_cuts, xsec, xsec_err


N_max = 200
recreate = []
# recreate = inDic.keys()
for n, d in inDic.iteritems():
    print '\n\n--> ' + d.sample

    outdir = os.path.join(d.candLoc, d.full_name)
    outyamlFile = os.path.join(outdir,'effMCgenerator.yaml')
    if os.path.isfile(outyamlFile) and not n in recreate:
        print 'Already present'
        dic = yaml.load(open(outyamlFile, 'r'))
        print dic
        continue

    fileList = d.MINIAOD_filelist
    if N_max > 0 and N_max < len(fileList):
        fileList = np.random.choice(fileList, N_max)
    elif len(fileList) == 0:
        print 'No MiniAODs found, skipping'
        continue

    N_gen, N_cuts, xsec, xsec_err = analyzeMINIAODs(fileList)
    s2 = np.square(xsec_err)
    num = np.sum(xsec/s2)
    den = np.sum(1./s2)
    xsecAvg = 1e3*num/den
    xsecAvg_err = 1e3*np.sqrt(1/den)
    chi2 = np.sum(np.square((xsec - xsecAvg*1e-3)/xsec_err))
    pval = rt.ROOT.Math.chisquared_cdf_c(chi2, len(xsec)-1)
    print 'Chi2: {:.1f}/{} ({:.1f}%)'.format(chi2, len(xsec)-1, pval*100)
    print 'Xsec: {:1.4e} +/- {:1.4e} fb ({:1.1e})'.format(xsecAvg, xsecAvg_err, xsecAvg_err/xsecAvg)
    d.xsec = [xsecAvg, xsecAvg_err]

    e, de = getEff(N_cuts, N_gen)
    print 'eff generator: {:1.3e} +/- {:1.3e} ({:1.1e})'.format(e,de, de/e)
    d.effGEN = [e, de]

    dump_dic = {}
    for k in ['xsec', 'effGEN']:
        aux = getattr(d, k)
        dump_dic[k] = [float(aux[0]), float(aux[1])]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(outyamlFile, 'w') as dumpF:
        dumpF.write(yaml.dump(dump_dic, default_flow_style=False, default_style=''))


# ntuplizer efficiency

for d in inDic.values():
    print '\n\n--> ' + d.sample

    if not os.path.isdir(d.ntuples_dir):
        continue
    cand_out_list = glob(os.path.join(d.ntuples_dir,'out/job*.out'))
    N_analyzed = 0
    N_trg = 0
    N_cand = 0
    print 'Analyzing {} ntuplizer job logs'.format(len(cand_out_list))
    pb = ProgressBar(maxEntry=len(cand_out_list))
    for ic, cand_out in enumerate(cand_out_list):
        pb.show(ic)
        eff_ln = []
        counters = []
        takingCounters = False
        with open(cand_out) as f:
            start_analyzing = False
            for line in f:
                if not start_analyzing and 'efficiency' not in line:
                    continue
                else:
                    start_analyzing = True
                if 'efficiency:' in line:
                    eff_ln.append(line)
                elif 'counters:' in line:
                    takingCounters = True
                elif takingCounters and line[:-1].isdigit():
                    counters.append(int(line[:-1]))
                elif takingCounters:
                    takingCounters = False

        if len(eff_ln) == 0:
            print bcolors.FAIL + "Warning: log file '%s' doesn't have any efficiencies!" % cand_out + bcolors.ENDC
            continue

        aux = re.search('[0-9]+/[0-9]+', eff_ln[0]).group(0)
        aux = aux.split('/')
        N_analyzed += int(aux[1])
        N_trg += int(aux[0])

        aux = re.search(': [0-9]+/', eff_ln[1]).group(0)
        N_cand += int(aux[2:-1])

        counters=np.array(counters)
        if not hasattr(d, 'counters'):
            d.counters = counters
        else:
            d.counters += counters

    d.nTotMINIAOD = N_analyzed
    d.nTotCAND = N_cand
    print 'Total MINIAOD:', N_analyzed
    print 'Total candidates:', N_cand

    e, de = getEff(N_trg, N_analyzed)
    d.effCAND_trg = e, de
    print 'eff candidates (trigger): {:1.3e} +/- {:1.3e} ({:1.1e})'.format(e,de, de/e)
    latexStr = '${:.2f} \pm {:.2f}$'.format(100*e, 100*de)

    e, de = getEff(N_cand, N_trg)
    d.effCAND_cand = e, de
    print 'eff candidates (cand): {:1.3e} +/- {:1.3e} ({:1.1e})'.format(e,de, de/e)
    latexStr += ' & ${:.2f} \pm {:.2f}$'.format(100*e, 100*de)


    e, de = getEff(N_cand, N_analyzed)
    d.effCAND = e, de
    print 'eff candidates: {:1.3e} +/- {:1.3e} ({:1.1e})'.format(e,de, de/e)
    latexStr += ' & ${:.2f} \pm {:.2f}$\\'.format(100*e, 100*de)
    print latexStr

    print 'Getting the total rates (if existing)'
    try:
        fCandLoc = glob(os.path.join(d.ntuples_dir,'out_CAND_*.root'))[0]
        fCand = ur.open(fCandLoc)
        Trate = fCand['p']['Trate']
        d.rate = {}
        for k in Trate.keys():
            r = Trate.array(k)[0]
            r *= 1e12 #GeV -> meV
            d.rate[str(k)] = r
        print 'Done'
    except:
        print bcolors.WARNING + 'Not found' + bcolors.ENDC

    dump_dic = {'nTotMINIAOD': int(d.nTotMINIAOD), 'nTotCAND': int(d.nTotCAND)}
    for k in ['effCAND', 'effCAND_trg', 'effCAND_cand']:
        aux = getattr(d, k)
        dump_dic[k] = [float(aux[0]), float(aux[1])]
    if hasattr(d, 'rate'):
        for k, v in d.rate.iteritems():
            dump_dic['rate_'+k] = float(v)
    with open(os.path.join(d.ntuples_dir,'effCAND.yaml'), 'w') as dumpF:
        dumpF.write(yaml.dump(dump_dic, default_flow_style=False, default_style=''))


t = PrettyTable()
t.field_names = ['Sample'] + [str(i) for i in range(d.counters.shape[0])]
for n, d in inDic.iteritems():
    eff = np.zeros((d.counters.shape[0], 2))
    eff[0] = d.effCAND_trg
    for i in range(d.counters[1:].shape[0]):
        eff[i+1] = getEff(d.counters[i+1], d.nTotMINIAOD)
    t.add_row([n] + ['{:.2f}'.format(100*e[0]) for e in eff])
    x = np.arange(eff.shape[0])
    p = plt.errorbar(x, eff[:, 0], eff[:,1], lw=0, elinewidth=5, label=n)

#     plt.plot(x[[0,-1]], 2*[d.effCAND[0]], '-', color=p[0].get_color())
#     plt.fill_between(x[[0,-1]], 2*[d.effCAND[0]-d.effCAND[1]], 2*[d.effCAND[0]+d.effCAND[1]], color=p[0].get_color(), alpha=0.2)
print t
plt.rcParams.update({'font.size': 20})
plt.xlabel('Counter')
plt.ylabel('Efficiency')
plt.legend(loc='best', numpoints=1)
plt.ylim(0.01,1.05)
plt.xlim(-1, eff.shape[0])
plt.grid(True, which='both')
plt.yscale('log')
plt.gcf().set_size_inches(10, 6)

t = PrettyTable()
t.field_names = ['Sample'] + [str(i) for i in range(d.counters.shape[0])]
for n, d in inDic.iteritems():
    eff = np.zeros((d.counters.shape[0], 2))
    eff[0] = d.effCAND_trg
    for i in range(d.counters[1:].shape[0]):
        eff[i+1] = getEff(d.counters[i+1], d.counters[i])
    t.add_row([n] + ['{:.2f}'.format(100*e[0]) for e in eff])
    x = np.arange(eff.shape[0])
    p = plt.errorbar(x, eff[:, 0], eff[:,1], fmt='o', lw=0, elinewidth=5, label=n)
print t
plt.rcParams.update({'font.size': 20})
plt.xlabel('Counter')
plt.ylabel('Efficiency')
plt.legend(loc='best', numpoints=1)
plt.ylim(0.2,1.05)
plt.xlim(-1, eff.shape[0])
plt.grid(True, which='both')
plt.yscale('log')
plt.gcf().set_size_inches(10, 6)


# Skim Efficiency

for p in order:
    if not p in inDic.keys():
        continue

    s = []
    for c in ['Low', 'Mid', 'High']:
        s.append(inDic[p].printSkimEffLatex(c+'_bare'))
    s = ' & '.join(s)
    s += '\\\\'
#     print p, s
    print s


# Comparison table

# Latex format

latexTable = r'''
\begin{tabular}{|c||c|c||cc|c||cc|cc|cc|}
 \hline
 Process &
 xsec [b] & $\varepsilon_\text{gen}$ &
 $\varepsilon_\text{trg}$ [\%]& $\varepsilon_\text{cand}$ [\%]&
 $\varepsilon_\text{ntp} = \varepsilon_\text{trg}\varepsilon_\text{cand}$ [\%] &
 $\varepsilon_\text{skim}^{low}$ [\%] & $\varepsilon_\text{ntp}\varepsilon_\text{skim}^{low}$ [\%] &
 $\varepsilon_\text{skim}^{mid}$ [\%] & $\varepsilon_\text{ntp}\varepsilon_\text{skim}^{mid}$ [\%] &
 $\varepsilon_\text{skim}^{high}$ [\%] & $\varepsilon_\text{ntp}\varepsilon_\text{skim}^{high}$ [\%] &
 \\
 \hline
 \hline
'''

procTraslation = {
    'mu_c0': r'$B\to D^*\mu\nu$ (hard $b\bar{b}$)',
    'muSoft_c0': r'$B\to D^*\mu\nu$ (soft QCD all)',
}

for n, ds in inDic.iteritems():
    name = procTraslation[n] if n in procTraslation.keys() else n
    fields = [name]
    xsec = '{:.3f}'.format(1e-12*ds.effMCgen['xsec'][0])
    eGen = '{:1.2e}'.format(ds.effMCgen['effGEN'][0])
    fields += [xsec, eGen]
    for k in ['effCAND_trg', 'effCAND_cand', 'effCAND']:
        s = '{:.2f}'.format(100*ds.effCand[k][0])
        fields.append(s)
    for k in ['Low', 'Mid', 'High']:
        s = '{:.2f}'.format(100*ds.getSkimEff(k+'_bare')[0])
        fields.append(s)
    latexTable += ' ' + ' & '.join(fields) + ' \\\\\n \hline\n'
latexTable += r'\end{tabular}' + '\n'
print latexTable

# Latex format transposed

latexTable = r'\begin{tabular}{|'+len(inDic.keys() + ['a'])*r'c|'+r'}\n\hline\n'

fieldNames = [
    r'Process',
    r'xsec [b]',
    r'$\varepsilon_\text{gen}$',
    r'$\varepsilon_\text{trg}$ [\%]',
    r'$\varepsilon_\text{cand}$ [\%]',
    r'$\varepsilon_\text{ntp} = \varepsilon_\text{trg}\varepsilon_\text{cand}$ [\%]',
    r'$\varepsilon_\text{skim}^{low}$ [\%]',
    r'$\varepsilon_\text{ntp}\varepsilon_\text{skim}^{low}$ [\%]',
    r'$\varepsilon_\text{skim}^{mid}$ [\%]',
    r'$\varepsilon_\text{ntp}\varepsilon_\text{skim}^{mid}$ [\%]',
    r'$\varepsilon_\text{skim}^{high}$ [\%]',
    r'$\varepsilon_\text{ntp}\varepsilon_\text{skim}^{high}$ [\%]'
]

procTraslation = {
    'mu_c0': r'$B\to D^*\mu\nu$ (hard $b\bar{b}$)',
    'muSoft_c0': r'$B\to D^*\mu\nu$ (soft QCD all)',
}

latexTable += ' '+fieldNames[0]
for n, ds in inDic.iteritems():
    name = procTraslation[n] if n in procTraslation.keys() else n
    latexTable += ' & ' + name
latexTable += ' \\\\\n \hline\n \hline\n'

latexTable += ' '+fieldNames[1]
for n, ds in inDic.iteritems():
    latexTable += ' & ' + '{:.3f}'.format(1e-12*ds.effMCgen['xsec'][0])
latexTable += ' \\\\\n \hline\n'

latexTable += ' '+fieldNames[2]
for n, ds in inDic.iteritems():
    latexTable += ' & ' + '{:1.2e}'.format(ds.effMCgen['effGEN'][0])
latexTable += ' \\\\\n \hline\n \hline\n'


for i, k in enumerate(['effCAND_trg', 'effCAND_cand', 'effCAND']):
    latexTable += ' '+fieldNames[3+i]
    for n, ds in inDic.iteritems():
        latexTable += ' & ' + '{:.2f}'.format(100*ds.effCand[k][0])
    latexTable += ' \\\\\n \hline\n'
latexTable += ' \hline\n'

for i, k in enumerate(['Low', 'Mid', 'High']):
    latexTable += ' '+fieldNames[6+2*i]
    for n, ds in inDic.iteritems():
        latexTable += ' & ' + '{:.2f}'.format(100*ds.getSkimEff(k+'_bare')[0])
    latexTable += ' \\\\\n \hline\n'
    latexTable += ' '+fieldNames[6+2*i+1]
    for n, ds in inDic.iteritems():
        latexTable += ' & ' + '{:.2f}'.format(100*ds.getSkimEff(k+'_bare')[0]*ds.effCand['effCAND'][0])
    latexTable += ' \\\\\n \hline\n \hline\n'

latexTable += r'\end{tabular}' + '\n'
print latexTable
