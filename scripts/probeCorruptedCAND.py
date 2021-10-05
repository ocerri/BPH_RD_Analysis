import sys, os, re
from glob import glob
import commands

sys.path.append('../lib')
sys.path.append('../analysis')
from progressBar import ProgressBar

import ROOT as rt
# rt.gErrorIgnoreLevel = rt.kError
rt.gErrorIgnoreLevel = rt.kFatal
# rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.ERROR)

fileDir = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu_wOC/'

fileList = glob(fileDir + '/out_CAND_*.root')
print 'Probing {} files'.format(len(fileList))

corrupted = []

pb = ProgressBar(maxEntry=len(fileList))
for i, fn in enumerate(fileList):
    pb.show(i)
    tree = rt.TChain('outA/Tevts')
    try:
        tree.Add(fn)
        for j in range(tree.GetEntry()):
            tree.GetEntry(j, 1)
            ev = tree
            aux = ev.pval_piK.size()
    except:
        print fn, 'corrupted'
        corrupted.append(fn)

print 'Corrupted files found:', len(corrupted)

for fn in corrupted:
    id = os.path.basename(fn)[9:-5]
    print id

    # newLoc = os.path.join(os.path.dirname(fn), 'corruptedFiles', os.path.basename(fn))
    # cmd = 'mv '+fn+' '+newLoc
    # os.system(cmd)
