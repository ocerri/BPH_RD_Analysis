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

jobsId = '875644'

cmd = 'condor_q {} -hold'.format(jobsId)
status, output = commands.getstatusoutput(cmd)

failedIds = []
for l in output.split('\n'):
    if re.match('[0-9]+\.[0-9]+[ ]+ocerri', l):
        failedIds.append(l.split(' ')[0].split('.')[1])
print 'Failed jobs:', len(failedIds)
logDir = '/storage/af/user/ocerri/BPH_RD_Analysis/scripts/tmp/B2DstMu_skimCAND_Bd_MuNuDst/out/'
for id in failedIds:
    fn = logDir + 'job_{}_{}.err'.format(id, jobsId)
    with open(fn, 'r') as f:
        lines = f.readlines()
        errorSolved = False
        for l in lines:
            if 'Error in <TBranchElement::GetBasket>: File:' in l:
                candName = l[l.find('/storage/af'): l.find('.root')+5]
                newLoc = os.path.join(os.path.dirname(candName), 'corruptedFiles', os.path.basename(candName))
                print id, ': moving', os.path.basename(candName)[9:-5]
                cmd = 'mv '+candName+' '+newLoc
                os.system(cmd)
                errorSolved = True
                break
        if not errorSolved:
            print lines
            raise
