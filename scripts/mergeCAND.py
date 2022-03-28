import sys, os, re
from glob import glob
import commands
from multiprocessing import Pool
import numpy as np

redF = 50

fileTemplate = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu_220326'
fileTemplate += '/out_CAND_*.root'

# Create merge directory
outdir = os.path.join(os.path.dirname(fileTemplate), 'merged')
if os.path.exists(outdir):
    print 'Output direcotry already exists'
    print outdir
    exit()

cmd = 'mkdir ' + outdir
os.system(cmd)

fileList = glob(fileTemplate)
nOutput = len(fileList)/redF
if len(fileList) % redF > 0:
    nOutput += 1
print 'Merging {} files into {} files'.format(len(fileList), nOutput)
# i = 0
# while i*redF < len(fileList):
#     iStart = i*redF
#     iStop = (i+1)*redF
#     # print 'Merging {}: {} - {}'.format(i, iStart, iStop)
#     # cmd = 'hadd {}/out_CAND_m{}.root '.format(outdir, i) + ' '.join(fileList[iStart: iStop])
#     # status, output = commands.getstatusoutput(cmd)
#
#     i += 1
#
# print i-1
# exit()


def mergeFun(i):
    iStart = i*redF
    iStop = (i+1)*redF
    print 'Merger {} ({} - {}): {}'.format(i, iStart, iStop, len(fileList[iStart: iStop]))
    cmd = 'hadd {}/out_CAND_m{}.root '.format(outdir, i) + ' '.join(fileList[iStart: iStop])
    status, output = commands.getstatusoutput(cmd)
    if status:
        print output
        raise
    return len(fileList[iStart: iStop])

p = Pool(20)
outputs = p.map(mergeFun, list(range(nOutput)))

print 'Merged {} files'.format(np.sum(outputs))
