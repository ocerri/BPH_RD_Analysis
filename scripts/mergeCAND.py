import sys, os, re
from glob import glob
import commands

redF = 50

fileTemplate = '/storage/af/group/rdst_analysis/BPhysics/data/cmsMC/CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples_B2DstMu_wOC'
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
print 'Merging {} files'.format(len(fileList))
i = 0
while i*redF < len(fileList):
    iStart = i*redF
    iStop = (i+1)*redF
    print 'Merging {}: {} - {}'.format(i, iStart, iStop)
    cmd = 'hadd {}/out_CAND_{}.root '.format(outdir, i) + ' '.join(fileList[iStart: iStop])
    status, output = commands.getstatusoutput(cmd)

    i += 1
