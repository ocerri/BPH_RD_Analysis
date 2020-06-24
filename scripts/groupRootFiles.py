from glob import glob
import argparse, os, re
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument ('inputFiles', help='Local path of the files to merge', nargs='+', type=str)
parser.add_argument ('--CMSSW_path', help='Path of the cmssw to be used', default='/storage/user/ocerri/CMSSW_10_2_3')
parser.add_argument ('-n', '--nInPerOut', help='Number of files to be grouped', default=10, type=int)
parser.add_argument ('--edm', help='Flag for edm files', default=False, action='store_true')
parser.add_argument ('--parallelType', choices=[None, 'pool'], default=None, help='Kind of parallelization')
args = parser.parse_args()

if len(args.inputFiles) == 1:
    args.inputFiles = glob(args.inputFiles[0])
nFileIn = len(args.inputFiles)
print 'Total number of files:', nFileIn

nFileOut = nFileIn / args.nInPerOut
if nFileIn % args.nInPerOut > 1:
    nFileOut += 1
print 'Merging in gropus of {} for a total of {}'.format(args.nInPerOut, nFileOut)

#Create file groups
inFileGroups = []
for i in range(nFileOut):
    i_start = i*args.nInPerOut
    if i == nFileOut-1:
        i_stop = len(args.inputFiles)
    else:
        i_stop = (i+1)*args.nInPerOut
    inFileGroups.append(args.inputFiles[i_start : i_stop])
print inFileGroups

os.path.
