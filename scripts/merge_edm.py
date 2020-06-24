from glob import glob
import argparse, os, re

parser = argparse.ArgumentParser()
parser.add_argument ('local_path', help='Local path of the files to merge', nargs='+')
parser.add_argument ('--CMSSW_path', help='absolute path of the cmssw to be used', default='/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/CMSSW_10_2_3')
parser.add_argument ('--N_file', help='Max number of files', default=-1, type=int)
args = parser.parse_args()

pwd = os.environ['PWD']

input_file_list = args.local_path
merge_list_path = pwd + '/' + os.path.dirname(input_file_list[0]) + '/merge_list.txt'
print merge_list_path
n_min = 9999999999
n_max = -1

f = open(merge_list_path, 'w')

N = 0
for fname in input_file_list:
    out = re.search('_[0-9]+.root', os.path.basename(fname))
    n = int(out.group(0)[1:-5])
    if n < n_min:
        n_min = n

    if n > n_max:
        n_max = n

    f.write('file:'+pwd+'/'+fname+'\n')
    N += 1
    if args.N_file>0 and N == N_file:
        break

f.close()

outname_root = pwd + '/' + fname.replace('_{}.root'.format(n), '_merged_{}-{}.root'.format(n_min, n_max))
print outname_root

cmd = 'cd {}/src; eval `scramv1 runtime -sh`;'.format(args.CMSSW_path)
cmd += ' edmCopyPickMerge inputFiles_load={} outputFile={}'.format(merge_list_path, outname_root)
print cmd
os.system(cmd)
