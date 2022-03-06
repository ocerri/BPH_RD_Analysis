#!/usr/bin/env python

import os, commands, json
import random, datetime
import argparse
parser = argparse.ArgumentParser(description='Dump condor jobs on hold', add_help=True)
parser.add_argument ('--output', '-o', default=os.path.join(os.environ['HOME'],'public_html/condor_held_'+os.environ['HOSTNAME'].split('.')[0]), type=str, help='Output directory.')
parser.add_argument ('--jobsStatus', default='held', choices=['run', 'held'], type=str, help='Job status.')
parser.add_argument ('--maxDump', default=10, type=int, help='Max numberf jobs to dump')
parser.add_argument ('--verbose', default=False, action='store_true', help='Run verbosity.')
args = parser.parse_args()

while args.output.endswith('/'):
        args.output = args.output[:-1]
verbose = args.verbose

if verbose:
	print 'Output:', args.output

if os.path.isdir(args.output):
	os.system('rm -rf '+args.output)
if os.path.isfile(args.output+'.log'):
	os.system('rm -rf '+args.output+'.log')

cmd = 'condor_q --'+args.jobsStatus+' -af ClusterId ProcId JobBatchName UserLog'
status, output = commands.getstatusoutput(cmd)

if status:
	if verbose:
		print '[ERROR] condor_q status', status
		print output
	else:
		with open(args.output+'.log', 'w') as f:
			f.write(output)
	exit()

os.makedirs(args.output)
os.system('cp {o}/../index.php {o}/'.format(o=args.output))

flog = open(os.path.join(args.output, 'log.log'), 'w')
flog.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
jobs = output.split('\n')
flog.write('Total number of jobs in status ' + args.jobsStatus + ': '+str(len(jobs)-1) + '\n')
flog.write('Copying files for max '+str(args.maxDump) + ' jobs\n')
flog.write('\n\nFull output:\n')
flog.write(50*'-' + '\n')
flog.write(output)
flog.write(50*'-' + '\n')

if not output:
	exit()
print 'Total jobs {}: {}'.format(args.jobsStatus, len(jobs))
if len(jobs) > args.maxDump:
    print 'Randomly selecting', args.maxDump, 'jobs to dump'
    jobs = random.sample(jobs, args.maxDump)

for job in sorted(jobs):
	cid, pid, name, logLoc = job.split(' ')
	outname = os.path.join(args.output, '__'.join([name, cid, pid]))
	target = logLoc.replace('.log', '')
	if verbose:
		print cid, pid, name
	for ext in ['log', 'err', 'out']:
		status, output = commands.getstatusoutput('cp -f {t}.{e} {o}.{e}'.format(t=target, o=outname, e=ext))
		if verbose:
			print '\t', ext, status

flog.close()
