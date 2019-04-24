import os, argparse, re
from glob import glob
import time
# import matplotlib.pyplot as plt
import numpy as np
from ShowProgressBar import ShowProgressBar
from prettytable import PrettyTable

def getTime_job(line):
    time_format = '%Y/%m/%d %H:%M:%S'
    out = re.search('\) .+ Job', line).group(0)[2:-4]
    return time.mktime(time.strptime('2019/'+out, time_format))

def getTime_record(line):
    time_format = '%d-%b-%Y %H:%M:%S'
    out = re.search(' at .+\.', line).group(0)[4:-1]
    return time.mktime(time.strptime(out, time_format))

parser = argparse.ArgumentParser()
parser.add_argument ('job_out_dir', help='location of the job_output directory')
parser.add_argument ('--clusterID', help='HTCondor cluster ID', default=-1)
parser.add_argument ('--N_triggered', help='Number of trigger events', default=-1)
args = parser.parse_args()

d = args.job_out_dir
if not d.endswith('/'):
    d += '/'

N_jobs = {}
N_jobs['sub'] = -1
N_req_evts_per_job = -1
with open(d + 'cfg/jobs.sub', 'r') as f:
    lns = f.readlines()
    last_line = lns[-1]
    if last_line.startswith('queue '):
        N_jobs['sub'] = int(last_line[6:])

    for l in lns:
        if l.startswith('arguments'):
            out = re.search('= [0-9]+ ', l).group(0)
            N_req_evts_per_job = int(out[2:-1])

glb_path = d + 'out/'
if args.clusterID >= 0:
    glb_path += '*.{}.*.log'.format(args.clusterID)
else:
    glb_path += '*.log'

fpath_list = glob(glb_path)
if len(fpath_list) < 1:
    print 'No files found in ', d
    if args.clusterID >= 0:
        print 'For clusterID', args.clusterID
    exit()

N_jobs['ukill'] = 0
N_jobs['done'] = 0
N_jobs['overtime'] = 0
N_jobs['unkown'] = 0
N_jobs['trouncated'] = 0
N_jobs['missing log'] = 0
N_jobs['missing out'] = 0
N_jobs['missing err'] = 0

total_job_time_completed = []
total_job_time_uncompleted = []
job_N_Pythia_selected = []
job_N_Pythia_generated = []
job_N_cmssw_sel = []

job_time_step = [[], [], [], [], []]

print 'Checking log, out and err files'
for ifn, fn in enumerate(fpath_list):
    ShowProgressBar(ifn, len(fpath_list))

    #Get the total time info from the log file
    with open(fn, 'r') as f:
        try:
            lns = f.readlines()
        except:
            N_jobs['missing log'] +=1
            continue

        sub_time = None
        start_time = None
        stop_time = None
        for il, l in enumerate(lns):
            if sub_time is None and 'Job submitted from host' in l:
                sub_time = getTime_job(l)
            elif start_time is None and 'Job executing on host' in l:
                start_time = getTime_job(l)
            elif stop_time is None:
                if 'Job terminated.' in l and 'Normal termination' in lns[il+1]:
                    N_jobs['done'] += 1
                    stop_time = getTime_job(l)
                    total_job_time_completed.append(stop_time-start_time)
                    break
                elif 'Job was aborted' in l:
                    stop_time = getTime_job(l)
                    if lns[il+1].endswith('wall time exceeded allowed max.\n'):
                        N_jobs['overtime'] += 1
                    elif 'via condor_rm' in lns[il+1]:
                        N_jobs['ukill'] += 1
                    total_job_time_uncompleted.append(stop_time-start_time)
                    break


        if stop_time is None:
            if l == '...\n'or len(l) < 2:
                N_jobs['trouncated'] += 1
            else:
                print '\nUnkwonw source of error in file:', fn
                print 'Last line:', l
                N_jobs['unkown'] += 1
            continue

        if 'Job was aborted' in l:
            continue

    #Get the pythia numbers from the out file
    with open(fn.replace('.log', '.out'), 'r') as f:
        try:
            lns = f.readlines()
        except:
            N_jobs['missing out'] +=1
            continue
        N_Pythia_selected = -1
        for il, l in enumerate(reversed(lns)):
            if '*-------  End PYTHIA Event and Cross Section Statistics' in l:
                N_Pythia_selected = int([x for x in lns[-1-il-2].split(' ') if x][5])
            elif '*-------  PYTHIA Event and Cross Section Statistics' in l:
                if 'Pythia::next():' in lns[-1-il-2]:
                    aux = int(re.search(': [0-9]+', lns[-1-il-2]).group(0)[2:])
                    job_N_Pythia_generated.append(aux + 500)
                    break
        job_N_Pythia_selected.append(int(N_Pythia_selected))

    #Get the signle step time & evts info from the err file
    with open(fn.replace('.log', '.err'), 'r') as f:
        try:
            lns = f.readlines()
        except:
            N_jobs['missing err'] +=1
            continue
        t_step = []
        for i in range(5):
            t_step.append([None, None])



        for il, l in enumerate(lns):
            if t_step[1][0] == None:
                if l.startswith('Begin processing the 1st record'):
                    t_step[1][0] = getTime_record(l)
            elif t_step[1][1] == None:
                if l.startswith('GenXsecAnalyzer'):
                    j=0
                    while not lns[il-3-j].startswith('Begin processing'):
                        j += 1
                    t_step[1][1] = getTime_record(lns[il-3-j][:-1])

            elif t_step[2][0] == None:
                if l.startswith('Begin processing the 1st record'):
                    t_step[2][0] = getTime_record(l)
            elif t_step[2][1] == None:
                if '  Closed file file:' in l:
                    j=0
                    while not lns[il-1-j].startswith('Begin processing the'):
                        j+=1
                    t_step[2][1] = getTime_record(lns[il-1-j])
            elif t_step[3][0] == None:
                if l.startswith('Begin processing the 1st record'):
                    t_step[3][0] = getTime_record(l)
            elif t_step[3][1] == None:
                if '  Closed file file:' in l:
                    j=0
                    while not lns[il-1-j].startswith('Begin processing the'):
                        j+=1
                    t_step[3][1] = getTime_record(lns[il-1-j])
            elif t_step[4][0] == None:
                if l.startswith('Begin processing the 1st record'):
                    t_step[4][0] = getTime_record(l)
            elif t_step[4][1] == None:
                if '  Closed file file:' in l:
                    j=0
                    while not lns[il-1-j].startswith('Begin processing the'):
                        j+=1
                    t_step[4][1] = getTime_record(lns[il-1-j])
                    N = re.search('Begin processing the [0-9]+', lns[il-1-j]).group(0)[21:]
                    job_N_cmssw_sel.append(int(N))
                    break

        for i in range(1,5):
            if (not t_step[i][1] is None) and (not t_step[i][0] is None):
                job_time_step[i].append(t_step[i][1] - t_step[i][0])
            else:
                job_time_step[i].append(0)
                job_N_cmssw_sel.append(0)

N_triggered = args.N_triggered
if N_triggered == -1:
    if os.path.isfile(d+'cfg/trigger_report.txt'):
        fin = open(d+'cfg/trigger_report.txt', 'r')
        N_triggered = int(fin.readlines()[0][22:])
    else:
        #Count number of triggered events
        cmd = 'cd /afs/cern.ch/user/o/ocerri/cernbox/BPhysics/CMSSW_10_2_3/src/ntuplizer/BPH_RD;'
        cmd += ' eval `scramv1 runtime -sh`;'
        cmd += ' python FWLite_master.py config/cfg_countBPH_trg.py -i "{}*SIM_[0-9]*.root"'.format(d)
        print 'Now running:\n'+cmd
        output = os.popen(cmd).read()

        for ln in output.split('\n'):
            if ln[:15] == 'Accepted events':
                N_triggered = int(ln[16:])
                print 'Trg evts:', N_triggered
                fout = open(d+'cfg/trigger_report.txt', 'w')
                fout.write('BPH triggered events: ' + ln[16:])
                fout.close()
                break

frep = open(d+'jobs_report.txt','w')
job_rep_table = PrettyTable()
for k, v in N_jobs.iteritems():
    job_rep_table.add_column(k, [v, '{:.1f}%'.format(100*float(v)/N_jobs['sub'])])
print '\nJobs report'
frep.write('\nJobs report\n')
print job_rep_table
frep.write(job_rep_table.get_string()+'\n')


print '\nEvent number report'
frep.write('\nEvent number report\n')
print 'Submitted: {}*{} = {:1.1e}'.format(N_jobs['sub'], N_req_evts_per_job, N_jobs['sub']*N_req_evts_per_job)
frep.write('Submitted: {}*{} = {:1.1e}\n'.format(N_jobs['sub'], N_req_evts_per_job, N_jobs['sub']*N_req_evts_per_job))
evt_N_table = PrettyTable()
evt_N_table.field_names = ['', 'Pythia Gen', 'Pythia sel', 'CMSSW sel', 'BPH trg']
evt_N_table.add_row(['Job avg',
                     '{:.0f}'.format(np.mean(job_N_Pythia_generated)),
                     '{:.0f}'.format(np.mean(job_N_Pythia_selected)),
                     '{:.0f}'.format(np.mean(job_N_cmssw_sel)), '{:.0f}'.format(float(N_triggered)/len(glob(d+'*SIM_[0-9]*.root')))
                     ])
evt_N_table.add_row(['Tot', np.sum(job_N_Pythia_generated), np.sum(job_N_Pythia_selected), np.sum(job_N_cmssw_sel), N_triggered])
evt_N_table.add_row(['Eff prev step',
                     '-',
                     '{:.2f}%'.format(100*np.sum(job_N_Pythia_selected)/float(np.sum(job_N_Pythia_generated))),
                     '{:.2f}%'.format(100*np.sum(job_N_cmssw_sel)/float(np.sum(job_N_Pythia_selected))),
                     '{:.2f}%'.format(100*N_triggered/float(np.sum(job_N_cmssw_sel)))
                    ])
N_req = float(N_jobs['sub']*N_req_evts_per_job)
evt_N_table.add_row(['Eff requested', '{:.2f}%'.format(100*np.sum(job_N_Pythia_generated)/N_req), '{:.2f}%'.format(100*np.sum(job_N_Pythia_selected)/N_req), '{:.2f}%'.format(100*np.sum(job_N_cmssw_sel)/N_req), '{:.2f}%'.format(100*N_triggered/N_req)])
print evt_N_table
frep.write(evt_N_table.get_string()+'\n')


job_time_step_a = np.array(job_time_step[1:]).T
total_job_time_completed = np.array(total_job_time_completed)
job_IO_time = total_job_time_completed - np.sum(job_time_step_a, axis=1)

time_table = PrettyTable()
time_table.field_names = ['', 'Tot', 'GEN-SIM', 'DIGI (+PU)', 'RECO', 'MINIAOD', 'Other']
time_table.add_row(['Tot [h]', '{:.1f}'.format(np.sum(total_job_time_completed)/3600.)]+
                   ['{:.1f}'.format(np.sum(job_time_step[i])/3600.) for i in range(1,5)]+
                   ['{:.1f}'.format(np.sum(job_IO_time)/3600.)]
                   )
time_table.add_row(['Job avg [h]', '{:.1f}'.format(np.mean(total_job_time_completed)/3600.)]+
                   ['{:.1f}'.format(np.mean(job_time_step[i])/3600.) for i in range(1,5)]+
                   ['{:.1f}'.format(np.mean(job_IO_time)/3600.)]
                   )
time_table.add_row(['Job 70% [h]', '{:.1f}'.format(np.percentile(total_job_time_completed, 70)/3600.)]+
                   ['{:.1f}'.format(np.percentile(job_time_step[i], 70)/3600.) for i in range(1,5)]+
                   ['{:.1f}'.format(np.percentile(job_IO_time, 70)/3600.)]
                   )
time_table.add_row(['Avg / trg evt [s]', '{:.1f}'.format(np.sum(total_job_time_completed)/float(N_triggered))]+
                   ['{:.1f}'.format(np.sum(job_time_step[i])/float(N_triggered)) for i in range(1,5)]+
                   ['{:.1f}'.format(np.sum(job_IO_time)/float(N_triggered))]
                   )

time_table.add_row(['Avg / req evt [s]', '{:.2f}'.format(np.sum(total_job_time_completed)/N_req)]+
                   ['{:.3f}'.format(np.sum(job_time_step[i])/N_req) for i in range(1,5)]+
                   ['{:.3f}'.format(np.sum(job_IO_time)/N_req)]
                   )


print '\nTime table'
frep.write('\nTime table\n')
print time_table
frep.write(time_table.get_string()+'\n')


frep.close()
