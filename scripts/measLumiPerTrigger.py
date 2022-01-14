#!/usr/bin/env python
# coding: utf-8

# Works only with work/CMSSW_10_2_3/src

# In[1]:


import os, re, json
import pickle
import commands
from glob import glob
import numpy as np
from prettytable import PrettyTable
import yaml


# In[2]:


os.system('export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH')


# # Input data

# In[3]:


# loc = os.environ['HOME'] + '/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/jobSubmission/tmp'
# loc += '/crab_ParkingBPH*_RDntuplizer_B2JpsiKst_200622'
# print loc
# loc = glob(loc)
# print len(loc)


# In[4]:


loc = os.environ['HOME'] + '/RDstAnalysis/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/jobSubmission/tmp'
loc += '/crab_ParkingBPH[1235]_*_RDntuplizer_B2DstMu_220110'
print loc
loc = glob(loc)
print len(loc)


# In[5]:


# loc = os.environ['HOME'] + '/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/jobSubmission/tmp'
# loc += '/crab_ParkingBPH*_RDntuplizer_combDmstMum_200611'
# print loc
# loc = glob(loc)
# print len(loc)


# In[6]:


nameTemplate = '/storage/af/group/rdst_analysis/BPhysics/data/cmsRD_mediumId_lostInnerHits/lumiReport/{}_brilcalcPerTrigger.yaml'
recreate = True

lumiInfo = {}
for main_dir in loc:
    print main_dir.split('/')[-1]
    s = main_dir.split('/')[-1]
    idx = s.find('Run2018')
    era = s[idx+len('Run2018')]
    idx = s.find('ParkingBPH')
    part = int(s[idx+len('ParkingBPH')])
    print era, part

    if not era in lumiInfo.keys():
        lumiInfo[era] = {}
        
    name = os.path.basename(main_dir)[5:]
    fname = nameTemplate.format(name)
    if os.path.isfile(fname) and not recreate:
        lumiDic = yaml.load(open(fname, 'rb'))
    else:
        cmd = 'brilcalc lumi -u /fb --precision 2f'
        cmd += ' --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json'
        cmd += ' --hltpath HLT_Mu\*_IP?_part' + str(part-1) + '_v?'
        cmd += ' -i {}/processedLumis.json'.format(main_dir + '/results')
        cmd += ' -c web'
        
        print(cmd)
        status, output = commands.getstatusoutput(cmd)
        print(output)
        aux = output.split('#Summary:')[1]
        aux = aux.split('\n')
        lumiDic = {}
        for l in aux:
            if not l.startswith('| HLT'):
                continue
            l = l[1:-1]
            content = l.split(' | ')
            trgPath = content[0].replace(' ', '')
            lumi = float(content[-1].replace(' ', ''))
            lumiDic[trgPath] = lumi
            
        lumi_str = yaml.dump(lumiDic, default_flow_style=False)
        f = open(fname, 'w')
        f.write(lumi_str)
        f.close()
            
    for trgPath, lumi in lumiDic.iteritems():
        if not '_part'+str(part-1) in trgPath:
            continue
        trgPath = trgPath.replace('_part' + str(part-1), '')          
        if lumi < 0.02:
            continue
        if not trgPath in lumiInfo[era].keys():
            lumiInfo[era][trgPath] = np.zeros(6)
        lumiInfo[era][trgPath][part-1] = lumi
        print trgPath, part, lumi
    print '\n'


# In[7]:


lumiInfo


# In[8]:


eras = list(np.sort(lumiInfo.keys()))
table = PrettyTable()

table.field_names = ['HLT path'] + eras

printed_paths = []
for e in eras:
    for t in lumiInfo[e]:
        t = t[:-3]
        if t in printed_paths:
            continue
        printed_paths.append(t)
        row = [t]
        for e in eras:
            if t in [path[:-3] for path in lumiInfo[e].keys()]:
                for path in lumiInfo[e].keys():
                    if t == path[:-3]:
                        row.append(str(lumiInfo[e][path]))
                        break
            else:
                row.append(' - ')
        table.add_row(row)
print table


# In[9]:


eras = list(np.sort(lumiInfo.keys()))
table = PrettyTable()

l = ['pT', 'IP']
for e in eras:
    l += [e, 'Tot ' + e]
l += ['Tot']
table.field_names = l

printed_paths = []
for e in eras:
    for t in lumiInfo[e]:
        tot_lumi = 0
        t = t[:-3]
        if t in printed_paths:
            continue
        printed_paths.append(t)
        row = [t[6:-4], t[-1]]
        for e in eras:
            if t in [path[:-3] for path in lumiInfo[e].keys()]:
                for path in lumiInfo[e].keys():
                    if t == path[:-3]:
                        row.append(str(lumiInfo[e][path]))
                        
                        s = np.sum(lumiInfo[e][path])
                        tot_lumi += s
                        row.append('{:.2f}'.format(s))
                        break
            else:
                row.append(' - ')
                row.append(' - ')
        row.append('{:.2f}'.format(tot_lumi))
        table.add_row(row)
print table

