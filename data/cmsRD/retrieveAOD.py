import commands
import random
import time

import signal

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

# fileList = open('/storage/af/user/ocerri/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/production/inputFiles_ParkingBPH1_Run2018D-05May2019promptD-v1_AOD.txt', 'r').readlines()
fileList = open('/storage/af/user/ocerri/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/production/inputFiles_ParkingBPH1_Run2018D-v1_RAW.txt', 'r').readlines()
random.shuffle(fileList)

for i, fn in enumerate(fileList):
    fn = fn[:-1]
    print '{}/{}: '.format(i+1, len(fileList))+fn
    cmd = 'xrdcp root://cmsxrootd.fnal.gov/{} ./'.format(fn)
    status, output = commands.getstatusoutput(cmd)
    print 'Returned status:', status
    if '[ERROR]' in output:
        for ln in output.split('\n'):
            if '[ERROR]' in ln:
                print ln
    else:
        print output
        print 'Transfer succedded!'
        break


# with timeout(seconds=30):
#     try:
#         pass
#     except TimeoutError:
#         pass
