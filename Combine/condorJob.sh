#!/bin/bash
directory=$1
command="${@:2}"

cd /storage/af/user/ocerri/work/CMSSW_10_2_13/src/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo "cd $directory"
cd $directory
echo "$command"
echo " "; echo " "; echo " "
echo "======================== JOB START ========================"
$command
exitcode=$?
echo " "; echo " "; echo " "
echo "======================== JOB DONE ========================"
echo "Exit code: $exitcode"
exit $exitcode
