# Analysis code for the R(D*) measurement with parked data in CMSSW

This code has to be run after sourceing CMSSW_10_2_3. It does not require compilation but needs the CMSSW environment to be active.

## Suggestions for installation

```
cd ~
mkdir RDstAnalysis
cd RDstAnalysis

cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
cd ../..

git clone https://github.com/ocerri/BPH_RD_Analysis
```

# Setting up Combine

To set up Combine you can run the following commands (see [https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/))

```
cd ~/RDstAnalysis
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
scramv1 b clean; scramv1 b -j12 # always make a clean build
```

# Install the CombineHarvester Tool

```
cd CMSSW_10_2_13/src
cmsenv
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
scram b -j12
```

Before running Combine, make sure to run:

```
source ~/RDstAnalysis/CMSSW_10_2_13; cmsenv
```

# Update .bashrc

You also want to update your ~/.bash_profile and add the following line:

```
ulimit -s unlimited
```

The reason is that at some point combine uses up too much stack space (maybe
printing to stdout?), and so this change is required to avoid a segfault.

# Instructions for reprocessing the data

First, make sure you have installed the ntuplizer using the instructions at
[https://github.com/ocerri/BPH_RDntuplizer](https://github.com/ocerri/BPH_RDntuplizer).

## Running the ntuplizer

```console
$ cd ~/RDstAnalysis/src/ntuplizer/BPH_RDntuplizer/jobSubmission
$ ./run-ntuplizer
```

This will add the ntuplizer jobs to an sqlite database at ~/state.db. Next, you can either submit them all via:

```console
$ submit-condor-jobs --max-jobs 100000
```

or you can set up a cron job like using `crontab -e`:

```
HOME=/storage/af/user/[user]
PATH=/usr/bin:$HOME/local/bin

0 * * * * source $HOME/cmsenv.sh; $HOME/RDstAnalysis/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/jobSubmission/submit-condor-jobs --max-jobs 10000 --loglevel notice --logfile $HOME/submit.log
```

where `cmsenv.sh` is a shell script like:

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export X509_USER_CERT=$HOME/.globus/usercert.pem
export X509_USER_KEY=$HOME/.globus/private/userkey.pem
export X509_USER_PROXY=/tmp/x509up_u${EUID}
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.14.2/Linux-x86_64/bin:${PATH}
export BOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/Boost/1.66.0-f50b5/x86_64-centos7-gcc7-opt/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc48-opt
cd $HOME/RDstAnalysis/CMSSW_10_2_3/
#echo `which cmsenv`
#cmsenv
eval `scramv1 runtime -sh`
cd -
```

If you set up the cron job make sure that your certificate is renewed via:

```
$ voms-proxy-init -voms cms -rfc -valid 192:0 --bits 2048
```

### Checking the status of the jobs

You can check the status of the jobs via:

```
$ condor_q
```

but you can also see via the database:

```
$ sqlite ~/state.db
sqlite> select state, count(state) from ntuplizer_jobs group by state;
RUNNING|9
SUCCESS|31084
```

When all the jobs are done, you now need to run the skimmer.

## Running the skimmer

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/scripts/
$ python B2JpsiKst_skimCAND_v1.py -d '.*' --cat low
$ python B2JpsiKst_skimCAND_v1.py -d '.*' --cat mid
$ python B2JpsiKst_skimCAND_v1.py -d '.*' --cat high
$ python B2DstMu_skimCAND_v1.py -d '.*' --cat low
$ python B2DstMu_skimCAND_v1.py -d '.*' --cat mid
$ python B2DstMu_skimCAND_v1.py -d '.*' --cat high
$ python B2DstMu_skimCAND_v1.py -d 'data' --applyCorr
```

## Running the getEfficiencies.py script

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/scripts/
$ ./generatorEfficiencies.py
```

## Setting up your website

The next few scripts all create plots which are meant to be viewable via a website. If you have a `public_html` folder, then you should set it up:

```console
$ mkdir -p ~/public_html/BPH_RDst
$ cp index.php ~/public_html/BPH_RDst
```

## Getting the trigger scale factors

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/scripts/
$ for t in "Mu7_IP4" "Mu9_IP6" "Mu12_IP6"; do ./triggerEfficiencies.py -v [version] -t $t -d RD --refIP BS; done
$ for t in "Mu7_IP4" "Mu9_IP6" "Mu12_IP6"; do ./triggerEfficiencies.py -v [version] -t $t -d MC --refIP BS; done
```

Next, you need to open up the ipython notebook:

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/analysis/
$ jupyter-notebook --no-browser --port 1234
```

Next, forward port 1234 to your local machine:

```console
$ ssh -NL 1234:localhost:1234 alatorre@login-2.hep.caltech.edu
```

Now you should be able to open the link printed out when you started the
jupyter notebook.

Next, you need to run the notebook `TriggerEfficienciesScaleFactors_v1.ipynb`
for each of the three categories by commenting out the top lines.

# B -> J/psi K\* calibration

To run the B -> J/psi K\* calibration:

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/analysis/
$ for cat in "low" "mid" "high"; do ./kinematicCalibration_Bd_JpsiKst.py -c $cat; done
```

## Running the fit

If you just want to run a single category fit:

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/Combine/
$ ./runCombine.py -c high
```

If you want to run the whole fit, it's better to submit it:

```console
$ cd ~/RDstAnalysis/BPH_RD_Analysis/Combine/
$ for cat in "low" "mid" "high" "comb"; do ./runCombine.py -c $cat -v [version] --submit; done
```

# Utilities

### Notebook extensions

To install notebook extensions:
``pip install jupyter_contrib_nbextensions``
``jupyter contrib nbextension install --user``
In case config parser have erros try to install the specific release
``pip install --user configparser==3.5.1``
