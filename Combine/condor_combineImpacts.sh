#!/bin/sh
ulimit -s unlimited
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /storage/user/ocerri/work/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
eval `scramv1 runtime -sh`
cd /storage/user/ocerri/BPhysics/Combine

if [ $1 -eq 0 ]; then
  combine -M MultiDimFit -n _paramFit_v7_B0pT --algo impact --redefineSignalPOIs r -P B0pT --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 1 ]; then
  combine -M MultiDimFit -n _paramFit_v7_B2DstCLNR0 --algo impact --redefineSignalPOIs r -P B2DstCLNR0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 2 ]; then
  combine -M MultiDimFit -n _paramFit_v7_B2DstCLNR1 --algo impact --redefineSignalPOIs r -P B2DstCLNR1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 3 ]; then
  combine -M MultiDimFit -n _paramFit_v7_B2DstCLNR2 --algo impact --redefineSignalPOIs r -P B2DstCLNR2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 4 ]; then
  combine -M MultiDimFit -n _paramFit_v7_B2DstCLNRhoSq --algo impact --redefineSignalPOIs r -P B2DstCLNRhoSq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 5 ]; then
  combine -M MultiDimFit -n _paramFit_v7_DstPi0Br --algo impact --redefineSignalPOIs r -P DstPi0Br --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 6 ]; then
  combine -M MultiDimFit -n _paramFit_v7_DstPi0Pi0Br --algo impact --redefineSignalPOIs r -P DstPi0Pi0Br --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 7 ]; then
  combine -M MultiDimFit -n _paramFit_v7_DstPipBr --algo impact --redefineSignalPOIs r -P DstPipBr --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 8 ]; then
  combine -M MultiDimFit -n _paramFit_v7_DstPipPi0Br --algo impact --redefineSignalPOIs r -P DstPipPi0Br --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 9 ]; then
  combine -M MultiDimFit -n _paramFit_v7_DstPipPimBr --algo impact --redefineSignalPOIs r -P DstPipPimBr --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 10 ]; then
  combine -M MultiDimFit -n _paramFit_v7_HcBr --algo impact --redefineSignalPOIs r -P HcBr --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 11 ]; then
  combine -M MultiDimFit -n _paramFit_v7_HcmixD02mu --algo impact --redefineSignalPOIs r -P HcmixD02mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 12 ]; then
  combine -M MultiDimFit -n _paramFit_v7_HcmixDs2mu --algo impact --redefineSignalPOIs r -P HcmixDs2mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 13 ]; then
  combine -M MultiDimFit -n _paramFit_v7_b2B0Had --algo impact --redefineSignalPOIs r -P b2B0Had --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 14 ]; then
  combine -M MultiDimFit -n _paramFit_v7_b2BpHad --algo impact --redefineSignalPOIs r -P b2BpHad --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 15 ]; then
  combine -M MultiDimFit -n _paramFit_v7_muBr --algo impact --redefineSignalPOIs r -P muBr --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 16 ]; then
  combine -M MultiDimFit -n _paramFit_v7_muonIdSF --algo impact --redefineSignalPOIs r -P muonIdSF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 17 ]; then
  combine -M MultiDimFit -n _paramFit_v7_pAddTk --algo impact --redefineSignalPOIs r -P pAddTk --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 18 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 19 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 20 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin10 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 21 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin11 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 22 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin12 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 23 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin13 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 24 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin14 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 25 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin15 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 26 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin16 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 27 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin17 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 28 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin18 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 29 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin19 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin19 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 30 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 31 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin20 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin20 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 32 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin21 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin21 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 33 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin22 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin22 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 34 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin23 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin23 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 35 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin24 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin24 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 36 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin25 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin25 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 37 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin26 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin26 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 38 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin27 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin27 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 39 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin28 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin28 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 40 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin29 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin29 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 41 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 42 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 43 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 44 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 45 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 46 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 47 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_m_mHad_bin9 --algo impact --redefineSignalPOIs r -P prop_binAddTk_m_mHad_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 48 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 49 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_DstPipPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_DstPipPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 50 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 51 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_Hc --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_Hc --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 52 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 53 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin0_tau --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin0_tau --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 54 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 55 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_DstPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_DstPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 56 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 57 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_DstPip --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_DstPip --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 58 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 59 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_Hc --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_Hc --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 60 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 61 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin10_tau --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin10_tau --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 62 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_DstPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_DstPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 63 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 64 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_DstPip --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_DstPip --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 65 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_DstPipPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_DstPipPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 66 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 67 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin11_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin11_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 68 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin12_DstPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin12_DstPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 69 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin12_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin12_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 70 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin12_DstPipPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin12_DstPipPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 71 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin12_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin12_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 72 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin12_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin12_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 73 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin13_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin13_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 74 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin13_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin13_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 75 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin13_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin13_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 76 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin13_tau --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin13_tau --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 77 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin14_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin14_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 78 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin14_Hc --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin14_Hc --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 79 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin14_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin14_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 80 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 81 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 82 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 83 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 84 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 85 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 86 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 87 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_DstPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_DstPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 88 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_DstPi0Pi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_DstPi0Pi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 89 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 90 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_Hc --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_Hc --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 91 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 92 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_mm_mHad_bin9_tau --algo impact --redefineSignalPOIs r -P prop_binAddTk_mm_mHad_bin9_tau --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 93 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 94 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 95 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin10 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 96 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin11 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 97 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin12 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 98 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin13 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 99 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin14 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 100 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin15 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 101 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin16 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 102 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin17 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 103 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin18 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 104 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin19 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin19 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 105 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 106 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin20 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin20 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 107 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin21 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin21 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 108 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin22 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin22 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 109 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin23 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin23 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 110 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin24 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin24 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 111 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin25 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin25 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 112 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin26 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin26 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 113 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin27 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin27 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 114 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin28 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin28 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 115 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin29 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin29 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 116 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 117 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin30 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin30 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 118 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin31 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin31 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 119 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin32 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin32 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 120 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin33 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin33 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 121 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin34 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin34 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 122 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 123 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 124 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 125 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 126 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 127 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_p_mHad_bin9 --algo impact --redefineSignalPOIs r -P prop_binAddTk_p_mHad_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 128 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 129 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 130 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin10 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 131 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin11 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 132 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin12 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 133 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin13 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 134 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin14 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 135 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin15 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 136 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin16 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 137 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin17 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 138 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin18 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 139 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin19 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin19 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 140 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 141 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin20 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin20 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 142 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin21 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin21 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 143 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin22 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin22 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 144 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin23 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin23 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 145 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin24 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin24 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 146 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin25 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin25 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 147 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin26 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin26 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 148 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_DstPip --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_DstPip --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 149 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_DstPipPi0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_DstPipPi0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 150 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_DstPipPim --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_DstPipPim --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 151 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_Hc --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_Hc --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 152 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_mu --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_mu --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 153 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin27_tau --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin27_tau --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 154 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin28 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin28 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 155 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin29 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin29 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 156 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 157 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 158 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 159 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 160 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 161 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 162 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mHad_bin9 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mHad_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 163 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 164 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 165 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin10 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 166 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin11 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 167 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin12 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 168 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin13 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 169 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin14 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 170 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin15 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 171 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin16 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 172 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin17 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 173 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin18 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 174 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin19 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin19 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 175 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 176 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin20 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin20 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 177 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin21 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin21 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 178 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin22 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin22 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 179 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin23 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin23 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 180 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 181 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 182 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 183 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 184 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 185 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 186 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pm_mVis_bin9 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pm_mVis_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 187 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin0 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 188 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin1 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 189 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin10 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 190 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin11 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 191 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin12 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 192 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin13 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 193 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin14 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 194 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin2 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 195 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin3 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 196 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin4 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 197 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin5 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 198 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin6 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 199 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin7 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 200 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin8 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 201 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binAddTk_pp_mHad_bin9 --algo impact --redefineSignalPOIs r -P prop_binAddTk_pp_mHad_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 202 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin0 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 203 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin1 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 204 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin10 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 205 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin11 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 206 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin12 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 207 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin13 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 208 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin14 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 209 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin15 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 210 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin16 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 211 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin2 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 212 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin3 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 213 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin4 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 214 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin5 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 215 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin6 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 216 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin7 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 217 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin8 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 218 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin0_bin9 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin0_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 219 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin0 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 220 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin1 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 221 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin10 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 222 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin11 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 223 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin12 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 224 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin13 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin13 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 225 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin14 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin14 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 226 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin15 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin15 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 227 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin16 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 228 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin17 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 229 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin18 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 230 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin2 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 231 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin3 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 232 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin4 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 233 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin5 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 234 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin6 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 235 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin7 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 236 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin8 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 237 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binEst_mu_q2bin1_bin9 --algo impact --redefineSignalPOIs r -P prop_binEst_mu_q2bin1_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 238 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin0 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 239 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin1 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 240 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin2 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 241 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin3 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 242 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin4 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 243 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin5 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 244 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin6 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 245 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin0_bin7 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin0_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 246 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin0 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 247 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin1 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 248 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin10 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin10 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 249 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin11 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin11 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 250 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin12 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin12 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 251 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin2 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 252 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin3 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin3 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 253 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin4 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin4 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 254 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin5 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin5 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 255 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin6 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin6 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 256 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin7 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin7 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 257 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin8 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin8 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 258 ]; then
  combine -M MultiDimFit -n _paramFit_v7_prop_binM2_miss_q2bin1_bin9 --algo impact --redefineSignalPOIs r -P prop_binM2_miss_q2bin1_bin9 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 259 ]; then
  combine -M MultiDimFit -n _paramFit_v7_trgSF --algo impact --redefineSignalPOIs r -P trgSF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi
if [ $1 -eq 260 ]; then
  combine -M MultiDimFit -n _paramFit_v7_xsecpp2bbXlumi --algo impact --redefineSignalPOIs r -P xsecpp2bbXlumi --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --X-rtd MINIMIZER_analytic -D data_obs --verbose -1 -m 120 -d cards/v7.root
fi