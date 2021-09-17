#!/usr/bin/env bash

declare -a processes=(
    #"CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToDstarTauNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BsToMuNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BsToTauNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BuToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #"CP_BsToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
)

output=$HOME/BPhysics/data/cmsMC

for process in "${processes[@]}"; do
    echo $process
    input_dir=$output/$process/ntuples/B2DstMu
    output_dir=$output/$process/skimmed/B2DstMu
    mkdir -p $output_dir
    ./skimmer.py $input_dir/out*.root --outdir $output_dir --applyCorr --selection B2DstMu
    sleep 1
done

./skimmer.py ~/BPhysics/data/cmsMC/CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/ntuples/B2JpsiKst/out*.root --outdir ~/BPhysics/data/cmsMC/CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/skimmed/B2JpsiKst --applyCorr --selection B2JpsiKst

#./skimmer.py ~/BPhysics/data/cmsMC/CP_General_BdToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen/ntuples/B2JpsiK/out*.root --outdir ~/BPhysics/data/cmsMC/CP_General_BdToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen/skimmed/B2JpsiK --applyCorr --selection B2JpsiKst
