#!/bin/bash

# if datacard folder already exists, change the folder name

# setting up environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits
eval `scramv1 runtime -sh`

process_id=$1

pdfs=("ua2" "ua2" "main" "ua2")
test_folder=("selftest" "bias_genmainfitua2" "bias_genaltfitua2")

mMed_values=(200 250 300 350 400 450 500 550)

for mMed in "${mMed_values[@]}"
  do
    echo "start with DC: mMed = $mMed"
    #python3 cli_boosted.py gen_datacards --bkg cutbased/merged_20240920/bkg_sel-cutbased.json --sig cutbased/smooth_20240920/SVJ_s-channel_mMed-${mMed}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json --outdir dc_20241028_cutbased
    for ((j=0; j<${#pdfs[@]}; j+=2))
    do
      echo "start with bias tests"
      
      #gentoys
      python3 cli_boosted.py gentoys dc_20241028_cutbased/dc_SVJ_s-channel_mMed-${mMed}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --outfile $process_id -t 30 --expectSignal 0.2 --rMin -2 --rMax 2 --pdf ${pdfs[j]} --outdir toysgen_mz${mMed}_gen${pdfs[j]}_siginj0p2
      #fittoys
      python3 cli_boosted.py fittoys dc_20241028_cutbased/dc_SVJ_s-channel_mMed-${mMed}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --toysFile toysgen_mz${mMed}_gen${pdfs[j]}_siginj0p2/higgsCombineObserveddc_SVJ_s-channel_mMed-${mMed}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth${process_id}.GenerateOnly.mH120.123456.root --expectSignal 0.2 --rMin -2 --rMax 2 --pdf ${pdfs[j+1]} --cminDefaultMinimizerStrategy=0 --outfile $process_id --outdir toysfit_mz${mMed}_gen${pdfs[j]}fit${pdfs[j+1]}_siginj0p2
      
      echo "done with bias tests!!!"
      mkdir ${test_folder[j/2]}_siginj0p2 
      mv toysgen_mz${mMed}_gen${pdfs[j]}_siginj0p2 ${test_folder[j/2]}_siginj0p2
      mv toysfit_mz${mMed}_gen${pdfs[j]}fit${pdfs[j+1]}_siginj0p2 ${test_folder[j/2]}_siginj0p2
      i=$((i+1))
    done
done
