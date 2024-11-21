#!/bin/bash

# if datacard folder already exists, change the folder name

# setting up environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits
eval `scramv1 runtime -sh`

process_id=$1

pdfs=("main" "ua2mod")
test_folder=("bias_genmainfitua2mod_4pars")

mMed_values=(200 250 300 350 400 450 500 550)
exp_sig=(0. 0. 0. 0. 0. 0. 0. 0.)

for ((s=0; s<${#mMed_values[@]}; s+=1))
  do
    echo "start with DC: mMed = ${mMed_values[s]}"
    for ((j=0; j<${#pdfs[@]}; j+=2))
    do
      echo "start with bias tests"
      rseed=$((1000+s)) # to produce different random seed for each job      
      #gentoys
      python3 cli_boosted.py gentoys dc_20241114_cutbased/dc_SVJ_s-channel_mMed-${mMed_values[s]}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --outfile $process_id -t 30 --expectSignal ${exp_sig[s]} --rMin -2 --rMax 2 --pdf ${pdfs[j]} --outdir toysgen_mz${mMed_values[s]}_gen${pdfs[j]}_siginj0p0
      #fittoys
      python3 cli_boosted.py fittoys dc_20241114_cutbased/dc_SVJ_s-channel_mMed-${mMed_values[s]}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --toysFile toysgen_mz${mMed_values[s]}_gen${pdfs[j]}_siginj0p0/higgsCombineObserveddc_SVJ_s-channel_mMed-${mMed_values[s]}_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth${process_id}.GenerateOnly.mH120.123456.root --expectSignal ${exp_sig[s]} --rMin -2 --rMax 2 --pdf ${pdfs[j+1]} --cminDefaultMinimizerStrategy=0 --outfile $process_id --outdir toysfit_mz${mMed_values[s]}_gen${pdfs[j]}fit${pdfs[j+1]}_siginj0p0 #+ -s $rseed
      
      echo "done with bias tests!!!"
      i=$((i+1))
    done
done

mkdir ${test_folder[j/2]}_siginj0p0
cp toysgen_mz${mMed_values[s]}_gen${pdfs[j]}_siginj0p0 ${test_folder[j/2]}_siginj0p0
cp toysfit_mz${mMed_values[s]}_gen${pdfs[j]}fit${pdfs[j+1]}_siginj0p0 ${test_folder[j/2]}_siginj0p0
