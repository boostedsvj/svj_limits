#!/bin/bash
# To run fully: ./run_bias_study.sh 
# To run with a specfic date: ./run_bias_study.sh -d 20240930 -f 

# Default values
datacards_date="20240919"  # Manually defined date for generating datacards
dc_date=$(date +%Y%m%d)    # Dynamically set today's date for the gentoys command
toys_date=$(date +%Y%m%d)  # Dynamically set today's date for the fittoys command
mDark_value="10"
run_only_fits=false        # Default to running all loops

# Parse command-line arguments
while getopts "d:f" opt; do
  case $opt in
    d)
      toys_date="$OPTARG"   # User-defined date for fittoys
      ;;
    f)
      run_only_fits=true    # Option to run only the fits loop
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Array of mMed values
mMed_values=(200 250 300 350 400 450 500 550)

# Outer loop: Generate the datacards (skip if only running fits)
if [ "$run_only_fits" = false ]; then
  for mMed in "${mMed_values[@]}"
  do
    # Generate datacards for the current mMed value with variable mDark and datacards_date
    python3 cli_boosted.py gen_datacards --bkg hists/merged_${datacards_date}/bkg_sel-bdt\=0.65.json --sig hists/smooth_${datacards_date}/SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.json -p main
  done

  # First loop: Run the 'gentoys' command
  for mMed in "${mMed_values[@]}"
  do
    # Run the 'gentoys' command with the current mMed value and variable mDark
    python3 cli_boosted.py gentoys dc_${dc_date}_bdt\=0.65/dc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.txt -t 300 --expectSignal 0 -s 1001 --pdf main
  done
fi

# Second loop: Run the 'fittoys' command (always run if run_only_fits is true)
for mMed in "${mMed_values[@]}"
do
  # Run the 'fittoys' command with the current mMed value and variable mDark
  python3 cli_boosted.py fittoys dc_${dc_date}_bdt\=0.65/dc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.txt --toysFile toys_${toys_date}/higgsCombineObserveddc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.GenerateOnly.mH120.1001.root --expectSignal 0 --rMax 5 --rMin -5 --pdf ua2 --cminDefaultMinimizerStrategy=0
done

