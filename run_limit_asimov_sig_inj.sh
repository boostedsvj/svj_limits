#!/bin/bash
# To run fully: ./run_bias_study.sh 
# To run with a specfic date: ./run_bias_study.sh -d 20240930 -f 

# Default values
datacards_date="20240919"  # Manually defined date for generating datacards
dc_date=$(date +%Y%m%d)    # Dynamically set today's date for the gentoys command
scan_date=$(date +%Y%m%d)  # Dynamically set today's date for the fittoys command
mDark_value="10"
#run_only_fits=false        # Default to running all loops
run_only_fits=false        # Default to running all loops

# Parse command-line arguments
while getopts "d:f" opt; do
  case $opt in
    d)
      scan_date="$OPTARG"   # User-defined date for fittoys
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

# Outer loop: Generate the datacards (skip if only running fits)
# Array of mMed values
mMed_values=(200 250 300 350 400 450 500 550)

# Outer loop: Generate the datacards (skip if only running fits)
if [ "$run_only_fits" = false ]; then
  for mMed in "${mMed_values[@]}"
  do
    # Generate datacards for the current mMed value with variable mDark and datacards_date
    python3 cli_boosted.py gen_datacards --bkg hists/merged_${datacards_date}/bkg_sel-bdt\=0.65.json --sig hists/smooth_${datacards_date}/SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.json
  done

  # Run the 'gentoys' command with the current mMed value and variable mDark
  python3 cli_boosted.py gentoys dc_${dc_date}_bdt\=0.65/dc_SVJ_s-channel_mMed-350_mDark-${mDark_value}_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt\=0.65_smooth.txt -t -1 --expectSignal 0.2 -s 1001

  # Likelihood scan for expected limits
  python3 cli_boosted.py likelihood_scan_mp dc_${dc_date}_bdt\=0.65/dc*mDark-10_rinv-0p3*bdt*smooth.txt --rMin 0.0 --rMax 5.0 --seed 1001 --asimov

fi

# Run the multiple file likelihood scan for an r range of 0 to 10
python3 cli_boosted.py likelihood_scan_mp dc_${dc_date}_bdt\=0.65/dc*mDark-10_rinv-0p3*bdt*smooth.txt --rMin 0.0 --rMax 5.0 --seed 1001 -t -1 --toysFile toys_${dc_date}/higgsCombineObserveddc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt=0.65_smooth.GenerateOnly.mH120.1001.root

# plot
python3 quick_plot.py brazil scans_sig_inj_test/higgsCombine*rinv-0p3*.root -o asimov_350sig_test_0p3_mdark10.pdf


