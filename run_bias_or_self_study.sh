#!/bin/bash
#==============================================================================
# run_bias_or_self_study.sh ---------------------------------------------------
#------------------------------------------------------------------------------
# Author(s): Brendan Regnery --------------------------------------------------
#------------------------------------------------------------------------------
# Basic functionality:
#   run a bias (main/alt toy gen + ua2 fit) or self test (ua2 gen + ua2 fit) 
#   default is a self test with no signal injected with options to change 
#------------------------------------------------------------------------------
# To run with default values: ./run_bias_or_self_study.sh 
# To run with specific options: 
# ./run_bias_study.sh -d 20240930 -f --sel cutbased --test_type bias --siginj 1.0 --mMed_values "200 300 500"

# Default values
hists_date="20240919"   # Date of the histograms used to make the data cards
dc_date=$(date +%Y%m%d)     # Today's date for gentoys command
toys_date=$(date +%Y%m%d)   # Today's date for fittoys command
toyfits_date=$(date +%b%d)  # Date format for toyfits directory as "Oct29"
mDark_value="10"
rinv_value="0p3"
run_only_fits=false         # Default to generate datacards, toys, and fit
sel="bdt=0.65"            # Default selection type
test_type="self"            # Default test type
siginj=0.0                  # Default injected signal value
mMed_values=(200 250 300 350 400 450 500 550)  # Default mMed values

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) toys_date="$2"; shift ;; 
        -f) run_only_fits=true ;; # option to only run toy fits on existing toys
        --sel) sel="$2"; shift ;;
        --test_type) test_type="$2"; shift ;;
        --hists_date) hists_date="$2"; shift ;;
        --siginj) siginj="$2"; shift ;;
        --mDark) mDark_value="$2"; shift ;;
        --rinv) rinv_value="$2"; shift ;; 
        --mMed_values) IFS=' ' read -r -a mMed_values <<< "$2"; shift ;;  # Parse mMed_values as array
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Set PDF option based on test_type
if [ "$test_type" == "bias" ]; then
    pdf_option="main"
else
    pdf_option="ua2"
fi

# This is used for the 'siginj' directory 0 means no signal, 1 means signal injected (not sig strength)
if (( $(echo "$siginj == 0" | bc -l) )); then
    inj=0
else
    inj=1
fi

# Generate the datacards (skip if only running fits)
if [ "$run_only_fits" = false ]; then
  for mMed in "${mMed_values[@]}"
  do
    # Generate datacards for the current mMed value with variable mDark and hists_date
    python3 cli_boosted.py gen_datacards \
      --bkg hists/merged_${hists_date}/bkg_sel-${sel}.json \
      --sig hists/smooth_${hists_date}/SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth.json
  done

  # First loop: Run the 'gentoys' command
  for mMed in "${mMed_values[@]}"
  do
    # Run the 'gentoys' command with the current mMed value and variable mDark
    python3 cli_boosted.py gentoys \
      dc_${dc_date}_${sel}/dc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth.txt \
      -t 300 \
      --expectSignal $siginj \
      -s 1001 \
      --pdf $pdf_option
  done
fi

# Fit the toys
for mMed in "${mMed_values[@]}"
do
  # Run the 'fittoys' command with the current mMed value and variable mDark
  python3 cli_boosted.py fittoys \
    dc_${dc_date}_${sel}/dc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth.txt \
    --toysFile toys_${toys_date}/higgsCombineObserveddc_SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth.GenerateOnly.mH120.1001.root \
    --expectSignal 0 \
    --rMax 5 \
    --rMin -5 \
    --pdf ua2 \
    --cminDefaultMinimizerStrategy=0
done

# Create the results directory that is expected by the plotting code
results_dir="${test_type}_test/siginj${inj}"
mkdir -p "$results_dir"

# Copy the ROOT files from the generated toyfits directory to the new results directory
toyfits_dir="toyfits_${toyfits_date}"
mv "$toyfits_dir"/*.root "$results_dir/"

