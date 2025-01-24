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
hists_date="20241115"   # Date of the histograms used to make the data cards
dc_date=$(date +%Y%m%d)     # Today's date for gentoys command
toys_date=$(date +%Y%m%d)   # Today's date for fittoys command
toyfits_date=$(date +%b%d)  # Date format for toyfits directory as "Oct29"
toy_seed=1001
mDark_value="10"
rinv_value="0p3"
run_only_fits=false         # Default to generate datacards, toys, and fit
skip_dc=false              # separate option to do just toys and fit
sel="bdt=0.67"            # Default selection type
test_type="self"            # Default test type
siginj="exp"                  # Default signal injected at exp limit strength
mMed_values=(200 250 300 350 400 450 500 550)  # Default mMed values
rmax=5

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) toys_date="$2"; shift ;; 
        -f) run_only_fits=true ;; # option to only run toy fits on existing toys
        --skip_dc) skip_dc=true ;;
        --sel) sel="$2"; shift ;;
        --test_type) test_type="$2"; shift ;;
        --dc_date) dc_date="$2"; shift ;;
        --hists_date) hists_date="$2"; shift ;;
        --siginj) siginj="$2"; shift ;;
        --rmax) rmax="$2"; shift ;;
        --mDark) mDark_value="$2"; shift ;;
        --rinv) rinv_value="$2"; shift ;; 
        --mMed_values) IFS=' ' read -r -a mMed_values <<< "$2"; shift ;;  # Parse mMed_values as array
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

get_signame(){
sig_name=SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth
}

# Set PDF option based on test_type
if [ "$test_type" == "bias" ]; then
    pdf_option="main"
else
    pdf_option="ua2"
fi

# This is used for naming the 'siginj' directory: 0 means no signal, 1 means signal injected (not sig strength)
if [ "$siginj" == 0 ]; then
    inj_dir=0
else
    inj_dir=1
fi

# Injected signal stregnth
declare -A sig_strength # declare associative array (bash >= 4.0)
if [ "$siginj" == "exp" ]; then
    sig_strength=( [200]=0.267 [250]=0.129 [300]=0.160 [350]=0.184 [400]=0.208 [450]=0.248 [500]=0.262 [550]=0.396 )
    # Cut-based values are needed!!
else
    sig_strength=( [200]=$siginj [250]=$siginj [300]=$siginj [350]=$siginj [400]=$siginj [450]=$siginj [500]=$siginj [550]=$siginj )
fi

# Generate the datacards (skip if only running fits)
if [ "$run_only_fits" == false ] && [ "$skip_dc" == false ]; then
  for mMed in "${mMed_values[@]}"
  do
    get_signame
    # Generate datacards for the current mMed value with variable mDark and hists_date
    (set -x; python3 cli_boosted.py gen_datacards \
      --bkg hists/merged_${hists_date}/bkg_sel-${sel}.json \
      --sig hists/smooth_${hists_date}/${sig_name}.json)
  done
fi

if [ "$run_only_fits" == false ]; then
  # First loop: Run the 'gentoys' command
  for mMed in "${mMed_values[@]}"
  do
    get_signame
	expect_sig=0
    if [ "$siginj" != 0 ]; then
	  expect_sig=${sig_strength[$mMed]}
	fi
    # Run the 'gentoys' command with the current mMed value and variable mDark
    (set -x; python3 cli_boosted.py gentoys \
      dc_${dc_date}_${sel}/dc_${sig_name}.txt \
      -t 300 \
      --expectSignal ${expect_sig} \
      -s ${toy_seed} \
      --pdf ${pdf_option})
  done
fi

# Fit the toys
for mMed in "${mMed_values[@]}"
do
  get_signame
  # Run the 'fittoys' command with the current mMed value and variable mDark
  # --range is used for likelihood fits and rMin/rMax for toy fits
  (set -x; python3 cli_boosted.py fittoys \
    dc_${dc_date}_${sel}/dc_${sig_name}.txt \
    --toysFile toys_${toys_date}/higgsCombineObserveddc_${sig_name}.GenerateOnly.mH120.${toy_seed}.root \
    --expectSignal 0 \
    --rMax ${rmax} \
    --rMin -${rmax} \
    --pdf ua2 \
    --cminDefaultMinimizerStrategy=0)
done

# Create the results directory that is expected by the plotting code
results_dir="${test_type}_test/siginj${inj_dir}"
mkdir -p "$results_dir"

# Copy the ROOT files from the generated toyfits directory to the new results directory
toyfits_dir="toyfits_${toyfits_date}"
mv "$toyfits_dir"/*.root "$results_dir/"
