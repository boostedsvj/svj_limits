#!/bin/bash
#==============================================================================
# run_limit_asimov_sig_inj.sh --------------------------------------------------
#------------------------------------------------------------------------------
# Author(s): Brendan Regnery --------------------------------------------------
#------------------------------------------------------------------------------
# Basic functionality:
#   Creates limits with expected 'asimov' values with an observed limit of an
#      asimov toy with signal injected at 350 GeV (same toy for all mass points)
#------------------------------------------------------------------------------
# To run fully: ./run_limit_asimov_sig_inj.sh 
# To run with a specfic date: ./run_limit_asimov_sig_inj.sh --mMed_values "300 350"

# Default values
hists_dir="hists"
hists_date="20241115"  # Date of the histograms used for making datacards
hists_date_anti=  # Date of the anti-tag CR histograms
dc_date=$(date +%Y%m%d)    # Dynamically set today's date 
scan_date=$(date +%Y%m%d)
toy_seed=1001
sel="bdt=0.67"
siginj=0.2
mInj=350
mDark_value="10"
rinv_value="0p3"
mMed_values=(200 250 300 350 400 450 500 550)
run_only_fits=false        # Default to running all loops
skip_dc=false              # separate option to do just toys and fit
only_inj=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) toys_date="$2"; shift ;; 
        -f) run_only_fits=true ;; # option to only run likelihood scan
        --skip_dc) skip_dc=true ;;
        --only_inj) only_inj=true ;; # only generates and fits 'observed' signal injected asimov toy
        --sel) sel="$2"; shift ;;
        --hists_dir) hists_dir="$2"; shift ;;
        --hists_date) hists_date="$2"; shift ;;
        --hists_date_anti) hists_date_anti="$2"; shift ;;
        --dc_date) dc_date="$2"; shift ;;
        --siginj) siginj="$2"; shift ;;
        --mInj) mInj="$2"; shift ;;
        --mDark) mDark_value="$2"; shift ;;
        --rinv) rinv_value="$2"; shift ;; 
        --mMed_values) IFS=' ' read -r -a mMed_values <<< "$2"; shift ;;  # Parse mMed_values as array
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -z "$hists_date_anti" ]; then
    hists_date_anti=hists_date
fi

get_signame(){
MASS=$1
if [ -z "$MASS" ]; then MASS=${mMed}; fi
sig_name=SVJ_s-channel_mMed-${MASS}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth
sig_name_anti=SVJ_s-channel_mMed-${MASS}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-anti${sel}_smooth
}

# if you want to generate cards and the asimov toy
if [ "$run_only_fits" == false ] && [ "$skip_dc" == false ]; then

  # Generate the datacards
  for mMed in "${mMed_values[@]}"
  do
    get_signame
    # Generate datacards for the current mMed value with variable mDark and hists_date
    (set -x; python3 cli_boosted.py gen_datacards \
      --regions ${sel} anti${sel} \
      --norm-type crtf \
      --bkg ${hists_dir}/merged_${hists_date}/bkg_sel-${sel}.json ${hists_dir}/merged_${hists_date_anti}/bkg_sel-anti${sel}.json \
      --sig ${hists_dir}/smooth_${hists_date}/${sig_name}.json ${hists_dir}/smooth_${hists_date_anti}/${sig_name_anti}.json)
  done
fi

if [ "$run_only_fits" == false ]; then
  # Generate an asimov toy with a signal at requested Zprime mass
  get_signame ${mInj}
  (set -x; python3 cli_boosted.py gentoys \
    dc_${dc_date}_${sel}/dc_${sig_name}.txt \
    -t -1 \
    --expectSignal ${siginj} \
    -s ${toy_seed})

  # Likelihood scan for expected limits
  # These become the 'asimov' files
  if [ "$only_inj" = false ]; then
    (set -x; python3 cli_boosted.py likelihood_scan_mp \
      dc_${dc_date}_${sel}/dc*mDark-${mDark_value}_rinv-${rinv_value}*${sel}*smooth.txt \
      --range 0.0 2.0 \
      --seed ${toy_seed} \
      --asimov)
  fi

fi

# Run the multiple file likelihood scan on the asimov toy with injected signal
# These become the 'observed' files
get_signame ${mInj}
(set -x; python3 cli_boosted.py likelihood_scan_mp \
  dc_${dc_date}_${sel}/dc*mDark-${mDark_value}_rinv-${rinv_value}*${sel}*smooth.txt \
  --range 0.0 2.0 \
  --seed ${toy_seed} \
  --toysFile toys_${dc_date}/higgsCombineObserveddc_${sig_name}.GenerateOnly.mH120.${toy_seed}.root \
  -t -1)

# plot
(set -x; python3 quick_plot.py brazil \
  scans_${dc_date}/higgsCombine*rinv-${rinv_value}*.root \
  -o asimov_${mInj}_sig_test_${rinv_value}_mdark${mDark_value}.pdf)
