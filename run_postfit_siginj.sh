#!/bin/bash
#==============================================================================
# run_postfit_sig_inj.sh --------------------------------------------------
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
dc_date=$(date +%Y%m%d)    # Dynamically set today's date 
scan_date=$(date +%Y%m%d)
toy_seed=123456
sel="bdt=0.67"
mDark_value="10"
rinv_value="0p3"
mMed_values=(200 350 550)
skip_dc=false              # separate option to do just toys and fit
only_inj=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) toys_date="$2"; shift ;; 
        --skip_dc) skip_dc=true ;;
        --only_inj) only_inj=true ;; # only generates and fits 'observed' signal injected asimov toy
        --sel) sel="$2"; shift ;;
        --hists_dir) hists_dir="$2"; shift ;;
        --hists_date) hists_date="$2"; shift ;;
        --dc_date) dc_date="$2"; shift ;;
        --mDark) mDark_value="$2"; shift ;;
        --rinv) rinv_value="$2"; shift ;; 
        --mMed_values) IFS=' ' read -r -a mMed_values <<< "$2"; shift ;;  # Parse mMed_values as array
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

get_signame(){
MASS=$1
if [ -z "$MASS" ]; then MASS=${mMed}; fi
sig_name=SVJ_s-channel_mMed-${MASS}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth
}

sig_strength=( [200]=0.267 [250]=0.129 [300]=0.160 [350]=0.184 [400]=0.208 [450]=0.248 [500]=0.262 [550]=0.396 )

# if you want to generate cards and the asimov toy
if [ "$skip_dc" == false ]; then

  # Generate the datacards
  for mMed in "${mMed_values[@]}"
  do

    mInj=${mMed}
    siginj=${sig_strength[$mMed]}

    get_signame
    # Generate datacards for the current mMed value with variable mDark and hists_date
    (set -x; python3 cli_boosted.py gen_datacards \
      --bkg ${hists_dir}/merged_${hists_date}/bkg_sel-${sel}.json \
      --sig ${hists_dir}/smooth_${hists_date}/${sig_name}.json)
    
    # Generate an asimov toy with a signal at requested Zprime mass
    get_signame ${mInj}
    (set -x; python3 cli_boosted.py gentoys \
      dc_${dc_date}_${sel}/dc_${sig_name}.txt \
      -t -1 \
      --expectSignal ${siginj} \
      -s ${toy_seed})

    # Fit the Asimov toy
    (set -x; python3 cli_boosted.py bestfit \
      dc_${dc_date}_${sel}/dc_${sig_name}.txt \
      -t -1 --tIsToyIndex --saveToys \
      --toysFile toys_${dc_date}/higgsCombineObserveddc_${sig_name}.GenerateOnly.mH120.${toy_seed}.root )

    # Plot the bestfit
    #bestfits_${dc_date}/higgsCombineObservedBestfit_dc_${sig_name}.MultiDimFit.mH120.root \
    (set -x; python3 quick_plot.py mtdist \
      bestfits_${dc_date}/higgsCombineObservedToy-1Bestfit_dc_${sig_name}.MultiDimFit.mH120.${toy_seed}.root \
      --outfile postfit_${mInj}.pdf )

    # Plot the bestfit
    #bestfits_${dc_date}/higgsCombineObservedBestfit_dc_${sig_name}.MultiDimFit.mH120.root \
    (set -x; python3 quick_plot.py mtdist \
      bestfits_${dc_date}/higgsCombineObservedToy-1Bestfit_dc_${sig_name}.MultiDimFit.mH120.${toy_seed}.root \
      --outfile postfit_${mInj}.png )

  done
fi

