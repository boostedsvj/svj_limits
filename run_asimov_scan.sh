#!/bin/bash
# Default values
hists_dir="hists"
hists_date="20241106" # Date of the histograms used for making datacards
scan_date=$(date +%Y%m%d)
toy_seed=1001
sel="cutbased"
mDark_values=(1 5 10)
rinv_values=(0 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 1)
mMed_values=(200 250 300 350 400 450 500 550)

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
  --sel)
    sel="$2"
    shift
    ;;
  --hists_dir)
    hists_dir="$2"
    shift
    ;;
  --hists_date)
    hists_date="$2"
    shift
    ;;
  --mDark_values)
    IFS=' ' read -r -a mMed_values <<<"$2"
    shift
    ;; # Parse mMed_values as array
  --rinv_values)
    IFS=' ' read -r -a mMed_values <<<"$2"
    shift
    ;; # Parse mMed_values as array
  --mMed_values)
    IFS=' ' read -r -a mMed_values <<<"$2"
    shift
    ;; # Parse mMed_values as array
  *)
    echo "Unknown parameter passed: $1"
    exit 1
    ;;
  esac
  shift
done

get_signame() {
  local mMed=$1
  local mDark=$2
  local rinv=$3
  echo SVJ_s-channel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth
}

# Generate the datacards
for mMed in "${mMed_values[@]}"; do
  for mDark in "${mDark_values[@]}"; do
    for rinv in "${rinv_values[@]}"; do
      get_signame $mMed $mDark $rinv
      # Generate datacards for the current mMed value with variable mDark and hists_date
      (
        set -x
        python3 cli_boosted.py gen_datacards \
          --bkg ${hists_dir}/merged_${hists_date}/bkg_sel-${sel}.json \
          --sig ${hists_dir}/smooth_${hists_date}/$(get_signame $mMed $mDark $rinv).json
      )&
    done
    wait
  done
done

# Generating the asimov toy as this is background only, we dont need to do additional selection
get_signame ${mInj}
(
  set -x
  python3 cli_boosted.py gentoys \
    dc_${scran_date}_${sel}/dc_$(get_signame 350 10 0p3).txt \
    -t -1 \
    --expectSignal 0.0 \
    -s ${toy_seed}
)

function get_result() {
  rfile="scans_${scan_date}/higgsCombinedc_$(get_signame $1 $2 $3)ScanAsimov.MultiDimFit.mH120.${toy_seed}.root"
  if [[ -f $rfile ]]; then
    python3 -c "import uproot ; print(uproot.open('$rfile')['limit'].arrays()[b'deltaNLL'][-1] > 1.5);"
  else
    echo "False"
  fi
}

for rmax in 2 3; do
  required_dcs=""
  for mMed in "${mMed_values[@]}"; do
    for mDark in "${mDark_values[@]}"; do
      for rinv in "${rinv_values[@]}" ; do
        if [ $(get_result $mMed $mDark $rinv) = "True" ]; then
          echo "Got good result for mMed-${mMed} mDark-${mDark} rinv-${rinv}"
        else
          required_dcs="${required_dcs} dc_${scan_date}_${sel}/dc_$(get_signame $mMed $mDark $rinv).txt"
        fi
      done
    done
  done
  echo $required_dcs
  python3 cli_boosted.py likelihood_scan_mp ${required_dcs} \
    --range 0.0 ${rmax} \
    --seed ${toy_seed} \
    --asimov
done
