#!/bin/bash

# Get the current date in the format YYYYMMDD
hists_date=20241115
dc_date=$(date +%Y%m%d)  
scan_date=$(date +%Y%m%d)  
sel="bdt=0.67"
old_datacards=false
skip_likelihood_scan=false
mMed_values=(200 250 300 350 400 450 500 550)
rinv_values=(0 0p1 0p2 0p3 0p4 0p5 0p6 0p7 0p8 0p9 1)
mDark_values=(1 5 10)

# Optionally, you can specify a different date as an argument
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) hists_date="$2"; shift ;; 
        --datacard_date) datacard_date="$2"; shift ;; 
        --sel) sel="$2"; shift ;;
        --old_datacards) old_datacards=true ;; 
        --skip_likelihood_scan) skip_likelihood_scan=true ;; 
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Generate the datacards
if [ "$old_datacards" = false ]; then
  for mDark_value in "${mDark_values[@]}"; do
    for rinv_value in "${rinv_values[@]}"; do
      for mMed in "${mMed_values[@]}"; do
        # Generate datacards for the current mMed value with variable mDark and hists_date
        python3 cli_boosted.py gen_datacards \
          --bkg hists/merged_${hists_date}/bkg_sel-${sel}.json \
          --sig hists/smooth_${hists_date}/SVJ_s-channel_mMed-${mMed}_mDark-${mDark_value}_rinv-${rinv_value}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-${sel}_smooth.json
      done
    done
  done
fi

# Likelihood scan for expected limits
# These become the 'asimov' files
if [ "$skip_likelihood_scan" = false ]; then
  python3 cli_boosted.py likelihood_scan_mp \
    dc_${dc_date}_${sel}/dc*mDark-*_rinv-*smooth.txt \
    --rMin 0.0 \
    --rMax 5.0 \
    --seed 1001 
    --asimov
fi

# Loop through all combinations of mz and rinv
for mDark in "${mDark_values[@]}"; do
  for mz in "${mMed_values[@]}"; do
    for rinv in "${rinv_values[@]}"; do
      # Replace 'p' with '.' for rinv when passing to the Python script
      rinv_float=$(echo "$rinv" | sed 's/p/./')
 
      # Check for observed and asimov files for the current mz and rinv
      asimov_file=$(ls scans_${scan_date}/higgsCombineAsimovdc_SVJ_s-channel_mMed-${mz}_mDark-${mDark}_rinv-${rinv}_alpha-peak_*Scan*.root 2>/dev/null)
      echo $(ls scans_${scan_date}/higgsCombineAsimovdc_SVJ_s-channel_mMed-${mz}_mDark-${mDark}_rinv-${rinv}_alpha-peak_*Scan*.root 2>/dev/null)
 
      # If both files exist, run the process_limits.py Python script
      if [[ -f "$asimov_file" ]]; then
        echo "Processing limits for mz=${mz}, rinv=${rinv_float}, mDark=${mDark} with observed and asimov files."
        python3 process_limits.py --mz "$mz" --rinv "$rinv_float" --mDark "$mDark" --date "$scan_date"
 
        # Run sigacceptance.py for the current mz and rinv combination
        #echo "Processing signal acceptance for mz=${mz}, rinv=${rinv_float}."
        #python3 sigacceptance.py --mz "$mz" --rinv "$rinv_float" 
          #--date "$scan_date"
      else
        echo "Skipping limits processing for mz=${mz}, rinv=${rinv_float}. Missing observed or asimov file in scans_${scan_date}."
      fi
    done
  done
done


