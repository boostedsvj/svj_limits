
#!/bin/bash
#==============================================================================
# all_sig_inj_asimov_limits.sh ------------------------------------------------
#------------------------------------------------------------------------------
# Author(s): Brendan Regnery --------------------------------------------------
#------------------------------------------------------------------------------
# Basic functionality:
#   Creates limits with expected 'asimov' values with an observed limit of an
#      asimov toy with signal injected at each z' mass
#      (same toy for all mass points)
#------------------------------------------------------------------------------
# To run : ./all_sig_inj_asimov_limits.sh 

mMed_values=(200 250 300 350 400 450 500 550)

# Injected signal stregnth
declare -A sig_strength # declare associative array (bash >= 4.0)
sig_strength=( [200]=0.271 [250]=0.136 [300]=0.165 [350]=0.189 [400]=0.210 [450]=0.249 [500]=0.265 [550]=0.397 )

for mMed in "${mMed_values[@]}"
do
  # run the full asimov expected limits only once
  if [ "$mMed" == 200 ]; then
    ./run_limit_asimov_sig_inj.sh --mInj ${mMed}  --siginj ${sig_strength[$mMed]} 
  # make only the new asimov toy with signal inject as the observed limit (and plot)
  else
    ./run_limit_asimov_sig_inj.sh --mInj ${mMed}  --siginj ${sig_strength[$mMed]}  --only_inj
  fi
done
