import os
import boosted_fits as bsvj
import quick_plot as tools
import ROOT as r
import numpy as np
import glob
import re  # Import regex for extracting mDark from filenames

# Pull mz, rinv, and date from arguments
mz = bsvj.pull_arg('--mz', type=float, required=True).mz
mDark = bsvj.pull_arg('--mDark', type=float, required=True).mDark
rinv = bsvj.pull_arg('--rinv', type=float, required=True).rinv
date = bsvj.pull_arg('--date', type=str, required=True).date

# Directory where the ROOT files are stored based on the date
input_directory = f"scans_{date}"

# Convert rinv to the format used in file names (e.g., 0.9 -> 0p9)
if rinv > 0.0 and rinv < 1.0 : rinv_str = str(rinv).replace('.', 'p')
else : rinv_str = str(int(rinv))

# Find matching rootfiles for observed and asimov based on mz and rinv
#observed_files = glob.glob(f"{input_directory}/higgsCombineObserveddc_SVJ_s-channel_mMed-{int(mz)}_mDark-{int(mDark)}_rinv-{rinv_str}_alpha-peak_*Scan*.root")
observed_files = glob.glob(f"{input_directory}/higgsCombineAsimovdc_SVJ_s-channel_mMed-550_mDark-10_rinv-1_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt=0.67_*Scan*.MultiDimFit.mH120.1001.root")
asimov_files = glob.glob(f"{input_directory}/higgsCombineAsimovdc_SVJ_s-channel_mMed-{int(mz)}_mDark-{int(mDark)}_rinv-{rinv_str}_alpha-peak_*Scan*.root")

# Ensure we have matching files for both observed and asimov
if len(observed_files) == 0:
    raise ValueError(f"No observed ROOT files found for mz={mz} and rinv={rinv_str} in {input_directory}")
if len(asimov_files) == 0:
    raise ValueError(f"No asimov ROOT files found for mz={mz} and rinv={rinv_str} in {input_directory}")

clean = bsvj.pull_arg('--clean', action='store_true').clean

# Create output directory based on the date if it doesn't exist
output_directory = f"cls_lim_{date}"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Make a new TTree (limit) in a new ROOT file, write CLS into it
output_file = f"{output_directory}/Limit_mz{int(mz)}_rinv{rinv}.root"
file = r.TFile.Open(output_file, "RECREATE")
tree = r.TTree("limit", "limit")

# Define the quantities for branches, now with only "limit" and "quantileExpected"
qtys = ["limit", "quantileExpected"] + ["trackedParam_{}".format(q) for q in ["mZprime", "mDark", "rinv", "xsec"]]
#keys = ["expected", "twosigma_down", "onesigma_down", "onesigma_up", "twosigma_up"]

# Define struct for ROOT
r.gROOT.ProcessLine("struct quantile_t { " + " ".join([f"Double_t {qty};" for qty in qtys]) + " };")
qobj = r.quantile_t()

# Create branches for each quantity
for qty in qtys:
    tree.Branch(qty, r.addressof(qobj, qty), f'{qty}/D')

# Get xsec from boosted_fits
xsec = bsvj.get_xs(mz)

# Process the observed and asimov files together
for observed, asimov in zip(observed_files, asimov_files):
    # Extract mDark from the observed filename using regex
    mDark_match = re.search(r'mDark-(\d+)', observed)
    if mDark_match:
        mDark = float(mDark_match.group(1))
    else:
        print(f"Warning: mDark not found in {observed}. Skipping this entry.")
        continue  # Skip this entry if mDark is not found

    obs, asimov = tools.extract_scans([observed, asimov], correct_minimum=True)
    if clean:
        obs = tools.clean_scan(obs)
        asimov = tools.clean_scan(asimov)
    
    cls = tools.get_cls(obs, asimov)
limit = tools.interpolate_95cl_limit(cls)


f = r.TFile(observed_files[0], 'READ')
t = f.Get("limit")
#N = t.GetEntries()

#for i in range(len(limit.observed)):
#    t.GetEntry(i)
qobj.trackedParam_mZprime = mz  # Get from input
qobj.trackedParam_mDark = mDark  # Extracted from filename
qobj.trackedParam_rinv = rinv  # Get from input
qobj.trackedParam_xsec = xsec  # From boosted_fits

# Loop over limits and corresponding quantile values
for lim, q in zip([limit.observed, limit.twosigma_down, limit.onesigma_down, limit.expected, limit.onesigma_up, limit.twosigma_up],
                  [-1, 0.025, 0.16, 0.50, 0.84, 0.975]):
    print("New Entries")
    print(lim)
    print(q)
    qobj.limit = lim
    qobj.quantileExpected = q  # Assign corresponding quantileExpected value
    tree.Fill()  # Fill the tree with current limit and quantile


file.Write()
file.Close()

