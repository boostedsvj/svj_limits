# Limits for the SVJ boosted analysis

## Setup 

1. Follow the `combine` instructions: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation .
Current results are using release `CMSSW_11_3_4`, tag v9.1.0.
Make sure to install `CombineHarvester` as well (bottom of that page).

2. Clone this repository and install dependencies:

```bash
cd $CMSSW_BASE/src
cmsenv
scram-venv
cmsenv
pip3 install git+https://github.com/boostedsvj/seutils
pip3 install git+https://github.com/boostedsvj/svj_ntuple_processing
git clone git@github.com:boostedsvj/svj_limits.git boosted/svj_limits
cd boosted/svj_limits
```

The code currently assumes Python 3.
For convenience, you can do:

```bash
alias python=python3
```

The commands below assume you are using this alias; if you're not, replace `python` with `python3`.


## Generating the datacards

You first need json files for signal, background, and (optionally) data; see https://github.com/boostedsvj/svj_uboost.

Then:

```bash
python cli_boosted.py gen_datacards --bkg merged_20240729/bkg_sel-cutbased.json --sig smooth_20240729/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json
for SIGNAL in cutbased/smooth_20240920/SVJ_s-channel_mMed-*_mDark-10_rinv-0p3*.json; do python cli_boosted.py gen_datacards --bkg cutbased/merged_20240920/bkg_sel-cutbased.json --sig $SIGNAL; done
```


## Running the likelihood scans

For all BDT working points and all signals, do simply:

```bash
python cli_boosted.py likelihood_scan dc_Dec07/*.txt
python cli_boosted.py likelihood_scan dc_Dec07/*.txt --asimov
```

Selecting BDT working points and particular signals is easily done via wildcard patterns to select the right datacards, e.g.:

```bash
python cli_boosted.py likelihood_scan dc_Dec07_minmt300/dc_mz*rinv0.3*bdt0p{0,3,5}*.txt --asimov --minmu -.5 --maxmu .5 -n 100
```

Note also the options `--minmu` and `--maxmu` which handle the range of the signal parameter to scan, and the option `-n` which controls the number of points in the range.


## Bias study

Generate toys:

```
python cli_boosted.py gentoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt -t 300 --expectSignal 0 -s 1001 --rMin -2 --rMax 3
for DC in dc_20240920_cutbased/*.txt; do python cli_boosted.py gentoys $DC -t 300 --expectSignal 0 -s 1001 --rMin -2 --rMax 3; done
```

Fit the toys:

```
python cli_boosted.py fittoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --toysFile toys_20240924/higgsCombineObserveddc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.GenerateOnly.mH120.1001.root --expectSignal 0 --rMin -2 --rMax 3
 
for MZ in 200 250 300 350 400 450 500 550; do python cli_boosted.py fittoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-$MZ*.txt --toysFile toys_20240924/higgsCombineObserveddc_SVJ_s-channel_mMed-$MZ*.root --expectSignal 0 --rMin -2 --rMax 3; done
```

## Nuisance impacts

The nuisance impacts show how much each systematic uncertainty impacts the fit, and whether any systematics are pulled or constrained.
The input is a datacard txt file.
Signal injection is required to see the effects of signal systematics (the uncertainty variations are checked at the likelihood minimum, which is `r=0` without signal injection).

```bash
python3 cli_boosted.py impacts dc.txt --nfits 16 --asimov --normRange 0.1 2.0 --rMin -10 --rMax 10 --robustFit 1 --expectSignal 0.2
```

## Plotting
Bkg fit 

```bash
python3 quick_plot.py bkgfit ua2 --bkg cutbased/merged_20240920/bkg_sel-cutbased.json --sig cutbased/smooth_20240920/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json --outfile fit_bkg.png
```

ΔNNL as a function of mu:

```bash
python quick_plot.py muscan scans_Dec07/*bdt0p3*Scan*.root
```

Fit to background distribution:
```bash
python3 quick_plot.py bkgfit ua2 --bkg merged_20240729/bkg_sel-cutbased.json --sig smooth_20240729/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json --outfile fit_bkg.png
```

MT histogram, with bkg-only fit and and sig+bkg fit:

```bash
python quick_plot.py mtdist scans_Dec07/higgsCombineObserved_dc_mz450_rinv0.3_bdt0p300Bestfit.MultiDimFit.mH120.root
```

Note you should use the `Bestfit`-tagged file, not `Scan`.
Apparently, the single snapshot stored in the `Scan` files is _not_ the best fit.


_Warning: Below here, readme outdated; need to check_

CLS:

```bash
python quick_plot.py cls scans_Dec07/higgsCombineObserved_dc_mz450_rinv0.3_bdt0p300.MultiDimFit.mH120.root scans_Dec07/higgsCombineAsimov_dc_mz450_rinv0.3_bdt0p300.MultiDimFit.mH120.root
```


Brazil band (relies on good interpolation; always check the CLs plots to double check!):

```bash
python quick_plot.py brazil scans_Dec07/higgsCombine*bdt0p3*.root
```
