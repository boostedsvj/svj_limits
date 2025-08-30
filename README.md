# Limits for the SVJ boosted analysis

This repo contains all the setup required for limit setting for the BSVJ search. The first section describes setup and some possible analyses,
then the last section of this README contains commands that operate high-level scripts for preforming the necessary statistical tests.

## Setup

We are using a patched version of `combine` v9.2.0. The original instructions are here: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation .

```bash
cmssw-el7
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
scram-venv
cmsenv
git clone git@github.com:boostedsvj/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit -b seed_fit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester && cd CombineHarvester && git checkout 610d8ded8d1f71ed4db4d634cd9e2d7b8ae7b560 && cd -
scram b -j 4
pip3 install git+https://github.com/boostedsvj/seutils
pip3 install git+https://github.com/boostedsvj/svj_ntuple_processing
pip3 install git+https://github.com/boostedsvj/rhalphalib
git clone git@github.com:boostedsvj/svj_limits.git boosted/svj_limits
cd boosted/svj_limits
```

## Analytic fit

### Generating the datacards

You first need json files for signal, background, and (optionally) data; see https://github.com/boostedsvj/svj_uboost.

Then:

```bash
python3 cli_boosted.py gen_datacards --bkg merged_20240729/bkg_sel-cutbased.json --sig smooth_20240729/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json
for SIGNAL in cutbased/smooth_20240920/SVJ_s-channel_mMed-*_mDark-10_rinv-0p3*.json; do python3 cli_boosted.py gen_datacards --bkg cutbased/merged_20240920/bkg_sel-cutbased.json --sig $SIGNAL; done
```


### Running the likelihood scans

For all BDT working points and all signals, do simply:

```bash
python3 cli_boosted.py likelihood_scan dc_Dec07/*.txt
python3 cli_boosted.py likelihood_scan dc_Dec07/*.txt --asimov
```

Selecting BDT working points and particular signals is easily done via wildcard patterns to select the right datacards, e.g.:

```bash
python3 cli_boosted.py likelihood_scan dc_Dec07_minmt300/dc_mz*rinv0.3*bdt0p{0,3,5}*.txt --asimov --minmu -.5 --maxmu .5 -n 100
```

Note also the options `--minmu` and `--maxmu` which handle the range of the signal parameter to scan, and the option `-n` which controls the number of points in the range.


### Bias study

Generate toys:

```bash
python3 cli_boosted.py gentoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt -t 300 --expectSignal 0 -s 1001 --rMin -2 --rMax 3
for DC in dc_20240920_cutbased/*.txt; do python3 cli_boosted.py gentoys $DC -t 300 --expectSignal 0 -s 1001 --rMin -2 --rMax 3; done
```

Fit the toys:

```bash
python3 cli_boosted.py fittoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.txt --toysFile toys_20240924/higgsCombineObserveddc_SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.GenerateOnly.mH120.1001.root --expectSignal 0 --rMin -2 --rMax 3
 
for MZ in 200 250 300 350 400 450 500 550; do python3 cli_boosted.py fittoys dc_20240920_cutbased/dc_SVJ_s-channel_mMed-$MZ*.txt --toysFile toys_20240924/higgsCombineObserveddc_SVJ_s-channel_mMed-$MZ*.root --expectSignal 0 --rMin -2 --rMax 3; done
```

### Nuisance impacts

The nuisance impacts show how much each systematic uncertainty impacts the fit, and whether any systematics are pulled or constrained.
The input is a datacard txt file.
Signal injection is required to see the effects of signal systematics (the uncertainty variations are checked at the likelihood minimum, which is `r=0` without signal injection).

```bash
python3 cli_boosted.py impacts dc.txt --nfits 16 --asimov --normRange 0.1 2.0 --rMin -10 --rMax 10 --robustFit 1 --expectSignal 0.2
```

### Plotting

Fit to background distribution:
```bash
python3 quick_plot.py bkgfit ua2 --bkg cutbased/merged_20240920/bkg_sel-cutbased.json --sig cutbased/smooth_20240920/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json --outfile fit_bkg.png
```

ΔNNL as a function of mu:
```bash
python3 quick_plot.py muscan scans_20241029/*bdt0p3*Scan*.root
```

Single parameter as a function of mu:
```bash
python3 quick_plot.py trackedparam bsvj_bkgfitua2_npars3_p1 scans_cutbased/*Scan*.root
```

All fit parameters as a function of mu:
```bash
python3 quick_plot.py trackedparams scans_cutbased/*Scan*.root
```

MT histogram, with bkg-only fit and and sig+bkg fit:

```bash
python3 quick_plot.py mtdist scans_20241029/higgsCombineObserved_dc_mz450_rinv0.3_bdt0p300Bestfit.MultiDimFit.mH120.root
```

Note you should use the `Bestfit`-tagged file, not `Scan`.
Apparently, the single snapshot stored in the `Scan` files is _not_ the best fit.


CLS:

```bash
python3 quick_plot.py cls scans_20241029/higgsCombineObserved_dc_mz450_rinv0.3_bdt0p300.MultiDimFit.mH120.root scans_Dec07/higgsCombineAsimov_dc_mz450_rinv0.3_bdt0p300.MultiDimFit.mH120.root
```


Brazil band (relies on good interpolation; always check the CLs plots to double check!):

```bash
python3 quick_plot.py brazil scans_20241029/higgsCombine*bdt0p3*.root #Expects '*Observed*' and '*Asimov*' files
```

### High-level scripts

These scripts are made to easily execute required statistical tests for the boosted SVJ search.

#### Bias and Self Tests

The bias and self tests can be run via a shell script which makes datacards, generates toys, and fits toys automatically. 
The shell script is filled with default settings to run a self test from 200 to 550 with rinv=0.3 and mdark=10 with no signal injected for the bdt based search, 
but these options can easily configured in the command line without needing to alter the script.

```bash
# Default settings
./run_bias_or_self_study.sh

# Run the bias test for the cut based search with a r_inj at expected limit strength for mZ 200
./run_bias_or_self_study.sh --sel cutbased --test_type bias --rinj exp --mMed_values "200"
```

Then to plot the results

```bash
python3 plot_bias_or_self_study.py --base_dir ./self_test --sel bdt=0.67 --test self
```

This plotting script can be configured for the bias test, self test, bdt search, cutbased search, and various injected signal levels.

Example from the BDT based search:
```bash
# Self test
# no signal injected
./run_bias_or_self_test.sh --rinj 0
# signal injected at expected limit
./run_bias_or_self_test.sh --rinj exp
# plot
python3 plot_bias_or_self_study.py --base-dir ./self_test --sel bdt=0.67 --test self

# Bias test
# no signal injected
./run_bias_or_self_study.sh --test_type bias --rinj 0
# signal injected at expected limit
./run_bias_or_self_study.sh --test_type bias --rinj exp
# plot
python3 plot_bias_or_self_study.py --base-dir ./bias_test --sel bdt=0.67 --test bias
```

#### Limits: Expected + Observed = Asimov toy with signal injected at 350

This shell script works in a similar way to the previous in that it is loaded with a number of default settings that can be re-configured on the 
command line. It defaults to running from 200 to 550 with rinv=0.3 and mdark=10 for the bdt based search. Signal injection mass and a strength
and mass can be configured.

```bash
# Running with signal injected at mZ' 200 with a strength of 0.271
./run_limit_asimov_sig_inj.sh --mInj 200 --rinj 0.271
# Or to run all signal masses (at expected limit strength)
./all_sig_inj_asimov_limits.sh
```

If you redo the tests multiple times, make sure to provide clean directories

