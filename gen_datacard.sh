#!/bin/bash

# Default values
sel="bdt=0.75"
skip_toy=false # option to reuse previously generated toy
skip_dc=false # option to do just diagnostics
skip_diagnostics=false # option to do just datacards
toy_date=$(date +%Y%m%d) # in case different from dc date in skip_toy case
dc_date=$(date +%Y%m%d) # manual date for diagnostic plot commands (in skip_dc case)
suff=
hists_dir="hists"
hists_date= # date of the histograms used to make the data cards
hists_date_anti= # date of the anti-tag CR histograms
toys=
toy_seed=1001
mMed=350
mDark="10"
rinv="0p3"
anti=

# arguments specific to datacard creation and diagnostics
tf_basis=Bernstein
tf_mc=
npar_mc=
npar_data=5 # default max for data ftest
ftest_data=true

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sel) sel="$2"; shift ;;
        --skip_dc) skip_dc=true ;;
        --skip_diagnostics) skip_diagnostics=true ;;
        --toy_date) dc_date="$2"; shift ;;
        --dc_date) dc_date="$2"; shift ;;
        --suff) suff="$2"; shift ;;
        --hists_dir) hists_dir="$2"; shift ;;
        --hists_date) hists_date="$2"; shift ;;
        --hists_date_anti) hists_date_anti="$2"; shift ;;
        --mMed) mMed="$2"; shift ;;
        --mDark) mDark="$2"; shift ;;
        --rinv) rinv="$2"; shift ;;
        --antiloose) anti=antiloose;;
        --tf_basis) tf_basis="$2"; shift ;;
        --tf_mc) tf_mc=true ;;
        --npar_mc) npar_mc="$2"; shift ;;
        --npar_data) npar_data="$2"; ftest_data=; shift ;; # select npar, disable f-test
        --npar_data_max) npar_data="$2"; shift ;; # change max value, still run f-test
        --toy_seed) toy_seed="$2"; shift ;;
        --toys) toys="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ "$skip_dc" == true ] && [ "$skip_diagnostics" == true ]; then
    echo "Error: no actions selected"
    exit 1
fi

# dependent values
: ${anti:=anti}
if [ "$sel" == "cutbased" ]; then
    # assign defaults if empty
    : ${hists_date:=20241101}
    if [ "$anti" == "antiloose" ]; then
        : ${hists_date_anti:=20250218}
    else
        : ${hists_date_anti:=20250206}
    fi
elif [ "$sel" == "cutbased_ddt" ]; then
    # assign defaults if empty
    : ${hists_date:=20250710}
    : ${hists_date_anti:="$hists_date"}
elif [ "$sel" == "cutbased_ddt=0.11" ] || [ "$sel" == "cutbased_ddt=0.12" ]; then
    # assign defaults if empty
    : ${hists_date:=20250715}
    : ${hists_date_anti:="$hists_date"}
elif [ "$sel" == "bdt=0.55" ] || [ "$sel" == "bdt=0.75" ]; then
    # assign defaults if empty
    : ${hists_date:=20250715}
    : ${hists_date_anti:="$hists_date"}
else
    echo "Unknown sel $sel"
    exit 1
fi
antisel=${anti}${sel}

# assemble arguments
toy_args="-s ${toy_seed} --expectSignal 0"
# todo: allow both npar_data and ftest
dc_args="--basis ${tf_basis} --basis-mc ${tf_basis} ${npar_data:+--npar ${npar_data}} ${ftest_data:+--ftest} ${tf_mc:+--tf-from-mc} ${toy_args} ${toys:+-t $toys} ${npar_mc:+--winner tf_mc $npar_mc} ${suff:+--suff $suff}"

# account for change in histogram production/naming
hists_name=
hists_name_anti=
if [ $hists_date -gt 20250401 ]; then
    hists_name=_mt
fi
if [ $hists_date_anti -gt 20250401 ]; then
    hists_name_anti=_mt
fi

# directories and names
signame_base=SVJ_s-channel_mMed-${mMed}_mDark-${mDark}_rinv-${rinv}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8

mergeddir=${hists_dir}/merged_${hists_date}
smoothdir=${hists_dir}/smooth_${hists_date}
signame=${signame_base}_sel-${sel}${hists_name}_smooth
bkgname=bkg_sel-${sel}${hists_name}
bkg=${mergeddir}/${bkgname}.json
bkgsmooth=${smoothdir}/${bkgname}_wide_smooth.json
sig=${smoothdir}/${signame}.json

antimergeddir=${hists_dir}/merged_${hists_date_anti}
antismoothdir=${hists_dir}/smooth_${hists_date_anti}
antisigname=${signame_base}_sel-${antisel}${hists_name_anti}_smooth
antibkgname=bkg_sel-${antisel}${hists_name_anti}
antibkg=${antimergeddir}/${antibkgname}.json
antibkgsmooth=${antismoothdir}/${antibkgname}_wide_smooth.json
antisig=${antismoothdir}/${antisigname}.json

region_args_base="--regions ${sel} ${antisel} --sig ${sig} ${antisig}"
region_args1="${region_args_base} --bkg ${bkgsmooth} ${antibkgsmooth}"
dc_dir=dc_${dc_date}_${sel}
signame_dc=${signame}${suff:+_$suff}
dc_name=dc_${signame_dc}.txt

data_toy_suff=${suff:+$suff_}simple
signame_toy=${signame}_${data_toy_suff}
data_toy_args="${toy_args} -t 1"
data_toy_file=toys_${toy_date}/higgsCombineObserveddc_${signame_toy}.GenerateOnly.mH120.${toy_seed}.root
dc_toy_name=dc_${signame_toy}.txt
region_args2="${region_args_base} --bkg ${bkg} ${antibkg} --data ${data_toy_file} ${data_toy_file}"

bf_dir=bestfits_${dc_date}
bf_file=${bf_dir}/higgsCombineObservedBestfit_dc_${signame_dc}.MultiDimFit.mH120.root

set -x

# 0. toy generation
if [ "$skip_toy" == false ]; then
    python3 cli_boosted.py gen_datacards --norm-type crsimple --suff ${data_toy_suff} ${region_args1}

    python3 cli_boosted.py gentoys ${dc_dir}/${dc_toy_name} ${data_toy_args}
fi

# 1. datacard generation
if [ "$skip_dc" == false ]; then
    python3 cli_boosted.py gen_datacards --norm-type rhalpha ${region_args2} ${dc_args}
fi

# 2. diagnostics
if [ "$skip_diagnostics" == false ]; then
    # run bestfit
    python3 cli_boosted.py bestfit ${dc_dir}/${dc_name} --range -1 1

    # hessian analysis
    python3 hessian.py -w ${bf_file}:w -f ${bf_file}:fit_mdf -s 0.1

    # make plots

    # SR vs. CR distributions
    python3 quick_plot.py bkgsrcr ${region_args2} -o ${dc_dir}/srcr_${sel}.png

    # TF fit(s)
    python3 quick_plot.py bkgtf ${region_args2} -o ${dc_dir}/tf_${signame_dc}.png --fit-data ${bf_file}:fit_mdf:w ${tf_mc:+--fit-mc ${dc_dir}/bkgfit_${signame_dc}.root:bkgfit}

    # postfit
    for channel in bsvj bsvjCR1; do
        python3 quick_plot.py mtdist ${bf_file} --channel ${channel} --outfile ${bf_dir}/bestfit_${channel}_${signame_dc}.png
    done

    # todo: move all plots into one folder?
fi
