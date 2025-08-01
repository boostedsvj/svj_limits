import os,sys
import argparse
from time import strftime
from copy import deepcopy
from dataclasses import dataclass
from boosted_fits import run_generic_command as run_cmd, blank_logger as logger

# default values
# format: X, antiX, antilooseX
# "" = use X value, None = no default histograms exist
hists_dates = {
    "cutbased": ("20241101", "20250218", "20250206"),
    "cutbased_ddt": ("20250710", "", ""),
    "cutbased_ddt=0.11": ("20250715", "", ""),
    "cutbased_ddt=0.12": ("20250711", "", ""),
    "bdt=0.55": ("20250715", "", ""),
    "bdt=0.75": ("20250711", "", ""),
}
def safe_len(val): return -1 if val is None else len(val)
hists_dates = {key : {"": val[0], "anti": val[0] if safe_len(val[1])==0 else val[1], "antiloose": val[0] if safe_len(val[2])==0 else val[2]} for key,val in hists_dates.items()}

allowed_basis = ["Bernstein", "Chebyshev"]
today = strftime('%Y%m%d')
change_date = 20250401
qty = "mt"

class ArgumentDefaultsRawHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
    """HelpFormatter that adds default values AND doesn't do line-wrapping"""
pass

def join_none(sep, arr):
    return sep.join(filter(None,arr))

@dataclass
class Signal:
    mMed: str
    mDark: str
    rinv: str

    def mMed_val(self) -> int:
        return int(self.mMed)
    def mDark_val(self) -> int:
        return int(self.mDark)
    def rinv_val(self) -> float:
        return float(self.rinv.replace('p','.'))

# run a command over potentially multiple signal arguments; variants:
# 1. "single": runs once even when multiple signals specified
# 2. "loop": runs in serial over signals in a loop
# 3. "mp": runs in parallel over signals using multiprocessing version of function
class Command:
    def __init__(self, pre, function, flags, function_mp=None, single=False):
        self.pre = pre
        self.function = function
        self.function_mp = function_mp
        self.flags = flags
        self.single = single

    def run(self, args, signals):
        # "single"
        if self.single:
            args_signal = fill_signal_args(args, signals[0])
            self.run_single(args_signal, self.function)
        # "loop"
        elif args.npool==0 or self.function_mp is None:
            for signal in signals:
                args_signal = fill_signal_args(args, signal)
                self.run_single(args_signal, self.function)
        # "mp"
        else:
            args_signal = None
            for signal in signals:
                args_tmp = fill_signal_args(args, signal)
                if args_signal is None: args_signal = args_tmp
                else:
                    for key in vars(args_signal).keys(): setattr(args_signal, key, join_none(" ",[getattr(args_signal, key), getattr(args_tmp, key)]))
            self.run_single(args_signal, self.function_mp)

    def run_single(self, args, fn):
        cmd_template = join_none(" ",[self.pre,fn,self.flags])
        run_cmd(cmd_template.format(**vars(args)), cmd_logger=logger)

# runs a list of commands
class StepRunner:
    def __init__(self, name, commands):
        self.name = name
        self.commands = commands

    def run(self, args, signals):
        for command in self.commands:
            command.run(args, signals)

steps = {}
steps['0'] = StepRunner('pseudodata', [
    Command("python3 cli_boosted.py", "gen_datacards", "--norm-type crsimple --suff {data_toy_suff} {region_args1} --sig {regions_sig}", single=True),
    Command("python3 cli_boosted.py", "gentoys", "{dc_dir}/{dc_toy_name} {data_toy_args}", single=True),
    Command("mv", "", "{data_toy_file_old} {data_toy_file}", single=True)
])
steps['1'] = StepRunner('datacard generation', [
    Command("python3 cli_boosted.py", "gen_datacards", "--norm-type rhalpha {region_args2} --sig {regions_sig} {dc_args}"),
])
steps['2'] = StepRunner('diagnostics', [
    # run bestfit
    Command("python3 cli_boosted.py", "bestfit", "--range -1 1 {dc_dir}/{dc_name}"),
    # hessian analysis
    Command("python3", "hessian.py", "-w {bf_file}:w -f {bf_file}:fit_mdf -s 0.1"),
    # make plots
    # SR vs. CR distributions
    Command("python3 quick_plot.py", "bkgsrcr", "{region_args2} --sig {regions_sig} -o {dc_dir}/srcr_{sel}.png"),
    # TF fit(s)
    Command("python3 quick_plot.py", "bkgtf", "{region_args2} --sig {regions_sig} -o {dc_dir}/tf_{signame_dc}.png --basis {tf_basis} --basis-mc {tf_basis} --fit-data {bf_file}:fit_mdf:w {fit_mc_arg}"),
    ] + [
    # postfit
    Command("python3 quick_plot.py", "mtdist", "{{bf_file}} --sel {0} --channel {1} --outfile {{bf_dir}}/bestfit_{1}_{{signame_dc}}.png".format(sel, channel))
        for sel, channel in [("{sel}", "bsvj"), ("{antisel}", "bsvjCR1")]
    ]
)
# todo: move all plots into one folder?

# special groups of steps
predefs = {
    'gen_datacard': ['0','1','2'],
    'gen_datacard_alt': ['0','1b','2b'],
    'closure': ['1'],
    'bias': ['1','1b'],
}

def fill_signal_args(args, signal):
    signal_args = deepcopy(args)

    # directories and names
    signame_base = f"SVJ_s-channel_mMed-{signal.mMed}_mDark-{signal.mDark}_rinv-{signal.rinv}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8"

    signame = f"{signame_base}_sel-{args.sel}{args.hists_name}_smooth"
    signal_args.sig = f"{args.smoothdir}/{signame}.json"

    antisigname = f"{signame_base}_sel-{args.antisel}{args.hists_name_anti}_smooth"
    signal_args.antisig = f"{args.antismoothdir}/{antisigname}.json"

    signal_args.regions_sig = f"{signal_args.sig} {signal_args.antisig}"
    signal_args.signame_dc = join_none("_",[signame, args.suff])
    signal_args.dc_name = f"dc_{signal_args.signame_dc}.txt"

    signal_args.signame_toy = f"{signame}_{args.data_toy_suff}"
    signal_args.data_toy_file_old = f"toys_{args.toy_date}/higgsCombineObserveddc_{signal_args.signame_toy}.GenerateOnly.mH120.{args.toy_seed}.root"
    signal_args.dc_toy_name = f"dc_{signal_args.signame_toy}.txt"

    signal_args.bf_file = f"{args.bf_dir}/higgsCombineObservedBestfit_dc_{signal_args.signame_dc}.MultiDimFit.mH120.root"
    signal_args.fit_mc_arg = f"--fit-mc {args.dc_dir}/bkgfit_{signal_args.signame_dc}.root:bkgfit" if args.tf_mc else ""

    return signal_args

# todo:
# option to swap out all pseudodata-related toy args for real data args
def derive_args(args_orig, alt=False):
    args = deepcopy(args_orig)

    # swap to alt versions for bias study
    if alt:
        args.tf_basis = allowed_basis[1] if args.tf_basis==allowed_basis[0] else allowed_basis[0]
        args.suff = join_none("_", [args.suff, 'alt'])

    # derived values
    args.anti = "antiloose" if args.antiloose else "anti"
    args.antisel = args.anti+args.sel
    args.hists_date = args.hists_date or hists_dates[args.sel][""]
    args.hists_date_anti = args.hists_date_anti or hists_dates[args.sel][args.anti]
    args.ftest = not args.npar_data
    args.npar_data = args.npar_data or args.npar_data_max

    # assemble arguments
    args.toy_args = f"-s {args.toy_seed} --expectSignal 0"
    args.dc_args = join_none(" ",[
        f"--basis-mc {args.tf_basis}",
        f"--tf-from-mc" if args.tf_mc else "",
        f"--npar-mc-max {args.npar_mc_max}",
        f"--winner tf_mc {args.npar_mc-1}" if args.npar_mc else "",
        f"--basis {args.tf_basis}",
        f"--npar {args.npar_data-1}",
        f"--ftest" if args.ftest else "",
        args.toy_args,
        f"-t {args.toys}" if args.toys else "",
        f"--suff {args.suff}" if args.suff else "",
    ])

    # account for change in histogram production/naming
    args.hists_name = f"_{qty}" if int(args.hists_date) > change_date else ""
    args.hists_name_anti = f"_{qty}" if int(args.hists_date_anti) > change_date else ""

    # directories and names
    args.mergeddir = f"{args.hists_dir}/merged_{args.hists_date}"
    args.smoothdir = f"{args.hists_dir}/smooth_{args.hists_date}"
    args.bkgname = f"bkg_sel-{args.sel}{args.hists_name}"
    args.bkg = f"{args.mergeddir}/{args.bkgname}.json"
    args.bkgsmooth = f"{args.smoothdir}/{args.bkgname}_wide_smooth.json"

    args.antimergeddir = f"{args.hists_dir}/merged_{args.hists_date_anti}"
    args.antismoothdir = f"{args.hists_dir}/smooth_{args.hists_date_anti}"
    args.antibkgname = f"bkg_sel-{args.antisel}{args.hists_name_anti}"
    args.antibkg = f"{args.antimergeddir}/{args.antibkgname}.json"
    args.antibkgsmooth = f"{args.antismoothdir}/{args.antibkgname}_wide_smooth.json"

    args.region_args_base = f"--regions {args.sel} {args.antisel}"
    args.region_args1 = f"{args.region_args_base} --bkg {args.bkgsmooth} {args.antibkgsmooth}"
    args.dc_dir = f"dc_{args.dc_date}_{args.sel}"

    args.data_toy_suff = join_none("_",[args.suff,"simple"])
    args.data_toy_args = f"{args.toy_args} -t 1"
    args.data_toy_file = f"toys_{args.toy_date}/higgsCombineObserveddc_bkg.GenerateOnly.mH120.{args.toy_seed}.root"
    args.region_args2 = f"{args.region_args_base} --bkg {args.bkg} {args.antibkg} --data {args.data_toy_file} {args.data_toy_file}"

    args.bf_dir = f"bestfits_{args.dc_date}"

    return args

if __name__=="__main__":
    desc = [
        "\n".join(["Steps:"]+[f"{key}. {val.name}" for key,val in sorted(steps.items())]),
        "(add 'b' to any step to use alt TF basis)",
        "",
        "\n".join(["Predefs:"]+[f"{key} = "+" ".join(val) for key,val in predefs.items()]),
        "",
        "\n".join(["Current defaults:", "* hists_dates: "]+[
            f"\t{key}: "+", ".join([
                f"{subkey} = {subval}" if subkey else f"{subval}" for subkey,subval in hists_dates[key].items()
            ]) for key in hists_dates.keys()
        ]),
    ]
    parser = argparse.ArgumentParser(
        formatter_class=ArgumentDefaultsRawHelpFormatter,
        description="\n".join(desc),
    )
    group_co = parser.add_argument_group("common")
    group_cx = group_co.add_mutually_exclusive_group()
    group_cx.add_argument("--steps", type=str, nargs='+', help="step(s) to run")
    group_cx.add_argument("--predef", type=str, default="", choices=list(predefs.keys()), help="predefined sequence to run")
    group_co.add_argument("--sel", type=str, required=True, help="selection name")
    group_co.add_argument("--antiloose", default=False, action="store_true", help="use antiloose region")
    group_co.add_argument("--suff", type=str, default="", help="suffix for dc")
    group_co.add_argument("--npool", type=int, default=0, help="number of cores for multiprocessing pool")
    group_co.add_argument("--hists-dir", type=str, default="hists", help="histogram directory")
    group_co.add_argument("--hists-date", type=str, default=None, help="date for signal region histograms")
    group_co.add_argument("--hists-date-anti", type=str, default=None, help="date for anti-signal region histograms")
    group_co.add_argument("--toy-date", type=str, default=today, help="date for toy folder (if skipping step 0)")
    group_co.add_argument("--dc-date", type=str, default=today, help="date for dc folder (if skipping step 1)")
    group_si = parser.add_argument_group("signal")
    group_sx = group_si.add_mutually_exclusive_group()
    group_sx.add_argument("--signal", dest="signals", metavar=("mMed","mDark","rinv"), type=str, default=["350","10","0p3"], nargs=3, help="signal parameters")
    group_sx.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    group_dc = parser.add_argument_group("datacard")
    group_dc.add_argument("--tf-basis", type=str, default=allowed_basis[0], choices=allowed_basis, help="transfer factor polynomial basis")
    group_dc.add_argument("--tf-mc", default=False, action="store_true", help="use TF from MC")
    group_dc.add_argument("--npar-mc", type=int, default=None, help="number of parameters for MC TF, if not using F-test")
    group_dc.add_argument("--npar-mc-max", type=int, default=6, help="max number of parameters for MC TF F-test")
    group_dc.add_argument("--npar-data", type=int, default=None, help="number of parameters for data TF, if not using F-test")
    group_dc.add_argument("--npar-data-max", type=int, default=5, help="max number of parameters for data TF F-test")
    group_dc.add_argument("--toys", type=int, default=0, help="number of toys for data TF F-test")
    group_dc.add_argument("--toy-seed", type=int, default=1001, help="random seed for toy generation")
    args = parser.parse_args()
    if args.predef: args.steps = predefs[args.predef]

    # signals
    signals = []
    if isinstance(args.signals,list):
        signals.append(Signal(*args.signals))
    else:
        with open(args.signals,'r') as sfile:
            for line in sfile:
                line = line.rstrip()
                if len(line)==0: continue
                signals.append(Signal(*line.split()))

    # execute requested steps in order
    for step in args.steps:
        alt = step.endswith('b')
        step = step.replace('b','')
        args_step = derive_args(args, alt=alt)
        steps[step].run(args_step, signals)
