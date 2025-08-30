import os,sys,re
import argparse
from time import strftime
from copy import deepcopy
from dataclasses import dataclass
from collections import defaultdict
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

def nat_sort(val):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', val)]

@dataclass
class Signal:
    mMed: str
    mDark: str
    rinv: str
    exp: str

    def mMed_val(self) -> int:
        return int(self.mMed)
    def mDark_val(self) -> int:
        return int(self.mDark)
    def rinv_val(self) -> float:
        return float(self.rinv.replace('p','.'))
    def exp_val(self) -> float:
        return float(self.exp)

def get_rinj(rinj, signal):
    if rinj<0:
        exp_val = signal.exp_val()
        return exp_val*abs(rinj)
    else:
        return rinj

def handle_signals_explim(args):
    if args.explim:
        args.explim = args.explim.format(**vars(args))
        # explim is the input file name (used in asimov_inj, self/bias, etc.)
        # explim_name is the output file name (created during likelihood scan for use in above steps)
        args.explim_name = args.explim[:]
        if not os.path.isfile(args.explim):
            logger.warning(f"explim file {args.explim} not found")
            args.explim = ""

    # signals
    explim_default = '0.2'
    explims = defaultdict(lambda: explim_default)
    if args.explim:
        with open(args.explim,'r') as efile:
            for line in efile:
                line = line.rstrip()
                props = line.split()
                explims[tuple(props[:-1])] = props[-1]
    else:
        logger.warning(f"explim file not provided; commands with --rinj -x will use default {explim_default}")
    def make_signal(props):
        # get expected limit strength
        props.append(explims[tuple(props)])
        return Signal(*props)
    signals = []
    if isinstance(args.signals,list):
        signals.append(make_signal(args.signals))
    else:
        with open(args.signals,'r') as sfile:
            for line in sfile:
                line = line.rstrip()
                if len(line)==0: continue
                props = line.split()
                signals.append(make_signal(props))

    return args, signals, explims

def get_signame(signal):
    return f"SVJ_s-channel_mMed-{signal.mMed}_mDark-{signal.mDark}_rinv-{signal.rinv}_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8"

# run a command over potentially multiple signal arguments; variants for broadcasting across multiple signals:
# 1. "single": runs once even when multiple signals specified (either independent of signal, or aggregates across all signals)
# 2. "loop": runs in serial over signals in a loop (default)
# 3. "mp": runs in parallel over signals using multiprocessing version of function
class Command:
    def __init__(self, pre, function, flags, cast='loop'):
        self.pre = pre
        self.function = function
        self.flags = flags
        self.cast = cast

    def get_flags(self, args):
        return self.flags.format(**vars(args))

    def run(self, args, signals, dryrun):
        # "single"
        if self.cast=='single':
            args_signal = fill_signal_args(args, signals[0])
            flags_signal = self.get_flags(args_signal)
            self.run_single(self.function, flags_signal, dryrun)
        # "loop"
        elif self.cast=='loop' or args.npool==0:
            for signal in signals:
                args_signal = fill_signal_args(args, signal)
                flags_signal = self.get_flags(args_signal)
                self.run_single(self.function, flags_signal, dryrun)
        # "mp"
        elif self.cast=='mp':
            # write out args to txt file
            argdir = 'args'
            os.makedirs(argdir, exist_ok=True)
            arg_filename = f"{argdir}/{self.function}_{strftime('%Y%m%d%H%M%S')}.txt"
            with open(arg_filename, 'w') as arg_file:
                for signal in signals:
                    args_signal = fill_signal_args(args, signal)
                    arg_file.write(args_signal.signame_base+' '+self.flags.format(**vars(args_signal))+'\n')

            # call general run_mp function
            self.run_single('run_mp', f'--npool {args.npool} {self.function} {arg_filename}', dryrun)
        else:
            raise RuntimeError(f"Unknown broadcast value {self.cast}")

    def run_single(self, fn, flags, dryrun):
        cmd = join_none(" ",[self.pre,fn,flags])
        if dryrun:
            logger.info(cmd)
        else:
            run_cmd(cmd, cmd_logger=logger)

# runs a list of commands
class StepRunner:
    def __init__(self, name, commands):
        self.name = name
        self.commands = commands

    def run(self, args, signals, dryrun):
        for command in self.commands:
            command.run(args, signals, dryrun)

steps = {}
steps['0'] = StepRunner('pseudodata', [
    Command("python3 cli_boosted.py", "gen_datacards", "--norm-type crsimple --suff {data_toy_suff} {region_args1} --sig {regions_sig}", cast='single'),
    Command("python3 cli_boosted.py", "gentoys", "{dc_dir}/{dc_toy_name} {data_toy_args}", cast='single'),
    Command("mv", "", "{data_toy_file_old} {data_toy_file}", cast='single'),
])
steps['1'] = StepRunner('datacard generation', [
    Command("python3 cli_boosted.py", "gen_datacards", "--norm-type rhalpha {region_args2} --sig {regions_sig} {dc_args}", cast='mp'),
])
steps['2'] = StepRunner('diagnostics', [
    # run bestfit
    Command("python3 cli_boosted.py", "bestfit", "{dc_dir}/{dc_name} --range -1 1", cast='mp'),
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
    # todo: plot toy-based f-test distribution
    # todo: plot f-test results vs. signal parameter(s)
    # todo: move all plots into one folder?
)
steps['3'] = StepRunner('bias toys', [
    Command("python3 cli_boosted.py", "gentoys", "{dc_dir}/{dc_name} {bias_toy_args} {rinj_arg}", cast='mp'),
    Command("mv", "", "{bias_toy_file_old} {bias_toy_file}"),
])
steps['3a'] = StepRunner('Asimov toy', [
    Command("python3 cli_boosted.py", "gentoys", "{dc_dir}/{scan_dc_name} {scan_toy_args} {rinj_arg}", cast='single'),
    Command("mv", "", "{scan_toy_file_old} {scan_toy_file}", cast='single'),
])
steps['4'] = StepRunner('likelihood scan', [
    Command("python3 cli_boosted.py", "likelihood_scan", "{dc_dir}/{dc_name} {scan_args}", cast='mp'),
    # dump expected limit signal strengths into a file
    Command("python3 quick_plot.py", "explim", "{all_scan_files} -o {explim_name}", cast='single'),
])
# todo: add "observed" likelihood command (using pseudodata/toy)
steps['5'] = StepRunner('likelihood plots', [
    Command("python3 quick_plot.py", "muscan", "{scan_files} -o {scan_dir}/muscan_{signame_dc}.png"),
    Command("python3 quick_plot.py", "cls", "{scan_files} -o {scan_dir}/cls_{signame_dc}.png"),
    Command("python3 quick_plot.py", "trackedparams", "{scan_files} -o {scan_dir}/{{}}_{signame_dc}.png"),
])
steps['6'] = StepRunner('Asimov injection test', [
    Command("python3 cli_boosted.py", "likelihood_scan", "{dc_dir}/{dc_name} {scan_inj_args}", cast='mp'),
])
steps['7'] = StepRunner('Asimov injection plots', [
    Command("python3 quick_plot.py", "brazil", "{all_scan_files} -o {scan_dir}/asimov__inj_{scan_inj_name_short}__sel-{sel}.png", cast='single'),
])
steps['8'] = StepRunner('bias fits', [
    Command("python3 cli_boosted.py", "fittoys", "{dc_dir}/{dc_name} {bias_fit_args} {bias_sig_args}", cast='mp'),
    Command("mkdir", "", "-p {bias_results_dir}", cast='single'),
    Command("mv", "", "{bias_fit_file} {bias_results_dir}/"),
])
steps['9'] = StepRunner('bias plots', [
    Command("python3", "plot_bias_or_self_study.py", "--base-dir {bias_results_basedir} --sel {sel} --test {bias_test_type} -s {btoy_seed} --explim {explim} --signals {signals}", cast='single'),
])
# todo: nuisance impacts

# special groups of steps
predefs = {
    'gen_datacard': ['0','1','2'],
    'gen_datacard_alt': ['0','1b','2b'],
    'likelihood': ['0','1','4','5'],
    'asimov_inj': ['0','1','3a','6','7'],
    'self': ['0','1','3','8','9'],
    'bias': ['0','1','1b','3b','8','9b'],
}

def fill_signal_args(args, signal):
    signal_args = deepcopy(args)

    # directories and names
    signal_args.signame_base = get_signame(signal)

    signal_args.signame = f"{signal_args.signame_base}_sel-{args.sel}{args.hists_name}_smooth"
    signal_args.sig = f"{args.smoothdir}/{signal_args.signame}.json"

    antisigname = f"{signal_args.signame_base}_sel-{args.antisel}{args.hists_name_anti}_smooth"
    signal_args.antisig = f"{args.antismoothdir}/{antisigname}.json"

    signal_args.regions_sig = f"{signal_args.sig} {signal_args.antisig}"
    signal_args.signame_dc = join_none("_",[signal_args.signame, args.suff])
    signal_args.dc_name = f"dc_{signal_args.signame_dc}.txt"

    signal_args.signame_dtoy = f"{signal_args.signame}_{args.data_toy_suff}"
    signal_args.data_toy_file_old = args.data_toy_file.replace("_bkg.",f"_{signal_args.signame_dtoy}.")
    signal_args.dc_toy_name = f"dc_{signal_args.signame_dtoy}.txt"

    signal_args.bf_file = f"{args.bf_dir}/higgsCombineObservedBestfit_dc_{signal_args.signame_dc}.MultiDimFit.mH120.root"
    signal_args.fit_mc_arg = f"--fit-mc {args.dc_dir}/bkgfit_{signal_args.signame_dc}.root:bkgfit" if args.tf_mc else ""

    signal_args.rinj_arg = f"--expectSignal {get_rinj(args.rinj,signal)}"
    signal_args.signame_btoy = join_none("_", [signal_args.signame, args.suff])
    signal_args.bias_toy_file_old = f"toys_{args.btoy_date}/higgsCombineObserveddc_{signal_args.signame_btoy}.GenerateOnly.mH120.{args.btoy_seed}.root"
    rinjname = "Exp" if args.rinj<0 else args.rinj
    signal_args.bias_toy_file = signal_args.bias_toy_file_old.replace(".GenerateOnly", f"_rinj{rinjname}.GenerateOnly")

    signal_args.bias_sig_args = f"--toysFile {signal_args.bias_toy_file} --expectSignal 0"
    signal_args.bias_fit_file = f"toyfits_{args.bfit_date}/higgsCombineObserveddc_{signal_args.signame_btoy}.FitDiagnostics.mH120.{args.btoy_seed}.root"

    signal_args._scan_toy_file_old = f"toys_{args.stoy_date}/higgsCombineAsimovdc_{signal_args.signame}.GenerateOnly.mH120.{args.stoy_seed}.root"
    signal_args._scan_toy_file = signal_args._scan_toy_file_old.replace(".GenerateOnly", f"_rinj{rinjname}.GenerateOnly")
    signal_args.scan_files = f"{args.scan_dir}/higgsCombinedc_{signal_args.signame}ScanObserved.MultiDimFit.mH120.{args.stoy_seed}.root"
    signal_args.scan_files += f" {signal_args.scan_files.replace('Observed','Asimov')}"

    return signal_args

# todo:
# option to swap out all pseudodata-related toy args for real data args
def derive_args(args_orig, signals, alt=False):
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
    args.ftoy_args = f"-s {args.ftoy_seed} --expectSignal 0"
    args.dc_args = join_none(" ",[
        f"--basis-mc {args.tf_basis}",
        f"--tf-from-mc" if args.tf_mc else "",
        f"--npar-mc-max {args.npar_mc_max}",
        f"--winner tf_mc {args.npar_mc-1}" if args.npar_mc else "",
        f"--basis {args.tf_basis}",
        f"--npar {args.npar_data-1}",
        f"--ftest" if args.ftest else "",
        args.ftoy_args,
        f"-t {args.ftoys}" if args.ftoys else "",
        f"--suff {args.suff}" if args.suff else "",
    ])
    args.bias_toy_args = f"-t {args.btoys} -s {args.btoy_seed}"
    args.bias_fit_args = f"{args.bias_toy_args} --range {args.brange[0]} {args.brange[1]} --cminDefaultMinimizerStrategy 0"
    args.scan_toy_args = f"--asimov -s {args.stoy_seed}"
    args.scan_args_base = f"--range {args.srange[0]} {args.srange[1]} -s {args.stoy_seed}"
    args.scan_args = f"{args.scan_args_base} --asimov"
    args.scan_inj_args = f"{args.scan_args_base} -t -1"

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
    args.data_toy_args = f"-s {args.dtoy_seed} --expectSignal 0 -t 1"
    args.data_toy_file = f"toys_{args.dtoy_date}/higgsCombineObserveddc_bkg.GenerateOnly.mH120.{args.dtoy_seed}.root"
    args.region_args2 = f"{args.region_args_base} --bkg {args.bkg} {args.antibkg} --data {args.data_toy_file} {args.data_toy_file}"

    args.bf_dir = f"bestfits_{args.dc_date}"

    args.bias_fits_dir = f"toyfits_{args.bfit_date}"
    args.bias_test_type = "bias" if alt else "self"
    args.bias_results_rinj = "1" if args.rinj!=0 else "0"
    args.bias_results_basedir = f"{args.bias_test_type}_test"
    args.bias_results_dir = f"{args.bias_results_basedir}/rinj{args.bias_results_rinj}"

    args.scan_dir = f"scans_{args.scan_date}"
    # signal-injected asimov toy is shared
    signal_inj = Signal(*args.siginj)
    siginj_args = fill_signal_args(args, signal_inj)
    args.scan_dc_name = siginj_args.dc_name
    args.scan_toy_file_old = siginj_args._scan_toy_file_old
    args.scan_toy_file = siginj_args._scan_toy_file
    args.scan_inj_args = join_none(" ",[args.scan_inj_args, f"--toysFile {args.scan_toy_file}"])
    args.scan_inj_name_short = join_none('_',args.siginj)

    # agglomerate all signal scan result files
    args.all_scan_files = []
    for signal in signals:
        signal_args = fill_signal_args(args, signal)
        args.all_scan_files.append(signal_args.scan_files)
    args.all_scan_files = join_none(" ",args.all_scan_files)

    return args

if __name__=="__main__":
    default_signal = ["350","10","0p3"]

    desc = [
        "\n".join(["Steps:"]+[f"{key}. {val.name}" for key,val in sorted(steps.items(), key=lambda item: nat_sort(item[0]))]),
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
    group_co.add_argument("--dryrun", default=False, action="store_true", help="just print commands, don't run anything")
    group_co.add_argument("--sel", type=str, required=True, help="selection name")
    group_co.add_argument("--antiloose", default=False, action="store_true", help="use antiloose region")
    group_co.add_argument("--suff", type=str, default="", help="suffix for dc")
    group_co.add_argument("--npool", type=int, default=0, help="number of cores for multiprocessing pool")
    group_co.add_argument("--hists-dir", type=str, default="hists", help="histogram directory")
    group_co.add_argument("--hists-date", type=str, default=None, help="date for signal region histograms")
    group_co.add_argument("--hists-date-anti", type=str, default=None, help="date for anti-signal region histograms")
    group_co.add_argument("--dtoy-date", type=str, default=today, help="date for pseudodata toy folder (if skipping step 0)")
    group_co.add_argument("--dtoy-seed", type=int, default=1001, help="random seed for pseudodata toy generation")
    group_co.add_argument("--dc-date", type=str, default=today, help="date for dc folder (if skipping step 1)")
    group_si = parser.add_argument_group("signal")
    group_sx = group_si.add_mutually_exclusive_group()
    group_sx.add_argument("--signal", dest="signals", metavar=("mMed","mDark","rinv"), type=str, default=default_signal, nargs=3, help="signal parameters")
    group_sx.add_argument("--signals", dest="signals", type=str, default="", help="text file w/ list of signal parameters")
    group_si.add_argument("--explim", type=str, default="explim_{sel}.txt", help="generated file with expected limit values per signal from Asimov scan")
    group_dc = parser.add_argument_group("datacard")
    group_dc.add_argument("--tf-basis", type=str, default=allowed_basis[0], choices=allowed_basis, help="transfer factor polynomial basis")
    group_dc.add_argument("--tf-mc", default=False, action="store_true", help="use TF from MC")
    group_dc.add_argument("--npar-mc", type=int, default=None, help="number of parameters for MC TF, if not using F-test")
    group_dc.add_argument("--npar-mc-max", type=int, default=6, help="max number of parameters for MC TF F-test")
    group_dc.add_argument("--npar-data", type=int, default=None, help="number of parameters for data TF, if not using F-test")
    group_dc.add_argument("--npar-data-max", type=int, default=5, help="max number of parameters for data TF F-test")
    group_dc.add_argument("--ftoys", type=int, default=0, help="number of toys for data TF F-test")
    group_dc.add_argument("--ftoy-seed", type=int, default=1005, help="random seed for F-test toy generation")
    group_ty = parser.add_argument_group("toys")
    group_ty.add_argument("--rinj", type=float, default=0, help="toy signal injection strength; if -x, strength = expected limit * x")
    group_sc = parser.add_argument_group("scan")
    group_sc.add_argument("--stoy-date", type=str, default=today, help="date for toy folder for Asimov test (if skipping step)")
    group_sc.add_argument("--stoy-seed", type=int, default=1002, help="random seed for toy generation for Asimov test")
    group_sc.add_argument("--scan-date", type=str, default=today, help="date for scans (if skipping step)")
    group_sc.add_argument("--srange", type=float, default=[0,2], nargs=2, help="likelihood scan r range")
    group_sc.add_argument("--siginj", metavar=("mMed","mDark","rinv"), type=str, default=default_signal, nargs=3, help="signal to inject for Asimov test")
    group_bi = parser.add_argument_group("bias")
    group_bi.add_argument("--btoys", type=int, default=0, help="number of toys for bias test")
    group_bi.add_argument("--btoy-seed", type=int, default=1003, help="random seed for toy generation for bias test")
    group_bi.add_argument("--btoy-date", type=str, default=today, help="date for toy folder for bias test (if skipping step)")
    group_bi.add_argument("--bfit-date", type=str, default=today, help="date for bias fit folder (if skipping step)")
    group_bi.add_argument("--brange", type=float, default=[-5,5], nargs=2, help="bias study r range")
    # todo: expose pointsRandProf options
    args = parser.parse_args()
    if args.predef: args.steps = predefs[args.predef]

    args, signals, explims = handle_signals_explim(args)
    args.siginj.append(explims[tuple(args.siginj)])

    # execute requested steps in order
    for step in args.steps:
        alt = step.endswith('b')
        step = step.replace('b','')
        args_step = derive_args(args, signals, alt=alt)
        steps[step].run(args_step, signals, args.dryrun)
