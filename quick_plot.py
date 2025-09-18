"""
Some scripts to quickly plot basic outcomes from combine scans
"""
from __future__ import print_function
import ROOT # type:ignore
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
from time import strftime
import imp
import argparse

# Add the directory of this file to the path so the boosted tools can be imported
import sys, os, os.path as osp, pprint, re, traceback, copy, fnmatch, shutil
from contextlib import contextmanager
sys.path.append(osp.dirname(osp.abspath(__file__)))
import svj_ntuple_processing as svj
import boosted_fits as bsvj
import run_rhalpha as rhalph
logger = bsvj.setup_logger('quickplot')

import numpy as np
import matplotlib as mpl
mpl.use('Agg') # in order to run in background / no-graphics environments
import matplotlib.pyplot as plt # type:ignore
import matplotlib.transforms as transforms

def set_mpl_fontsize(small=22, medium=28, large=32, legend=None):
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=medium if legend is None else legend)    # legend fontsize
    plt.rc('figure', titlesize=large)  # fontsize of the figure title
set_mpl_fontsize()

cms_yellow = '#fffc54'
cms_green = '#74f94c'

def get_color_cycle():
    from itertools import cycle
    colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    return colors

scripter = bsvj.Scripter()

def cmd_exists(executable):
    if shutil.which(executable):
        return True
    else:
        return False

BATCH_MODE = False
def batch_mode(flag=True):
    global BATCH_MODE
    BATCH_MODE = bool(flag)

DEBUG = False
def debug(flag=True):
    global DEBUG
    DEBUG = bool(flag)


@contextmanager
def quick_ax(figsize=(12,12), outfile='test.png'):
    try:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
        yield ax
    finally:
        plt.savefig(outfile, bbox_inches='tight')
        if not(BATCH_MODE) and cmd_exists('imgcat'): os.system('imgcat ' + outfile)


def name_from_combine_rootfile(rootfile, strip_obs_asimov=False, sel=True):
    param_tex = {
        'mz': r'm_{\mathrm{Z}^{\prime}}',
        'mdark': r'm_{\mathrm{dark}}',
        'rinv': r'r_{\mathrm{inv}}',
    }
    meta = svj.metadata_from_path(rootfile)
    if meta['sample_type']=='sig':
        names = ['']
        if not strip_obs_asimov:
            names = [x for x in ['Observed','Asimov'] if x in rootfile]+names
        name = names[0]
        name += ' ($'+', '.join(val+'='+str(meta[key]) for key,val in param_tex.items())+'$)'
    else:
        name = osp.basename(rootfile).rsplit('.',3)[0].replace('higgsCombine','')
        if strip_obs_asimov:
            name = name.replace('Observed_','').replace('Asimov_','')
        name = name.replace('.MultiDimFit','')

    if sel:
        match = re.search(r'sel-([^_]*_ddt[0-9=.]*|[^_]*)', rootfile)
        if match: name = name+'\n'+match.group(1)
    return name

def namespace_to_attrdict(args):
    return bsvj.AttrDict(**vars(args))

def get_mz(path):
    return svj.metadata_from_path(path)['mz']

def get_rinv(path):
    return svj.metadata_from_path(path)['rinv']

def get_bdt_str(path):
    return re.search(r'bdt([p\d]+)', osp.basename(path)).group(1)

def get_bdt(path):
    return float(get_bdt_str(path).replace('p', '.'))

def organize_rootfiles(rootfiles, split_bdt_wps=False):
    """
    Takes an unordered set of rootfiles, and splits them up logically
    by obs/asimov, mz, and bdt
    """
    rootfiles_orig = rootfiles[:]
    rootfiles = []
    for rootfile in rootfiles_orig:
        if not osp.exists(rootfile):
            logger.warning(f"organize_rootfiles: {rootfile} not found, skipping")
            continue
        rootfiles.append(rootfile)

    if split_bdt_wps:
        bdt_wps = set(get_bdt(f) for f in rootfiles)
        logger.info('Found the following bdt wps: %s', bdt_wps)
        out = []
        for bdt_wp in sorted(bdt_wps):
            out.append(organize_rootfiles([r for r in rootfiles if get_bdt(r)==bdt_wp], False))
        for row in out:
            print(row)
        return out

    rootfiles.sort(key=get_mz)

    obs_rootfiles = [f for f in rootfiles if 'Observed' in osp.basename(f)]
    asimov_rootfiles = [f for f in rootfiles if 'Asimov' in osp.basename(f)]
    if len(obs_rootfiles)==0:
        obs_rootfiles = asimov_rootfiles[:]
    assert [get_mz(f) for f in obs_rootfiles] == [get_mz(f) for f in asimov_rootfiles]
    return obs_rootfiles, asimov_rootfiles


class Scan(object):
    """
    Poor man's pd.Dataframe to contain the data of a single likelihood scan
    """
    def __init__(self):
        self.df = {}

    def __getitem__(self, where):
        new = copy.copy(self)
        new.df = {key: val[where] for key, val in self.df.items()}
        return new

    @property
    def n(self):
        for arr in self.df.values():
            return len(arr)

def extract_scans(rootfiles, correct_minimum=False, keep_all=False):
    if isinstance(rootfiles, str): rootfiles = [rootfiles]
    scans = []

    for rootfile in rootfiles:
        with bsvj.open_root(rootfile) as tf:
            limit = tf.Get('limit')

            cscan = Scan()

            keys = {
                'r' : 'mu',
                'deltaNLL' : 'dnll',
                'quantileExpected' : 'quantile',
                'iToy' : 'itoy'
                }
            # Add the tracked params
            listofbranches = limit.GetListOfBranches()
            for i_branch in range(listofbranches.GetEntries()):
                branch = listofbranches[i_branch].GetName()
                if branch.startswith('trackedParam_'):
                    keys[branch] = branch.replace('trackedParam_','')
                elif branch.startswith('trackedError_'):
                    keys[branch] = branch.replace('trackedError_','error_')

            # Read in the values from the TTree
            for _ in limit:
                for tree_key, target_key in keys.items():
                    cscan.df.setdefault(target_key, [])
                    cscan.df[target_key].append( getattr(limit, tree_key) )
            # Turn into numpy arrays
            for key, val in cscan.df.items(): cscan.df[key] = np.array(val)

            # Split per toy
            for i_toy in set(cscan.df['itoy']):
                scan = cscan[cscan.df['itoy'] == i_toy]

                # Take out the bestfit
                is_bestfit = scan.df['quantile'] == -1
                assert is_bestfit.sum() == 1
                i_bestfit = is_bestfit.argmax()
                scan.bestfit = scan[i_bestfit]
                scan = scan[~is_bestfit]

                if correct_minimum:
                    i_minimum = np.argmin(scan.df['dnll'])
                    min_dnll = scan.df['dnll'][i_minimum]
                    scan.bestfit = scan[i_minimum]
                    logger.warning(
                        'Shifting curve by {:.4f} and setting bestfit to {:.4f}'
                        .format(min_dnll, scan.bestfit.df['mu']))
                    scan.df['dnll'] -= min_dnll

                scans.append(scan)

    # for debugging
    if keep_all:
        return scans

    # keep lowest likelihood for each r value
    min_scans = []
    for scan in scans:
        inds = np.lexsort((scan.df['dnll'], scan.df['mu']))
        scan = scan[inds]
        # remove duplicates
        _, uniq_inds = np.unique(scan.df['mu'], return_index=True)
        scan = scan[uniq_inds]
        min_scans.append(scan)

    return min_scans


def clean_scan(scan):
    dnll = scan.df['dnll']

    imin = np.argmin(dnll)

    # Filter left branch: Forbid positive slope
    keep_left = [imin]
    for i in range(imin, -1, -1):
        if dnll[i] < dnll[keep_left[-1]]: continue
        keep_left.append(i)

    # Filter left branch: Forbid negative slope
    keep_left = [imin]
    for i in range(imin, -1, -1):
        while dnll[i] < dnll[keep_left[-1]]:
            keep_left.pop()
        keep_left.append(i)

    # Filter right branch: Forbid negative slope
    keep_right = [imin]
    for i in range(imin+1, len(dnll)):
        while dnll[i] < dnll[keep_right[-1]]:
            keep_right.pop()
        keep_right.append(i)

    keep = np.array(keep_left[::-1] + keep_right[1:])
    logger.warning('Filtering out {} points due to bad slope'.format(len(dnll) - len(keep)))
    return scan[keep]


def apply_ranges(ax, do_x=True, do_y=True):
    if do_x:
        xmin = bsvj.read_arg('--xmin', type=float).xmin
        xmax = bsvj.read_arg('--xmax', type=float).xmax
        if xmax is not None: ax.set_xlim(right=xmax)
        if xmin is not None: ax.set_xlim(left=xmin)
    if do_y:
        ymin = bsvj.read_arg('--ymin', type=float).ymin
        ymax = bsvj.read_arg('--ymax', type=float).ymax
        if ymax is not None: ax.set_ylim(top=ymax)
        if ymin is not None: ax.set_ylim(bottom=ymin)

@scripter
def muscan():
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    correctminimum = bsvj.pull_arg('--correctminimum', action='store_true').correctminimum
    include_dots = bsvj.pull_arg('--include-dots', action='store_true').include_dots
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='muscan.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean

    with quick_ax(outfile=outfile) as ax:
        min_mu = 1e6
        max_mu = -1e6
        draw_str = '-'
        alpha=1.
        if clean:
            draw_str += '-'
            alpha=.2
        if include_dots: draw_str += 'o'

        for name,scan in get_scans(rootfiles, correctminimum).items():
            mu = scan.df['mu']
            dnll = scan.df['dnll']
            min_mu = min(min_mu, np.min(mu))
            max_mu = max(max_mu, np.max(mu))
            line = ax.plot(mu, dnll, draw_str, label=None if clean else name, alpha=alpha)[0]
            if clean:
                cln = clean_scan(scan)
                ax.plot(cln.df['mu'], cln.df['dnll'], label=name, color=line.get_color())

        ax.plot([min_mu, max_mu], [.0, .0], color='lightgray')
        ax.plot([min_mu, max_mu], [.5, .5], label='$1\sigma$')
        ax.plot([min_mu, max_mu], [1., 1.], label='$2\sigma$')
        ax.set_xlabel('$\mu$')
        ax.set_ylabel('$\Delta NLL$')
        apply_ranges(ax)
        ax.legend(framealpha=0.)

def get_scans(rootfiles, correct_minimum=False, keep_all=False):
    if isinstance(rootfiles, str): rootfiles = [rootfiles]
    scans = {}
    for rootfile in rootfiles:
        if not osp.exists(rootfile):
            logger.warning(f"get_scans: {rootfile} not found, skipping")
            continue
        name = name_from_combine_rootfile(rootfile)
        stmp = extract_scans(rootfile, correct_minimum, keep_all)
        if len(stmp)==1: scans[name] = stmp[0]
        else: logger.warning(f"{len(stmp)}!=1 scans found in {rootfile}, skipping")
    return scans

def plot_trackedparam(scans, param, outfile, clean, error=False):
    drwstr = '-'
    alpha = 1.
    if clean:
        drwstr += '-'
        alpha = .2
        if error: alpha = .6

    def error_band(ax, scan, param):
        ax.fill_between(scan.df['mu'], scan.df[param]-scan.df['error_'+param], scan.df[param]+scan.df['error_'+param], alpha=0.2)

    with quick_ax(outfile=outfile) as ax:
        for name,scan in scans.items():
            ax.plot(scan.df['mu'], scan.df[param], drwstr, label=name, alpha=alpha)
            if clean:
                cln = clean_scan(scan)
                ax.plot(cln.df['mu'], cln.df[param], label=name)
                if error: error_band(ax, cln, param)
            else:
                if error: error_band(ax, scan, param)

        ax.set_xlabel('$\mu$')
        ax.set_ylabel(param)
        apply_ranges(ax)
        ax.legend(framealpha=0.)

@scripter
def trackedparam():
    param = bsvj.pull_arg('param', type=str).param
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean
    error = bsvj.pull_arg('--error', action='store_true').error

    plot_trackedparam(get_scans(rootfiles), param, outfile, error)

@scripter
def trackedparams():
    # automatically plot any param whose error is also stored
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='{}.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean

    scans = get_scans(rootfiles)
    scan = scans[next(iter(scans))]
    pars_to_scan = [par for par in scan.df if ('error_'+par in scan.df or 'theta' in par) and par!='r']
    for par in pars_to_scan:
        logger.info(f"Plotting {par}")
        plot_trackedparam(scans, par, outfile.format(par), clean, error='error_'+par in scan.df)

@scripter
def debugparams():
    # plot params vs. likelihood for specified r value
    # expects input from MultiDimFit w/ --debugRandIteration option
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='debugparams.png').outfile
    r_debug = bsvj.read_arg('-r', '--mu', type=float).mu
    eps_debug = bsvj.read_arg('-e', '--eps', type=float, default=0.0001).eps
    sigma_debug = bsvj.read_arg('-s', '--sigma', type=float, default=4).sigma

    scans = get_scans(rootfiles, keep_all=True)
    scan = scans[next(iter(scans))]
    # only consider params with stored errors, but exclude r and auto-generated normalization params
    # i.e. only show params actually varied during fit
    pars_to_scan = [par for par in scan.df if 'error_'+par in scan.df and par!='r' and (par.startswith('bsvj') or 'bsvj' not in par)]

    # based on https://stackoverflow.com/a/60417813
    y_axis_dist = 0.2
    def plot_with_y_axis(scan, ax, param, idx, colors):
        canvas = ax
        # generate a new y axis and move to left side
        if idx>0:
            # advance cycle
            canvas = ax.twinx()
            canvas.spines['left'].set_position(("axes", -y_axis_dist*idx))
            canvas.set_frame_on(True)
            canvas.patch.set_visible(False)
            for nsp,sp in canvas.spines.items():
                if nsp!='left': sp.set_visible(False)
        else:
            canvas.set_xlabel('$\Delta NLL$')
        # todo: show error bars?
        pcolor = next(colors)
        canvas.scatter(scan.df['dnll'], scan.df[param], marker='o', facecolors='none', color=pcolor)
        apply_ranges(canvas, do_y=False)
        canvas.set_ylabel(param)
        canvas.yaxis.label.set_color(pcolor)
        canvas.tick_params(axis='y', colors=pcolor)
        canvas.yaxis.set_label_position('left')
        canvas.yaxis.set_ticks_position('left')

    for name,scan in scans.items():
        oname = outfile.replace(".png","_"+name.split()[0]+".png")
        # pick desired r value
        scan = scan[abs(scan.df['mu']-r_debug)<eps_debug]
        # remove extreme outliers
        scan = scan[abs(scan.df['dnll']-np.mean(scan.df['dnll']))<sigma_debug*np.std(scan.df['dnll'])]
        colors = get_color_cycle()
        with quick_ax(outfile=oname) as ax:
            for ipar,par in enumerate(pars_to_scan):
                plot_with_y_axis(scan, ax, par, ipar, colors)

def get_toy(file):
    toys = file.Get("toys")
    if toys==None: return None
    # pick the first RooDataSet
    for key in toys.GetListOfKeys():
        obj = toys.Get(key.GetName())
        if isinstance(obj, ROOT.RooDataSet):
            logger.info(f"Using toy {key.GetName()}")
            return obj
    return None

@scripter
def mtdist():
    import warnings
    warnings.filterwarnings('ignore', r'The value of the smallest subnormal')

    rootfile = bsvj.pull_arg('rootfile', type=str).rootfile
    only_sig = bsvj.pull_arg('--onlysig', action='store_true').onlysig
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='muscan.png').outfile
    bkg_names = bsvj.pull_arg('--bkg', type=str, default=['roomultipdf','bkg'], nargs='*').bkg
    sig_name = bsvj.pull_arg('--sig', type=str, default='sig').sig
    ch_name = bsvj.pull_arg('--channel', type=str, default='bsvj').channel
    sel_name = bsvj.pull_arg('--sel', type=str, default=None).sel
    title = bsvj.pull_arg('--title', type=str, default=None).title
    show_chi2 = bsvj.pull_arg('--chi2', action='store_true').chi2

    from scipy.interpolate import make_interp_spline # type:ignore

    print(rootfile)
    with bsvj.open_root(rootfile) as f:
        ws = bsvj.get_ws(f)
        toy = get_toy(f)

    mt = ws.var('mt')
    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[1:]+mt_binning[:-1])
    mt_bin_widths = mt_binning[1:] - mt_binning[:-1]

    mu_prefit = ws.var('r').getVal()

    # Get the data histogram
    data = ws.data('data_obs')
    data_label = 'Data'
    # check for a toy
    if toy is not None:
        data = toy
        data_label = 'Data (toy)'
    y_data = bsvj.roodataset_values(data,channel=ch_name)[1]

    # Get histogram from generated toy
    errs_data = np.sqrt(y_data)
    logger.info(f'Prefit data # entries = {y_data.sum():.2f}, should match with datacard')

    # helper to break out of RooMultiPdf
    def get_pdf(ws, pdf_name):
        pdf = ws.pdf(pdf_name)
        if hasattr(pdf,'getCurrentPdf'):
            pdf = pdf.getCurrentPdf()
        return pdf

    # check the background name/type
    for bkg_name in ['roomultipdf','bkg']:
        bkg_name_shape = f'shapeBkg_{bkg_name}_{ch_name}'
        bkg_name_norm = f'n_exp_final_bin{ch_name}_proc_{bkg_name}'
        bkg_pdf = get_pdf(ws,bkg_name_shape)
        if bkg_pdf:
            break

    # Get the prefit background histogram
    y_bkg_init = bsvj.pdf_values(bkg_pdf, mt_bin_centers)
    bkg_norm_init = ws.function(bkg_name_norm).getVal()
    y_bkg_init *= bkg_norm_init
    logger.info(f'Prefit bkg norm = {y_bkg_init.sum():.2f}, should match with datacard')

    # calculate chi square for prefit bkg fit
    _wrapped_prefit = bsvj.PDF()
    _wrapped_prefit.pdf = bkg_pdf
    _wrapped_prefit.n_pars = bkg_pdf.getParameters(data).getSize()
    print(_wrapped_prefit.n_pars)
    chi2_prefit_vf = bsvj.get_chi2_viaframe(mt, _wrapped_prefit, data)
    chi2_prefit = chi2_prefit_vf['chi2']
    ndf_prefit = chi2_prefit_vf['ndf']

    # signal info
    sig_name_shape = f'shapeSig_{ch_name}_{sig_name}'
    sig_name_norm = f'n_exp_bin{ch_name}_proc_{sig_name}'
    sig_name_norm_final = f'n_exp_final_bin{ch_name}_proc_{sig_name}'
    has_systematics = not(bool(ws.embeddedData(sig_name_shape)))
    logger.info(f'Datacard {"has" if has_systematics else "does not have"} systematics')
    if has_systematics: sig_name_shape = f'{sig_name_shape}_morph'

    # Get the pre-fit signal histogram
    if has_systematics:
        # Datacard with systematics
        # The signal histogram is saved only as a pdf, and reconstructing what
        # the signal should look like at mu=1, systs=0 is a little more complicated
        # Get the PDF and normalization separately
        # Temporarily ensure mu=1
        ws.var('r').setVal(1.0)
        sig = ws.pdf(sig_name_shape)
        norm_init = ws.function(sig_name_norm).getVal()
        y_sig = norm_init * bsvj.pdf_values(sig, mt_bin_centers)
        ws.var('r').setVal(mu_prefit)
    else:
        # Datacard without systematics: Just get the datahist
        sig = ws.embeddedData(sig_name_shape)
        y_sig = bsvj.roodataset_values(sig,channel=ch_name)[1]
    logger.info(f'Prefit signal norm = {y_sig.sum():.2f}, should match with datacard')

    # __________________________________
    # Load snapshot - everything is final fit values from this point onward
    snapshot = 'MultiDimFit'
    if 'FitDiagnostics' in rootfile: snapshot = 'fit_s'
    ws.loadSnapshot(snapshot)

    # Best fit mu value
    mu = ws.var('r').getVal()

    # Final-fit bkg
    bkg_pdf_final = get_pdf(ws,bkg_name_shape)
    y_bkg = bsvj.pdf_values(bkg_pdf_final, mt_bin_centers)
    bkg_norm = ws.function(bkg_name_norm).getVal()
    y_bkg *= bkg_norm
    logger.info(f'Initial bkg norm: {bkg_norm_init:.2f}; Final bkg norm: {bkg_norm:.2f}')

    # Compute bkg + mu * sig
    if has_systematics:
        # Use the shape pdf
        sig_final = ws.pdf(sig_name_shape)
        norm = ws.function(sig_name_norm_final).getVal()
        logger.info(f'Initial signal norm: {norm_init:.2f}; Postfit signal norm: {norm:.2f}')
        # mu should be already included for post fit signal, right?
        y_sig_postfit = norm * bsvj.pdf_values(sig_final, mt_bin_centers)
        y_sb = y_bkg + y_sig_postfit
    else:
        # No shape changes, just multiply signal by signal strength
        sig_final = sig
        y_sig_postfit = mu*y_sig
        y_sb = y_bkg + y_sig_postfit

    # calculate chi square for s+b
    _wrapped_postfit = bsvj.PDF()
    bkg_frac = ROOT.RooRealVar("bkg_frac", "bkg_frac", sum(y_bkg) / sum(y_sb))
    _wrapped_postfit.pdf = ROOT.RooAddPdf("postfit", "postfit", bkg_pdf_final, sig_final, bkg_frac)
    _wrapped_postfit.n_pars = bkg_pdf_final.getParameters(data).getSize() + 1
    chi2_postfit_vf = bsvj.get_chi2_viaframe(mt, _wrapped_postfit, data)
    chi2_sb = chi2_postfit_vf['chi2']
    ndf_sb = chi2_postfit_vf['ndf']

    fig, (ax, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16), sharex=True
        )
    plt.subplots_adjust(hspace=0.05)

    ax2.plot([mt_binning[0], mt_binning[-1]], [0,0], c='black')

    # petroff 6-color scheme
    colors = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
    cnames = ["blue", "orange", "red", "mauve", "gray", "purple"]
    petroff = {cname:color for cname,color in zip(cnames,colors)}

    if only_sig:
        ax.step(mt_binning[:-1], y_sig, where='post', c=petroff["orange"], label=r'$S_{prefit}$')
        ax.step(mt_binning[:-1], y_sig_postfit, where='post', c=petroff["red"], label=f'$S_{{fit}}$ ($\mu={mu:.3f}$)')
        ax.step(mt_binning[:-1], y_sig_postfit/mu, where='post', c=petroff["purple"], label=r'$S_{fit}$ ($\mu=1$)')
        ax2.plot([mt_binning[0], mt_binning[-1]], [1,1], c='black')
        ax2.step(mt_binning[:-1], y_sig_postfit/y_sig, where='post', c=petroff["red"])
        ax2.step(mt_binning[:-1], y_sig_postfit/y_sig/mu, where='post', c=petroff["purple"])
    else:
        ax.errorbar(
            mt_bin_centers, y_data,
            xerr=.5*mt_bin_widths, yerr=errs_data,
            fmt='o', c='black', label=data_label,
            )

        class RangeChecker():
            def __init__(self):
                self.ymin = -1
                self.ymax = 1
            def __call__(self, lines):
                line = lines[0]
                self.ymin = np.minimum(self.ymin, np.min(line.get_ydata()))
                self.ymax = np.maximum(self.ymax, np.max(line.get_ydata()))
        checker = RangeChecker()

        # set order here to get correct order in legend

        mt_fine = np.linspace(mt_binning[0], mt_binning[-1], 100) # For fine plotting
        spl = make_interp_spline(mt_bin_centers, y_bkg, k=3)  # type of this is BSpline
        y_bkg_fine = spl(mt_fine)
        ax.plot(mt_fine, y_bkg_fine, label=r'$B_{\mathrm{fit}}$', c=petroff["blue"])
        _ = ax2.step(mt_binning[:-1], (y_data - y_bkg) / np.sqrt(y_data), where='post', c=petroff["blue"])
        checker(_)

        ax.step(mt_binning[:-1], np.abs(y_sig_postfit), where='post', label=r'$S_{{\mathrm{{fit}}}}$ ($\mu_{{\mathrm{{fit}}}}={0:.2f}$)'.format(mu), c=petroff["red"])
        _ = ax2.step(mt_binning[:-1], y_sig_postfit / np.sqrt(y_data), where='post', c=petroff["red"])
        checker(_)

        ax.step(mt_binning[:-1], y_sb, where='post', c=petroff["mauve"], label=r'$B_{\mathrm{fit}}+S_{\mathrm{fit}}$'+(f' [{chi2_sb:.1f}/{ndf_sb}]' if show_chi2 else ''))
        _ = ax2.step(mt_binning[:-1], (y_data - y_sb) / np.sqrt(y_data), where='post', c=petroff["mauve"])
        checker(_)

        ax.step(mt_binning[:-1], y_bkg_init, where='post', c=petroff["gray"], linestyle='--', label=r'$B_{\mathrm{prefit}}$'+(f' [{chi2_prefit:.1f}/{ndf_prefit}]' if show_chi2 else ''))
        _ = ax2.step(mt_binning[:-1], (y_data - y_bkg_init) / np.sqrt(y_data), where='post', c=petroff["gray"], linestyle='--')
        checker(_)

        ax.step(mt_binning[:-1], y_sig, where='post', label=r'$S_{\mathrm{prefit}}$ ($\mu=1$)', c=petroff["orange"], linestyle='--')
        ax2.step(mt_binning[:-1], y_sig / np.sqrt(y_data), where='post', c=petroff["orange"], linestyle='--')
        # do not check range

    if title is None:
        title = name_from_combine_rootfile(rootfile, sel=(sel_name is None))
        if sel_name: title = '\n'.join([title,sel_name])
    leg = ax.legend(framealpha=0.0, fontsize=22, title=title, ncol=1 if only_sig else 2)
    leg._legend_box.align = "left"
    ax.set_ylabel('$N_{\mathrm{events}}$')
    ax2.set_xlabel(r'$m_{\mathrm{T}}$ [GeV]')
    ax.set_yscale('log')
    ax2.set_ylabel('(data - fit) / $\sqrt{\mathrm{data}}$')
    if only_sig: ax2.set_ylabel('postfit / prefit')

    # axis ranges
    if only_sig:
        ax.set_ylim(1., None)
    else:
        ax.set_ylim(np.power(10, np.floor(np.log10(np.min(y_data[y_data>0])))), None)
        ax2.set_ylim(1.1*checker.ymin, 1.1*checker.ymax)
    apply_ranges(ax)

    plt.savefig(outfile, bbox_inches='tight')
    if not(BATCH_MODE) and cmd_exists('imgcat'): os.system('imgcat ' + outfile)



def get_cls(obs, asimov):
    from scipy.stats import norm # type:ignore
    quantiles = np.array([0.025, 0.16, 0.50, 0.84, 0.975])

    # Keep only scan points where both obs and asimov have a mu
    keep_obs = np.isin(obs.df['mu'], asimov.df['mu'])
    keep_asimov = np.isin(asimov.df['mu'], obs.df['mu'])
    obs = obs[keep_obs]
    asimov = asimov[keep_asimov]

    # Filter out duplicates
    obs = obs[np.unique(obs.df['mu'], return_index=True)[1]]
    asimov = asimov[np.unique(asimov.df['mu'], return_index=True)[1]]

    np.testing.assert_array_equal(obs.df['mu'], asimov.df['mu'])

    # 1007.1727 Section 2.5: Alternative test statistic \tilde{q}_{\mu} for upper limits, Eq. (16)
    dnll_obs = obs.df['dnll']
    q_obs = []
    for i, mu in enumerate(obs.df['mu']):
        if mu < obs.bestfit.df['mu']:
            dnll_obs_min = np.min(dnll_obs[:i+1])
            dnll_obs_constrained = dnll_obs[i] - dnll_obs_min
        else:
            dnll_obs_constrained = dnll_obs[i]
        q_obs.append(2.*max(dnll_obs_constrained, 0.))
    q_obs = np.array(q_obs)
    assert q_obs.shape == (obs.n,)

    q_A = 2. * asimov.df['dnll']
    q_A[q_A < 0.] = 0.  # Set negative values to 0

    assert np.all(  ((q_obs >= 0.) & (q_obs <= q_A)) | (q_obs > q_A)  )

    # CL_{s} = p_{s+b}/p_{b}
    # from 1007.1727:
    # \sigma_{A}^2 = \mu^2/q_A : Eq. (31) (Section 3.2)
    # p_{\mu} = 1 - F(\tilde{q}_{\mu}|\mu) : Eq. (67) (Section 3.7)
    # below: take \tilde{q}_{\mu} -> q_obs
    # s+b: \mu^{\prime} = \mu
    #      CDF from Eq. (66) (Section 3.7), substituting Eq. (31)
    #      F(q_obs|\mu) = \Phi(\sqrt{q_obs}) for 0 < q_obs <= q_A
    #                   = \Phi((q_obs+q_A)/(2\sqrt{q_A})) for q_obs > q_A
    #   b: \mu^{\prime} = 0
    #      CDF from Eq. (65) (Section 3.7), substituting Eq. (31)
    #      F(q_obs|0)   = \Phi(\sqrt{q_obs} - \sqrt{q_A}) for 0 < q_obs <= q_A
    #                   = \Phi((q_obs - q_A)/(2\sqrt{q_A})) for q_obs > q_A
    sb = np.where(
        q_obs <= q_A,
        1. - norm.cdf( np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs+q_A) , np.sqrt(q_A)) )
        )
    b = np.where(
        q_obs <= q_A,
        norm.cdf( np.sqrt(q_A)-np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs-q_A) , np.sqrt(q_A)) )
        )
    s = sb / b

    # expected: take q_obs -> q_A, therefore always in q_obs <= q_A case
    # follow 1007.1727 Section 4.3 Eq. (89): take \mu -> \mu +/- N\sigma
	# \tilde{q}_{\mu} = (\mu - \hat{\mu})^2 / \sigma^2 : Eq. (62) (Section 3.7)
	# therefore \sqrt{q} -> \sqrt{q} +/- N
	# then use Eq. (66), (65) as above
	# s+b: F = \Phi(\sqrt{q} +/- N)
	#   b: F = \Phi(+/-N)
	# and N = ppf(q), so Phi(ppf(q)) = q
    s_exp = { q : (1.-norm.cdf(np.sqrt(q_A) - norm.ppf(q))) / q for q in quantiles}

    return bsvj.AttrDict(s=s, b=b, sb=sb, q_obs=q_obs, q_A=q_A, obs=obs, asimov=asimov, s_exp=s_exp)


def interpolate_95cl_limit(cls):
    mu = cls.obs.df['mu']
    def interpolate(cl, thing):
        select = ((cl < .99) & (cl > .001) & (mu>0))
        if select.sum() == 0:
            logger.error('0.01<cl<0.20 & mu>0 yields NO scan points; can\'t interpolate %s', thing)
            return None

        order = np.argsort(cl[select])
        try:
            if DEBUG:
                with quick_ax() as ax:
                    ax.set_title(thing)
                    ax.plot(cl[select][order], mu[select][order])
                    ax.plot([.05, .05], [min(mu[select][order]), max(mu[select][order])])
                    ax.set_xlabel('cl')
                    ax.set_ylabel('mu')
            res = np.interp(.05, cl[select][order], mu[select][order])
        except ValueError as e:
            logger.error('Interpolation failed for %s: %s', thing, e)
            res = None
        return res

    d = bsvj.AttrDict()
    d['twosigma_down'] = interpolate(cls.s_exp[0.975], 'twosigma_down')
    d['onesigma_down'] = interpolate(cls.s_exp[0.84], 'onesigma_down')
    d['expected'] = interpolate(cls.s_exp[0.5], 'expected')
    d['onesigma_up'] = interpolate(cls.s_exp[0.16], 'onesigma_up')
    d['twosigma_up'] = interpolate(cls.s_exp[0.025], 'twosigma_up')
    d['observed'] = interpolate(cls.s, 'observed')
    d['twosigma_success'] = (d['twosigma_down'] is not None) and (d['twosigma_up'] is not None)
    d['onesigma_success'] = (d['onesigma_down'] is not None) and (d['onesigma_up'] is not None)
    return d


def safe_divide(a, b):
    return np.divide(a, b, out=np.zeros_like(a), where=b!=0)


class LimitObj:
    def __init__(self):
        self.rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
        self.clean = bsvj.pull_arg('--clean', action='store_true').clean

        # compute cls curves and limits
        # sorted by signal parameters
        self.results = {}
        for observed, asimov in zip(*organize_rootfiles(self.rootfiles)):
            meta = svj.metadata_from_path(asimov)
            key = (meta['mz'], meta['mdark'], meta['rinv'])

            obs, asi = extract_scans([observed, asimov], correct_minimum=True)
            if self.clean:
                obs = clean_scan(obs)
                asi = clean_scan(asi)

            result = {}
            result['cls'] = get_cls(obs, asi)
            result['limit'] = interpolate_95cl_limit(result['cls'])
            result['observed'] = observed
            result['asimov'] = asimov
            self.results[key] = result


@scripter
def explim():
    limits = LimitObj()
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='explim.txt').outfile

    with open(outfile, 'w') as ofile:
        for key,result in limits.results.items():
            ofile.write(' '.join(str(x) for x in [key[0], key[1], str(key[2]).replace('.','p'), result['limit'].expected])+'\n')


@scripter
def cls():
    limits = LimitObj()
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile

    for key,result in limits.results.items():
        cls = result['cls']
        limit = result['limit']

        logger.info(
            'Limit {} {} {}: 2sd={:.4f} 1sd={:.4f} exp={:.4f} 1su={:.4f} 2su={:.4f} obs={:.4f}'
            .format(
                key[0], key[1], key[2],
                limit.twosigma_down,
                limit.onesigma_down,
                limit.expected,
                limit.onesigma_up,
                limit.twosigma_up,
                limit.observed
                )
            )

        with quick_ax(outfile=outfile) as ax:

            mu = cls.obs.df['mu']
            mu_best = cls.obs.bestfit.df['mu']

            ax.plot([], [], ' ', label=name_from_combine_rootfile(result['observed'], True))
            ax.plot([mu[0], mu[-1]], [.05, .05], label='95%', c='purple')
            ax.plot(mu, cls.s, label='s', c='black')
            ax.plot(mu, cls.b, label='b', c='blue')
            ax.plot(mu, cls.sb, label='sb', c='red')
            ax.plot(
                [mu_best, mu_best], [0., 1.05],
                c='#f88379', label=r'$\mu_{{best}}={:.2f}$'.format(mu_best), alpha=.8
                )

            # Expected
            ax.fill_between(mu, cls.s_exp[0.975], cls.s_exp[0.84], color=cms_yellow, alpha=0.25)
            ax.fill_between(mu, cls.s_exp[0.84], cls.s_exp[0.16], color=cms_green, alpha=0.25)
            ax.fill_between(mu, cls.s_exp[0.16], cls.s_exp[0.025], color=cms_yellow, alpha=0.25)
            ax.plot(mu, cls.s_exp[0.5], c='black', linestyle='--', label=r'$s_{exp}$')

            # Limit points
            s = 45
            if limit.twosigma_success:
                ax.scatter([limit.twosigma_down, limit.twosigma_up], [.05, .05], c='xkcd:dark yellow', s=s)
            if limit.onesigma_success:
                ax.scatter([limit.onesigma_down, limit.onesigma_up], [.05, .05], c=cms_green, s=s)
            if limit.expected is not None: ax.scatter([limit.expected], [.05], c='black', s=s)
            if limit.observed is not None: ax.scatter([limit.observed], [.05], c='black', s=s)

            ax.legend()
            ax.set_xlim(0.)
            ax.set_ylim(0., 1.05)
            apply_ranges(ax)
            ax.set_xlabel(r'$\mu$')
            ax.set_ylabel('CL')


@scripter
def brazil():
    limits = LimitObj()
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile

    points = []
    for key,result in limits.results.items():
        points.append(bsvj.AttrDict(
            mz = key[0],
            xsec = bsvj.get_xs(key[0]),
            limit = result['limit'],
            cls = result['cls']
        ))

    print(
        '{:<5s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} | {:>8s}'
        .format(
            'mz',
            '2s down',
            '1s down',
            'exp',
            '1s up',
            '2s up',
            'obs'
            )
        )

    def format(nr, w=8):
        if nr is not None:
            return '{:+{w}.3f}'.format(nr, w=w)
        else:
            return '{:>{w}s}'.format('err', w=w)

    for p in points:
        print(
            '{:<5.0f} {} {} {} {} {} | {}'
            .format(
                p.mz,
                format(p.limit.twosigma_down),
                format(p.limit.onesigma_down),
                format(p.limit.expected),
                format(p.limit.onesigma_up),
                format(p.limit.twosigma_up),
                format(p.limit.observed)
                )
            )

    fig = plt.figure(figsize=(12,10))
    ax = fig.gca()

    with quick_ax(figsize=(12,10), outfile=outfile) as ax:
        ax.plot([],[],label='95% CL upper limits (cut-based)',color='white')
        ax.plot([],[],label=r'$m_{\mathrm{dark}}=10$ GeV, $r_{\mathrm{inv}}=0.3$',color='white')

        ax.fill_between(
            [p.mz for p in points if p.limit.twosigma_success],
            [p.limit.twosigma_down*p.xsec for p in points if p.limit.twosigma_success],
            [p.limit.twosigma_up*p.xsec for p in points if p.limit.twosigma_success],
            color=cms_yellow, label='95% expected'
            )
        ax.fill_between(
            [p.mz for p in points if p.limit.onesigma_success],
            [p.limit.onesigma_down*p.xsec for p in points if p.limit.onesigma_success],
            [p.limit.onesigma_up*p.xsec for p in points if p.limit.onesigma_success],
            color=cms_green, label='68% expected'
            )

        ax.plot(
            [p.mz for p in points if p.limit.expected is not None],
            [p.limit.expected*p.xsec for p in points if p.limit.expected is not None],
            c='blue', linestyle='--', label='Median expected'
            )
        ax.plot(
            [p.mz for p in points if p.limit.observed is not None],
            [p.limit.observed*p.xsec for p in points if p.limit.observed is not None],
            c='black', linestyle='-', label='Observed',marker='o'
            )
        ax.plot(
            [p.mz for p in points if p.limit.observed is not None],
            [p.xsec for p in points if p.limit.observed is not None],
            c='magenta', linestyle='-',label='Theoretical'
        )
        #ax.text(1.5,0.7,'95% CL upper limits (cut-based)',fontsize=10)
        #ax.text(1.5,0.5, r'$m_{dark}$=10 GeV, $r_{inv}$=0.3',fontsize=10)
        ax.set_xlabel(r'$m_{\mathrm{Z}^{\prime}}$ [GeV]')
        ax.set_ylim(0.1,50)
        ax.grid(True)
        #ax.set_ylabel(r'$\mu$')
        ax.set_ylabel(r'$\sigma B$ [pb]')
        ax.set_yscale('log')
        apply_ranges(ax)
        ax.legend(framealpha=0.0)


@scripter
def bkgfit():
    """
    Bkg fit plots
    """
    jsons = bsvj.get_jsons()
    regions = bsvj.pull_arg('--regions', type=str, nargs='+').regions
    pdftype = bsvj.pull_arg('pdftype', type=str, choices=bsvj.known_pdfs()).pdftype
    linscale = bsvj.pull_arg('--lin', action='store_true').lin
    scipyonly = bsvj.pull_arg('--scipyonly', action='store_true').scipyonly
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    gof_type = bsvj.pull_arg('--gof-type', type=str, default='rss', choices=bsvj.choices('gof')).gof_type
    asimov = bsvj.pull_arg('--asimov', default=False, action="store_true").asimov

    input = bsvj.InputData(regions, "theta", **jsons, asimov=asimov)

    bin_centers = .5*(input.mt_array[:-1]+input.mt_array[1:])
    bin_width = input.regions[0].mt[1] - input.regions[0].mt[0]

    pdfs = bsvj.pdfs_factory(pdftype, input.regions[0].mtvar, input.regions[0].bkg_th1, name=pdftype)

    for pdf in pdfs:
        if scipyonly:
            pdf.res = bsvj.fit_scipy_robust(pdf.expression, pdf.th1, cache=None)
            # Fill in the fitted parameters
            for p, val in zip(pdf.parameters, pdf.res.x):
                # Make sure the newly fitted value is actually in range
                if val < p.getMin(): p.setMin(val - 0.1*abs(val))
                if val > p.getMax(): p.setMax(val + 0.1*abs(val))
                p.setVal(val)
        else:
            pdf.res = bsvj.fit(pdf)


    if not scipyonly:
        # Make sure pdfs are really fitted
        pdf = pdfs[0]
        res_par_set = pdf.res.floatParsFinal()
        np.testing.assert_almost_equal(
            [p.getVal() for p in pdf.parameters],
            [res_par_set.at(i).getVal() for i in range(pdf.n_pars)]
            )

    # Make sure evaluation makes sense
    y_pdf_eval = bsvj.eval_expression(pdf.expression, [bin_centers] + [p.getVal() for p in pdf.parameters])
    y_pdf_eval /= y_pdf_eval.sum()
    np.testing.assert_almost_equal(y_pdf_eval, pdf.evaluate(bin_centers), decimal=2)

    # Do the fisher test and mark the winner pdf
    gofs, n_bins = bsvj.gof_bkgfit(input.regions[0].mtvar, input.regions[0].data_datahist, pdfs, gof_type=gof_type)
    winner = bsvj.do_fisher_test(gofs, n_bins)
    pdfs[winner].is_winner = True

    bkg_vals = np.asarray([input.regions[0].bkg_th1.GetBinContent(i+1) for i in range(input.regions[0].bkg_th1.GetNbinsX())])
    bkg_errs = np.asarray([input.regions[0].bkg_th1.GetBinError(i+1) for i in range(input.regions[0].bkg_th1.GetNbinsX())])

    figure, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16), sharex=True)

    ax.plot([], [], ' ', label=f'{pdftype}, {input.regions[0].metadata["selection"]}')
    ax.step(input.regions[0].mt[:-1], bkg_vals, where='post', label=r'BKG', c='b')
    ax2.plot([input.regions[0].mt[0], input.regions[0].mt[-1]], [0.,0.], c='gray')

    fine_mt_axis = np.linspace(input.regions[0].mt[0], input.regions[0].mt[-1], 100)
    for pdf in pdfs:
        par_vals = [p.getVal() for p in pdf.parameters]

        y_pdf = pdf.evaluate(bin_centers)
        if abs(1. - y_pdf.sum()) > 0.01: logger.error('PDF norm is off from 1.:', y_pdf.sum())

        if getattr(pdf, 'is_winner', False):
            logger.warning('par vals: %s', par_vals)
            logger.warning('y_pdf pre norm: %s (norm=%s)', y_pdf, y_pdf.sum())

        y_pdf *= bkg_vals.sum()

        if getattr(pdf, 'is_winner', False):
            logger.warning('y_pdf post norm: %s (norm=%s)', y_pdf, y_pdf.sum())

        chi2_vf = bsvj.get_chi2_viaframe(input.regions[0].mtvar, pdf, input.regions[0].data_datahist)
        chi2 = chi2_vf['chi2']

        label = (
            '{}, $\\chi^{{2}}={:.5f}$: ['.format(pdf.n_pars, chi2)
            + ', '.join(['{:.2f}'.format(v) for v in par_vals]) + ']'
            )
        if getattr(pdf, 'is_winner', False): label += ' WINNER'

        y_pdf_fine = bsvj.eval_expression(pdf.expression, [fine_mt_axis] + par_vals)
        bin_scale = bin_width / (fine_mt_axis[1]-fine_mt_axis[0])
        y_pdf_fine = y_pdf_fine / y_pdf_fine.sum() * sum(bkg_vals) * bin_scale
        line = ax.plot(fine_mt_axis, y_pdf_fine, label=label)[0]

        pulls = (bkg_vals - y_pdf) / bkg_errs
        ax2.scatter(bin_centers, pulls, color=line.get_color())

    ax.legend(fontsize=18, framealpha=0.0)
    ax.set_ylabel('$N_{events}$')

    ax2.set_ylabel(r'(bkg - pdf) / $\Delta$bkg')
    ax2.set_xlabel(r'$m_{T}$ (GeV)')
    if not linscale: ax.set_yscale('log')
    plt.savefig(outfile, bbox_inches='tight')
    if not(BATCH_MODE) and cmd_exists('imgcat'): os.system('imgcat ' + outfile)

def get_objs(file_and_objs):
    fsplit = file_and_objs.split(':')
    fname = fsplit[0]
    f = ROOT.TFile.Open(fname)
    if len(fsplit[1:])==1: return f.Get(fsplit[1])
    else: return [f.Get(oname) for oname in fsplit[1:]]

# common operations to evaluate tf fit
def get_tf_fit(fitresult, tf_name, tf_th1, mtscaled, bkg_eff=1.0, basis='Bernstein'):
    import rhalphalib as rl
    npar = len([f for f in fitresult.floatParsFinal() if tf_name in f.GetName()])-1
    tf_fn = rl.BasisPoly(tf_name, (npar,), ["mt"], basis=basis)
    tf_fn.update_from_roofit(fitresult)
    tf_fn_vals, tf_fn_band = tf_fn(mtscaled, nominal=True, errorband=True)
    # multiply bkg_eff into final values
    tf_fn_vals = bkg_eff * tf_fn_vals
    tf_fn_band = (bkg_eff * tf_fn_band[0], bkg_eff * tf_fn_band[1])
    chi2 = bsvj.get_tf_chi2(tf_th1, tf_fn_vals)
    ndf = len(tf_fn_vals) - npar - 1
    return {'tf_fn': tf_fn, 'tf_fn_vals': tf_fn_vals, 'tf_fn_band': tf_fn_band, 'chi2': chi2, 'ndf': ndf, 'npar': npar, 'fitresult': fitresult}

# replace in string starting from last match
def rreplace(s, old, new, count=1):
    return (s[::-1].replace(old[::-1], new[::-1], count))[::-1]

# plot a given tf and (optional) fit
def plot_tf(outfile, mt, tf, fit=None, title="", label="MC", ylabel="TF", suff="", canvas=None):
    if suff:
        outfile = rreplace(outfile,'.',f'_{suff}.',1)

    colors = get_color_cycle()
    if canvas is None: # Ploting existing items
        figure, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16), sharex=True)
        pcolor = next(colors)
        ax.errorbar(mt['pts'], tf['arr']['vals'], yerr=tf['arr']['errs'], label=label, color=pcolor)
        ax.set_ylabel(ylabel)
        xlabel = r'$m_{\mathrm{T}}$ [GeV]'
        pcolor = next(colors)
        ax2.set_ylabel(r'(TF - fit) / $\Delta$TF')
        ax2.set_xlabel(xlabel)
    else:
        figure, (ax, ax2) = canvas
        _ = next(colors)
        _ = next(colors)
        pcolor = next(colors) # Moving the color label
        print(outfile, "Updating canvas")

    ax.plot(mt['pts'], fit['tf_fn_vals'], label=f"$fit^{{{suff}}}$ ($\\mathrm{{n}} = {fit['npar']+1}$, $\\chi^2/\\mathrm{{ndf}} = {fit['chi2']:.1f}/{fit['ndf']}$)", color=pcolor)
    ax.fill_between(mt['pts'], fit['tf_fn_band'][0], fit['tf_fn_band'][1], alpha=0.2, color=pcolor)
    leg_args = {'fontsize': 18, 'framealpha': 0.0}
    if title: leg_args['title'] = title
    ax.legend(**leg_args)
    # pulls in lower panel
    pulls = (tf['arr']['vals'] - fit['tf_fn_vals']) / tf['arr']['errs']
    ax2.plot(mt['range'], [0.,0.], c='gray')
    ax2.scatter(mt['pts'], pulls, color=pcolor)
    apply_ranges(ax)
    figure.savefig(outfile, bbox_inches='tight')
    return figure, (ax, ax2) # Returning the plot containers so that it can be updated

@scripter
def bkgtf():
    """
    Bkg tf plots
    """
    jsons = bsvj.get_jsons()
    regions = bsvj.pull_arg('--regions', type=str, nargs='+').regions
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='tf.png').outfile
    fit_mc_file = bsvj.read_arg('--fit-mc', type=str, default=None).fit_mc
    fit_data_file = bsvj.read_arg('--fit-data', type=str, default=None).fit_data
    basis = bsvj.pull_arg('--basis', default='Bernstein').basis
    basis_mc = bsvj.pull_arg('--basis-mc', default='Bernstein').basis_mc
    asimov = bsvj.pull_arg('--asimov', default=False, action="store_true").asimov
    verbose = bsvj.pull_arg('-v','--verbose', default=False, action="store_true").verbose
    closure = bsvj.pull_arg('--closure', default=False, action="store_true").closure

    # input histograms always required: used to get x-axis
    title = name_from_combine_rootfile(jsons['sigfiles'][0])
    input = bsvj.InputData(regions, "rhalpha", **jsons, asimov=asimov)
    mt = {}
    mt['pts'] = input.mt_array[:-1] + 0.5*np.diff(input.mt_array)
    mt['scaled'] = (mt['pts'] - min(mt['pts']))/(max(mt['pts']) - min(mt['pts']))
    mt['range'] = [input.mt_array[0], input.mt_array[-1]]

    # TF from MC
    tf_mc = {}
    tf_mc['bkg_eff'] = input.regions[0].bkg_datahist.sum(False) / input.regions[1].bkg_datahist.sum(False)
    if verbose: print('mc_bkg_eff', tf_mc['bkg_eff'])
    # use ROOT to propagate errors when dividing
    tf_mc['th1'] = bsvj.get_tf_th1(input.regions)
    tf_mc['arr'] = bsvj.th1_to_hist(tf_mc['th1'])
    if verbose: print('tf_mc_th1', tf_mc['arr']['vals'].tolist())

    # polynomial fit from workspace
    fit_mc = None
    if fit_mc_file is not None:
        fitresult_mc = get_objs(fit_mc_file)
        fit_mc = get_tf_fit(fitresult_mc, 'tf_mc', tf_mc['th1'], mt['scaled'], tf_mc['bkg_eff'], basis=basis_mc)
        if verbose: print('fit_mc', fit_mc['tf_fn_vals'].tolist())
        if verbose: print('chi2_mc', fit_mc['chi2'], fit_mc['ndf'])

    # plot TF from MC
    mc_canvas = plot_tf(outfile, mt, tf_mc, fit_mc, ylabel=f'$TF_{{\\mathrm{{MC}}}}$ ({regions[0]} / {regions[1]})', suff='mc', title=title)

    # TF from data: everything comes from postfit file
    fit_data = None
    if fit_data_file is not None:
        fitresult_data, ws_data = get_objs(fit_data_file)
        # postfit shapes from workspace, based on mtdist()
        snapshot = 'MultiDimFit'
        if 'FitDiagnostics' in fit_data_file: snapshot = 'fit_s'
        ws_data.loadSnapshot(snapshot)
        postfit_regions = []
        for region in input.regions:
            channel = region.bin_name
            if closure:
                # closure/sanity check: postfit bkg ratio should be exactly equal to TF by construction
                bkg_name_shape = f'shapeBkg_bkg_{channel}'
                bkg_name_norm = f'n_exp_final_bin{channel}_proc_bkg'
                bkg_pdf = ws_data.pdf(bkg_name_shape)
                bkg_norm = ws_data.function(bkg_name_norm).getVal()
                bkg_hist = {
                    'vals': bkg_norm * bsvj.pdf_values(bkg_pdf, mt['pts']),
                    'errs': bkg_norm * bsvj.pdf_errors(bkg_pdf, fitresult_data, mt['pts']),
                    'binning': input.mt_array,
                    'metadata': {},
                }
                bkg_th1 = bsvj.Histogram(bkg_hist).th1(channel)
                postfit_regions.append(bkg_th1)
            else:
                # data - r*sig
                data_set = ws_data.data('data_obs')
                data_values = bsvj.roodataset_values(data_set,channel=channel,vars=region.mtvar)
                data_th1 = bsvj.Histogram({
                    'vals': data_values[1],
                    'errs': [bsvj.PoissonErrorUp(v) for v in data_values[1]], # not automatically populated in roodataset
                    'binning': input.mt_array,
                    'metadata': {}
                }).th1(f'data_{channel}')
                sig_name_shape = f'shapeSig_{channel}_sig_morph'
                sig_name_norm_final = f'n_exp_final_bin{channel}_proc_sig'
                mu = ws_data.var('r').getVal()
                sig_pdf = ws_data.pdf(sig_name_shape)
                sig_norm = ws_data.function(sig_name_norm_final).getVal()
                sig_th1 = bsvj.Histogram({
                    'vals': sig_norm * bsvj.pdf_values(sig_pdf, mt['pts']),
                    'errs': sig_norm * bsvj.pdf_errors(sig_pdf, fitresult_data, mt['pts']), # returns 0, so errors not propagated
                    'binning': input.mt_array,
                    'metadata': {},
                }).th1(f'sig_{channel}')
                data_th1.Add(sig_th1,-1)
                postfit_regions.append(data_th1)
        tf_data = {}
        tf_data['bkg_eff'] = postfit_regions[0].Integral() / postfit_regions[1].Integral()
        if verbose: print('data_bkg_eff', tf_data['bkg_eff'])
        tf_data['th1'] = bsvj.get_tf_th1(postfit_regions)
        label_data = r"$B_{\mathrm{fit}}" if closure else "$\mathrm{Data} - S_{\mathrm{fit}}$ "+"($\mu_{{\mathrm{{fit}}}}={0:.2f}$)".format(mu)

        # comparison of FD and MDF
        if verbose:
            print("par fitresult(FD) workspace(MDF)")
            print("\n".join([f"{f.GetName()} {f.getVal()} {ws_data.function(f.GetName()).getVal()}" for f in fitresult_data.floatParsFinal() if 'tf' in f.GetName()]))

        # if MC TF also used, divide it out first
        fit_mc_nuis = np.array([f.getVal() for f in fitresult_data.floatParsFinal() if 'tf_mc' in f.GetName()])
        if verbose: print('fit_mc_nuis',fit_mc_nuis.tolist())
        if fit_mc and len(fit_mc_nuis)>0:
            # reconstruct MC TF fit from decorrelated parameters
            paramfile = fit_mc_file.split(':')[0].replace("/bkgfit_","/mctf_").replace(".root",".npy")
            fit_mc_nominal = np.load(paramfile)
            if verbose: print('fit_mc_nominal',fit_mc_nominal.tolist())
            decofile = paramfile.replace("/mctf_","/deco_")
            decoVector = np.load(decofile)
            if verbose: print('decoVector',decoVector.tolist())
            fit_mc_parvalues = np.full(fit_mc_nominal.shape, None)
            for i in range(fit_mc_nominal.size):
                coef = decoVector[:, i]
                order = np.argsort(np.abs(coef))
                fit_mc_parvalues[i] = np.sum(coef[order] * fit_mc_nuis[order]) + fit_mc_nominal[i]
            if verbose: print('fit_mc_parvalues',fit_mc_parvalues.tolist())
            fit_mc['tf_fn'].set_parvalues(fit_mc_parvalues)
            fit_mc_vals = tf_mc['bkg_eff'] * fit_mc['tf_fn'](mt['scaled'], nominal=True)
            if verbose: print('fit_mc_vals', fit_mc_vals.tolist())

            # plot combined TF without dividing out MC
            # todo: propagate all uncertainties (currently just data TF uncertainties)
            tf_comb = {}
            tf_comb['bkg_eff'] = fit_mc_vals # bkg_eff multiplies tf_fn_vals in get_tf_fit(), so include entire MC TF in this case
            tf_comb['th1'] = tf_data['th1'].Clone()
            tf_comb['arr'] = bsvj.th1_to_hist(tf_comb['th1'])
            if verbose: print('tf_comb_th1', tf_comb['arr']['vals'].tolist())
            fit_comb = get_tf_fit(fitresult_data, 'tf_data', tf_comb['th1'], mt['scaled'], tf_comb['bkg_eff'], basis=basis)
            if verbose: print('fit_comb', fit_comb['tf_fn_vals'].tolist())
            if verbose: print('chi2_comb', fit_comb['chi2'], fit_comb['ndf'])
            plot_tf(outfile, mt, tf_comb, fit_comb, ylabel=f'$\\mathrm{{TF}}_{{\\mathrm{{comb}}}}$ ({regions[0]} / {regions[1]})', suff='comb', label=label_data, title=title)

            # divide out values
            for i in range(tf_data['th1'].GetNbinsX()):
                tf_data['th1'].SetBinContent(i+1, tf_data['th1'].GetBinContent(i+1)/fit_mc_vals[i])
                tf_data['th1'].SetBinError(i+1, tf_data['th1'].GetBinError(i+1)/fit_mc_vals[i])
            if verbose: print('tf_data_res_th1', bsvj.th1_to_hist(tf_data['th1'])['vals'].tolist())
            # bkg_eff already included in MC TF
            tf_data['bkg_eff'] = 1.0
            suff_data = 'data_res'

            # Postfit MC-only TF w/ uncertainties
            tf_post = {} # Creating the new ite for plotting
            tf_post['bkg_eff'] = fit_mc_vals
            tf_post['th1'] = tf_mc['th1'].Clone() # Key distinction compared with above
            tf_post['arr'] = bsvj.th1_to_hist(tf_post['th1'])
            if verbose: print('tf_post_th1', tf_post['arr']['vals'].tolist())
            fit_post = get_tf_fit(fitresult_mc, 'tf_mc', tf_post['th1'], mt['scaled'], tf_post['bkg_eff'], basis=basis_mc)
            plot_tf(outfile, mt, tf_post, fit_post, ylabel=f'$TF_{{\\mathrm{{MC}}}}^{{\\mathrm{{postfit}}}}$ ({regions[0]} / {regions[1]})', suff='mcpost', title=title, canvas=mc_canvas)
            # plot_tf(outfile, mt, tf_post, fit_post, ylabel=f'$TF_{{\\mathrm{{MC}}}}^{{\\mathrm{{postfit}}}}$ ({regions[0]} / {regions[1]})', suff='mcpost', title=title)
        else:
            if verbose: print('tf_data_th1', bsvj.th1_to_hist(tf_data['th1'])['vals'].tolist())
            tf_data['bkg_eff'] = tf_mc['bkg_eff']
            suff_data = 'data'

        tf_data['arr'] = bsvj.th1_to_hist(tf_data['th1'])
        fit_data = get_tf_fit(fitresult_data, 'tf_data', tf_data['th1'], mt['scaled'], tf_data['bkg_eff'], basis=basis)
        if verbose: print('fit_data', fit_data['tf_fn_vals'].tolist())
        if verbose: print('chi2_data', fit_data['chi2'], fit_data['ndf'])

        escape = lambda x: x.replace('_','\\_')
        plot_tf(outfile, mt, tf_data, fit_data, ylabel=f'$\\mathrm{{TF}}_{{\\mathrm{{{escape(suff_data)}}}}}$ ({regions[0]} / {regions[1]})', suff=suff_data, label=label_data, title=title)


@scripter
def ftest_toys():
    ftest_dump = bsvj.pull_arg('--results_dump', type=str).results_dump
    outpre = bsvj.pull_arg('-o', '--outpre', type=str).outpre
    dump = imp.load_source('ftest_dump', ftest_dump)
    winner = dump.winner
    nbins = dump.nbins
    results = dump.results

    # Running the toys items
    range_max = len(results)
    range_i = list(range(range_max))
    range_j = lambda i: list(range(i+1, range_max+1))
    range_i, range_j, range_max = bsvj.range_toys(results)

    # Extracting the information
    ftest_toys = {}
    ftest_data = {}
    ftest_pval = {}
    for i in range_i:
        for j in range_j:
            if (i,j) not in results: continue
            results_dict = bsvj.extract_results_toys(results,i,j)
            n1, n2 = results_dict["n1"], results_dict["n2"]
            toys1, toys2 = results_dict["gof1"]["toys"], results_dict["gof2"]["toys"]
            data1, data2 = results_dict["gof1"]["data"][0], results_dict["gof2"]["data"][0]
            ftest_toys[(n1,n2)] = [bsvj.fisher_metric(toys1[x], toys2[x], n1, n2, nbins) for x in toys1.keys()]
            ftest_data[(n1,n2)] = bsvj.fisher_metric(data1, data2, n1, n2, nbins)
            ftest_pval[(n1,n2)] = bsvj.compute_fisher_toys(results_dict["gof1"], results_dict["gof2"], n1, n2, nbins)

    winner = None
    for (n1,n2) in ftest_toys.keys(): # Creating the per-parameter comparion f-test plots
        f_toys = np.array(ftest_toys[(n1,n2)])
        f_toys = f_toys[f_toys > 0] # Only plotting stuff larger than 0
        f_data = ftest_data[(n1, n2)]
        p_val = ftest_pval[(n1,n2)]
        if p_val > 0.05 and winner is None:
            winner = n1
        outfile = f"{outpre}_fstat-{n1}vs{n2}.png"
        with quick_ax(outfile=outfile) as ax:
            h_val, bins = np.histogram(f_toys, bins=40) # Getting the bin values required to normalized the F-distribution plot
            ax.hist(f_toys, bins=40, histtype='step', label="Toys") # Tooys results
            x = np.linspace(np.min(f_toys), np.max(f_toys), 500) # F-distribution
            f = np.array([ROOT.TMath.FDist(_x, n2-n1, nbins-n2) for _x in x]) * (len(f_toys) * (bins[1]-bins[0]))
            ax.plot(x,f,color='r', label=f"F-dist, ndf=({n2-n1}, {nbins-n2})")
            ax.vlines([f_data], ymin=0, ymax=3, color='b' )
            ax.plot([],[],color='b', label=f'Observed ({f_data:.4f})')
            ax.plot([],[],color='none', label=f"p-value: {p_val:.4f}")
            ax.set_xlabel("F-test statistics")
            ax.set_ylabel("Number of toys")
            ax.set_ylim(top=np.max(h_val)*1.5)
            ax.legend(title=f"F-test for fit orders: ({n1}, {n2})")

    # Plotting the files for the winner evaluation
    for i,j in [(i,j) for i in range_i for j in range_j]:
        if (i,j) not in results: continue
        results_dict = bsvj.extract_results_toys(results,i,j)
        n1 = results_dict["n1"]
        n2 = results_dict["n2"]
        toys1, data1 = list(results_dict["gof1"]["toys"].values()), results_dict["gof1"]["data"][0]
        toys2 = list(results_dict["gof2"]["toys"].values())
        outfile = f"{outpre}_gof-npar{n1}.png"
        with quick_ax(outfile=outfile) as ax:
            h_val, bins = np.histogram(toys1, bins=40) # Getting the bin values required to normalized the F-distribution plot
            ax.hist(toys1, bins=40, histtype='step', label=f"Toys (n = {n1})") # Toys results
            ax.hist(toys2, bins=40, histtype='step', label=f"Toys (n = {n2})") # Toys results
            ax.vlines([data1], ymin=0, ymax=3, color='b' )
            ax.plot([],[],color='b', label=f'Observed = {data1:.4f} (n={n1})')
            ax.set_xlabel("Goodness of fit")
            ax.set_ylabel("Number of toys")
            ax.set_ylim(top=np.max(h_val)*1.5)
            ax.legend(title=f"{n1}-parameter fit")


@scripter
def ftest_scan():
    ftest_dir= bsvj.pull_arg('--results_dir', type=str).results_dir
    sel = bsvj.pull_arg('--sel', type=str).sel
    signals = bsvj.pull_arg("--signals", dest="signals", type=str, default="").signals
    outdir = bsvj.pull_arg('-o', '--outdir', type=str).outdir

    with open(signals,'r') as sfile:
        signals = [rhalph.Signal(*line.split(), 0) for line in sfile]
        ftest_dump_list = [f'{ftest_dir}/{rhalph.get_signame(s)}_sel-{sel}_mt_smooth_ftest-results.py' for s in signals]

    # Aggregarating the result into a signal file
    result = {
        (sig.mMed, sig.mDark, sig.rinv): imp.load_source('ftest_dump', ftest_dump).winner
        for sig, ftest_dump in zip(signals, ftest_dump_list)
        if os.path.exists(ftest_dump)
    }
    # Scanning verse mp
    for mDark in set(sig[1] for sig in result.keys()):
        with quick_ax(outfile=f"{outdir}/{sel}_ftest_scan_vs_mMed_mDark={mDark}.png") as ax:
            for rinv in sorted(set(sig[2] for sig in result.keys())):
                plot_points = np.array([(float(sig[0]), npar) for sig, npar in result.items() if sig[1]==mDark and sig[2] == rinv])
                if(len(plot_points) == 0): continue
                shift = float(rinv.replace('p', '.')) * 0.2
                mMed = plot_points[:,0]
                npar = plot_points[:, 1] + shift
                ax.plot(mMed, npar, marker='o', label="$r_{inv}$ = " + rinv.replace("p", "."))
            ax.set_ylabel('Chosen number of parameters')
            ax.set_xlabel("$m_{X}$ [GeV]")
            ax.legend(title="$m_{dark}$ = " + mDark + " GeV")


def plot_hist(th1, ax, **kwargs):
    hist = bsvj.th1_to_hist(th1)
    def get_kwargs(orig, keys):
        return {k : orig.get(k, None) for k in keys}
    step_keys = ['where','label','color']
    ax.step(hist['binning'][:-1], hist['vals'], **get_kwargs(kwargs, step_keys))
    fill_keys = ['where','alpha','color']
    fill_dict = get_kwargs(kwargs, fill_keys)
    fill_dict['step'] = fill_dict.pop('where')
    ax.fill_between(hist['binning'][:-1], hist['vals']-hist['errs'], hist['vals']+hist['errs'], **fill_dict)

@scripter
def bkgsrcr():
    """
    Bkg SR vs. CR plots
    """
    jsons = bsvj.get_jsons()
    regions = bsvj.pull_arg('--regions', type=str, nargs='+').regions
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='tf.png').outfile
    asimov = bsvj.pull_arg('--asimov', default=False, action="store_true").asimov

    input = bsvj.InputData(regions, "rhalpha", **jsons, asimov=asimov)

    with quick_ax(outfile=outfile) as ax:
        colors = get_color_cycle()
        pcolor = next(colors)
        plot_hist(input.regions[1].bkg_th1, ax, where='post', label=regions[1], alpha=0.2, color=pcolor)
        pcolor = next(colors)
        plot_hist(input.regions[0].bkg_th1, ax, where='post', label=regions[0], alpha=0.2, color=pcolor)
        ax.legend(fontsize=18, framealpha=0.0)
        ax.set_xlabel(r'$m_{\mathrm{T}}$ [GeV]')
        ax.set_ylabel(f'Number of events')
        ax.set_yscale('log')
        apply_ranges(ax)

@scripter
def hist():
    """
    Quickly plot histograms from a file
    """
    infile = bsvj.pull_arg('infile', type=str).infile
    histograms = bsvj.pull_arg('histograms', type=str, nargs='*').histograms

    if infile.endswith('.json'):
        with open(infile, 'r') as f:
            input = json.load(f, cls=bsvj.Decoder)

        if len(histograms) == 0:
            bsvj.ls_inputdata(input)

        else:
            with quick_ax() as ax:
                for pat in histograms:
                    for h in fnmatch.filter(input.keys(), pat):
                        hist = input[h]
                        ax.step(hist.binning[:-1], hist.vals, where='post', label=h)
                ax.legend()
                ax.set_yscale('log')
                ax.set_xlabel('MT (GeV)')

    elif infile.endswith('.root'):
        with bsvj.open_root(infile) as tf:
            w = tf.Get('SVJ')
            if w:
                if len(histograms)==0:
                    w.Print()
                else:
                    with quick_ax() as ax:
                        for h in histograms:
                            rdh = w.data(h)
                            x, y, dy = bsvj.roodataset_values(rdh)
                            ax.step(x, y, where='post', label=h)
                        ax.legend()
                        ax.set_yscale('log')
                        ax.set_xlabel('MT (GeV)')


if __name__ == '__main__':
    batch_mode(bsvj.pull_arg('-b', '--batch', action='store_true').batch)
    debug(bsvj.pull_arg('-d', '--debug', action='store_true').debug)
    fontsize = bsvj.read_arg('--fontsize', type=int, default=18).fontsize
    set_mpl_fontsize(legend=fontsize)
    scripter.run()
