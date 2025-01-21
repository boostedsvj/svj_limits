"""
Some scripts to quickly plot basic outcomes from combine scans
"""
from __future__ import print_function
import ROOT # type:ignore
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
from time import strftime
import argparse

# Add the directory of this file to the path so the boosted tools can be imported
import sys, os, os.path as osp, pprint, re, traceback, copy, fnmatch, shutil
from contextlib import contextmanager
sys.path.append(osp.dirname(osp.abspath(__file__)))
import svj_ntuple_processing as svj
import boosted_fits as bsvj
logger = bsvj.setup_logger('quickplot')

import numpy as np
import matplotlib as mpl
mpl.use('Agg') # in order to run in background / no-graphics environments
import matplotlib.pyplot as plt # type:ignore

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


def name_from_combine_rootfile(rootfile, strip_obs_asimov=False):
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
    pars_to_scan = [par for par in scan.df if 'error_'+par in scan.df and par!='r']
    for par in pars_to_scan:
        logger.info(f"Plotting {par}")
        plot_trackedparam(scans, par, outfile.format(par), clean, error=True)

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
        from itertools import cycle
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = cycle(prop_cycle.by_key()['color'])
        with quick_ax(outfile=oname) as ax:
            for ipar,par in enumerate(pars_to_scan):
                plot_with_y_axis(scan, ax, par, ipar, colors)

@scripter
def mtdist():
    rootfile = bsvj.pull_arg('rootfile', type=str).rootfile
    only_sig = bsvj.pull_arg('--onlysig', action='store_true').onlysig
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='muscan.png').outfile
    #toyrootfile = bsvj.pull_arg('--toyrootfile', type=str).toyrootfile

    from scipy.interpolate import make_interp_spline # type:ignore

    with bsvj.open_root(rootfile) as f:
        ws = bsvj.get_ws(f)

    mt = ws.var('mt')
    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[1:]+mt_binning[:-1])
    mt_bin_widths = mt_binning[1:] - mt_binning[:-1]

    mu_prefit = ws.var('r').getVal()

    # Get the prefit background histogram
    y_bkg_init = bsvj.pdf_values(ws.pdf('shapeBkg_roomultipdf_bsvj'), mt_bin_centers)
    pdf_raw_norm_prefit = np.sum(y_bkg_init)
    bkg_norm_init = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()
    y_bkg_init *= bkg_norm_init
    logger.info(f'Prefit bkg norm = {y_bkg_init.sum():.2f}, should match with datacard')

    has_systematics = not(bool(ws.embeddedData('shapeSig_sig_bsvj')))
    logger.info(f'Datacard {"has" if has_systematics else "does not have"} systematics')

    # Get the pre-fit signal histogram
    if has_systematics:
        # Datacard with systematics
        # The signal histogram is saved only as a pdf, and reconstructing what
        # the signal should look like at mu=1, systs=0 is a little more complicated
        # Get the PDF and normalization separately
        # Temporarily ensure mu=1
        ws.var('r').setVal(1.0)
        sig = ws.pdf('shapeSig_bsvj_sig_morph')
        norm_init = ws.function('n_exp_binbsvj_proc_sig').getVal()
        y_sig = norm_init * bsvj.pdf_values(sig, mt_bin_centers)
        ws.var('r').setVal(mu_prefit)
    else:
        # Datacard without systematics: Just get the datahist
        sig = ws.embeddedData('shapeSig_sig_bsvj')
        y_sig = bsvj.roodataset_values(sig)[1]
    logger.info(f'Prefit signal norm = {y_sig.sum():.2f}, should match with datacard')

    # Get the data histogram
    data = ws.data('data_obs')
    y_data = bsvj.roodataset_values(data)[1]

    # Get histogram from generated toy
    #with bsvj.open_root(toyrootfile) as f:
    #  toy = f.Get("toys/toy_1")
    #data = ROOT.RooDataSet(toy,'mt')
    #y_data = bsvj.roodataset_values(data)[1]
    errs_data = np.sqrt(y_data)
    logger.info(f'Prefit data # entries = {y_data.sum():.2f}, should match with datacard')

    # __________________________________
    # Load snapshot - everything is final fit values from this point onward
    ws.loadSnapshot('MultiDimFit')

    # Best fit mu value
    mu = ws.var('r').getVal()

    # Final-fit bkg
    bkg = ws.pdf('shapeBkg_roomultipdf_bsvj')
    y_bkg = bsvj.pdf_values(bkg, mt_bin_centers)
    logger.warning('y_bkg_postfit: %s', y_bkg)
    pdf_raw_norm_postfit = np.sum(y_bkg)
    bkg_norm = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()
    y_bkg *= bkg_norm
    logger.info(f'Initial bkg norm: {bkg_norm_init:.2f}; Final bkg norm: {bkg_norm:.2f}')

    # Compute bkg + mu * sig
    if has_systematics:
        # Use the shape pdf
        sig = ws.pdf('shapeSig_bsvj_sig_morph')
        norm = ws.function('n_exp_final_binbsvj_proc_sig').getVal()
        logger.info(f'Initial signal norm: {norm_init:.2f}; Postfit signal norm: {norm:.2f}')
        # mu should be already included for post fit signal, right?
        y_sig_postfit = norm * bsvj.pdf_values(sig, mt_bin_centers)
        y_sb = y_bkg + y_sig_postfit
    else:
        # No shape changes, just multiply signal by signal strength
        y_sig_postfit = mu*y_sig
        y_sb = y_bkg + y_sig_postfit

    fig, (ax, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16)
        )

    ax.plot([], [], ' ', label=name_from_combine_rootfile(rootfile))
    ax2.plot([mt_binning[0], mt_binning[-1]], [0,0], c='gray')

    if only_sig:
        ax.step(mt_binning[:-1], y_sig, where='post', c='purple', label=r'$S_{prefit}$')
        ax.step(mt_binning[:-1], y_sig_postfit, where='post', c='cyan', label=f'$S_{{postfit}}$ ($\mu={mu:.3f}$)')
        ax.step(mt_binning[:-1], y_sig_postfit/mu, where='post', c='red', label=r'$S_{postfit}$ ($\mu=1$)')
        ax2.plot([mt_binning[0], mt_binning[-1]], [1,1], c='purple')
        ax2.step(mt_binning[:-1], y_sig_postfit/y_sig, where='post', c='cyan')
        ax2.step(mt_binning[:-1], y_sig_postfit/y_sig/mu, where='post', c='red')
    else:
        ax.errorbar(
            mt_bin_centers, y_data,
            xerr=.5*mt_bin_widths, yerr=errs_data,
            fmt='o', c='black', label='Data'
            )
        # logger.warning('data (roodst):  %s', y_data)

        ax.step(mt_binning[:-1], y_bkg_init, where='post', c='purple', label=r'$B_{prefit}$')
        ax2.step(
            mt_binning[:-1], (y_bkg_init - y_data) / np.sqrt(y_data), where='post', c='purple',
            )

        mt_fine = np.linspace(mt_binning[0], mt_binning[-1], 100) # For fine plotting
        spl = make_interp_spline(mt_bin_centers, y_bkg, k=3)  # type of this is BSpline
        y_bkg_fine = spl(mt_fine)
        ax.plot(mt_fine, y_bkg_fine, label=r'$B_{fit}$', c='b')
        ax2.step(
            mt_binning[:-1], (y_bkg - y_data) / np.sqrt(y_data), where='post', c='b',
            )

        ax.step(
            mt_binning[:-1], y_sb, where='post', c='r',
            label=r'$B_{{fit}}+\mu_{{fit}}$S ($\mu_{{fit}}$={0:.1f})'.format(mu)
            )
        ax2.step(
            mt_binning[:-1], (y_sb - y_data) / np.sqrt(y_data), where='post', c='r',
            )

        ax.step(mt_binning[:-1], y_sig, where='post', label=r'S ($\mu$=1)', c='g')

    ax.legend(framealpha=0.0, fontsize=22)
    ax.set_ylabel('$N_{events}$')
    ax.set_xlabel(r'$m_{T}$ (GeV)')
    ax.set_yscale('log')
    ax2.set_ylabel('(pdf - data) / sqrt(data)', fontsize=18)
    if only_sig: ax2.set_ylabel('postfit / prefit', fontsize=18)
    ax.set_ylim(1., None)

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


@scripter
def cls():
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean

    for observed, asimov in zip(*organize_rootfiles(rootfiles)):
        obs, asimov = extract_scans([observed, asimov], correct_minimum=True)
        if clean:
            obs = clean_scan(obs)
            asimov = clean_scan(asimov)

        cls = get_cls(obs, asimov)
        limit = interpolate_95cl_limit(cls)

        logger.info(
            'Limit {}: 2sd={:.4f} 1sd={:.4f} exp={:.4f} 1su={:.4f} 2su={:.4f} obs={:.4f}'
            .format(
                get_mz(observed),
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

            ax.plot([], [], ' ', label=name_from_combine_rootfile(observed, True))
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
    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean

    obs_rootfiles, asimov_rootfiles = organize_rootfiles(rootfiles)

    points = []
    for obs_rootfile, asimov_rootfile in zip(obs_rootfiles, asimov_rootfiles):
        mz = get_mz(obs_rootfile)
        xsec = bsvj.get_xs(mz)
        assert mz == get_mz(asimov_rootfile)
        obs, asimov = extract_scans([obs_rootfile, asimov_rootfile])
        if clean:
            obs = clean_scan(obs)
            asimov = clean_scan(asimov)
        cls = get_cls(obs, asimov)
        limit = interpolate_95cl_limit(cls)
        points.append(bsvj.AttrDict(
            mz = mz,
            xsec = xsec,
            limit = limit,
            cls = cls
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
        ax.plot([],[],label=r'$m_{dark}$=10 GeV, $r_{inv}$=0.3',color='white')

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
        ax.set_xlabel(r'$m_{Z\prime}$ (GeV)')
        ax.set_ylim(0.1,50)
        ax.grid(True)
        #ax.set_ylabel(r'$\mu$')
        ax.set_ylabel(r'$\sigma \times BR$(pb)')
        ax.set_yscale('log')
        apply_ranges(ax)
        ax.legend(framealpha=0.0)


@scripter
def bkgfit():
    """
    Bkg fit plots
    """
    jsons = bsvj.get_jsons()
    pdftype = bsvj.pull_arg('pdftype', type=str, choices=bsvj.known_pdfs()).pdftype
    linscale = bsvj.pull_arg('--lin', action='store_true').lin
    scipyonly = bsvj.pull_arg('--scipyonly', action='store_true').scipyonly
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    gof_type = bsvj.pull_arg('--gof-type', type=str, default='chi2', choices=['chi2','rss']).gof_type

    input = bsvj.InputData(**jsons)

    mt = bsvj.get_mt(input.mt[0], input.mt[-1], input.n_bins, name='mt')
    bin_centers = .5*(input.mt_array[:-1]+input.mt_array[1:])
    bin_width = input.mt[1] - input.mt[0]
    bkg_hist = input.bkg['bkg']
    bkg_th1 = bkg_hist.th1('bkg')
    if input.data is not None:
        data_th1 = input.data['data'].th1('data')
    else:
       data_th1 = bkg_th1
    data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), data_th1, 1.)

    pdfs = bsvj.pdfs_factory(pdftype, mt, bkg_th1, name=pdftype)

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
    winner = bsvj.do_fisher_test(mt, data_datahist, pdfs, gof_type=gof_type)
    pdfs[winner].is_winner = True

    bkg_hist.vals = np.array(bkg_hist.vals)
    bkg_hist.shape = bkg_hist.vals / (bkg_hist.vals.sum()*bin_width)

    figure, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16))

    ax.plot([], [], ' ', label=f'{pdftype}, {input.metadata["selection"]}')
    ax.step(input.mt[:-1], bkg_hist.vals, where='post', label=r'BKG', c='b')
    ax2.plot([input.mt[0], input.mt[-1]], [0.,0.], c='gray')


    fine_mt_axis = np.linspace(input.mt[0], input.mt[-1], 100)
    for pdf in pdfs:
        par_vals = [p.getVal() for p in pdf.parameters]

        y_pdf = pdf.evaluate(bin_centers)
        if abs(1. - y_pdf.sum()) > 0.01: logger.error('PDF norm is off from 1.:', y_pdf.sum())

        if getattr(pdf, 'is_winner', False):
            logger.warning('par vals: %s', par_vals)
            logger.warning('y_pdf pre norm: %s (norm=%s)', y_pdf, y_pdf.sum())

        y_pdf *= bkg_hist.vals.sum()
        # y_pdf *= fit_norm(y_pdf, bkg_hist.vals) # Should be close to 1.0

        if getattr(pdf, 'is_winner', False):
            logger.warning('y_pdf post norm: %s (norm=%s)', y_pdf, y_pdf.sum())

        chi2_vf = bsvj.get_chi2_viaframe(mt, pdf, data_datahist)
        chi2 = chi2_vf['chi2']

        label = (
            '{}, $\\chi^{{2}}={:.5f}$: ['.format(pdf.n_pars, chi2)
            + ', '.join(['{:.2f}'.format(v) for v in par_vals]) + ']'
            )
        if getattr(pdf, 'is_winner', False): label += ' WINNER'

        y_pdf_fine = bsvj.eval_expression(pdf.expression, [fine_mt_axis] + par_vals)
        bin_scale = bin_width / (fine_mt_axis[1]-fine_mt_axis[0])
        y_pdf_fine = y_pdf_fine / y_pdf_fine.sum() * sum(bkg_hist.vals) * bin_scale
        line = ax.plot(fine_mt_axis, y_pdf_fine, label=label)[0]

        pulls = (y_pdf - bkg_hist.vals) / bkg_hist.errs
        ax2.scatter(bin_centers, pulls, color=line.get_color())


    ax.legend(fontsize=18, framealpha=0.0)
    ax.set_ylabel('$N_{events}$')

    ax2.set_ylabel(r'(pdf - bkg) / $\Delta$bkg')
    ax2.set_xlabel(r'$m_{T}$ (GeV)')
    if not linscale: ax.set_yscale('log')
    plt.savefig(outfile, bbox_inches='tight')
    if not(BATCH_MODE) and cmd_exists('imgcat'): os.system('imgcat ' + outfile)


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
