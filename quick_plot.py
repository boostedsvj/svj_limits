"""
Some scripts to quickly plot basic outcomes from combine scans
"""
from __future__ import print_function
import ROOT # type:ignore
from time import strftime
import argparse

# Add the directory of this file to the path so the boosted tools can be imported
import sys, os, os.path as osp, pprint, re, traceback, copy
from contextlib import contextmanager
sys.path.append(osp.dirname(osp.abspath(__file__)))
import boosted_fits as bsvj
logger = bsvj.setup_logger('quickplot')

import numpy as np
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
    """
    Checks if a command can be found on the system path.
    Not a very smart implementation but does the job usually.
    See https://stackoverflow.com/a/28909933/9209944 .
    """
    return any(os.access(os.path.join(path, executable), os.X_OK) for path in os.environ["PATH"].split(os.pathsep))

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


def name_from_combine_rootfile(rootfile, strip_obs_asimov=False):#, strip_pdfname=False):
    name = osp.basename(rootfile).rsplit('.',5)[0].replace('higgsCombine','')
    name = name.replace('_rinv0','')
    name = name.replace('dc','')
    print(name)
    #name = name.replace 
    if strip_obs_asimov:
        name = name.replace('Observed_','').replace('Asimov_','')
    '''if strip_pdfname:
        name = name.replace('main_','').replace('ua2','')'''
    name = name.replace('.MultiDimFit','')
    return name

def namespace_to_attrdict(args):
    return bsvj.AttrDict(**vars(args))

def get_mz(path):
    return int(re.search(r'mz(\d+)', osp.basename(path)).group(1))

def get_rinv(path):
    return int(re.search(r'rinv([\d\.]+)', osp.basename(path)).group(1))

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
 
    # mzs = {get_mz(rootfile) for rootfile in rootfiles}

    # for rootfile in rootfiles:
    #     mz = get_mz(rootfile)
    #     bdt = get_bdt_str(rootfile)

    # out = {}

    # for rootfile in rootfiles:
    #     mz = get_mz(rootfile)
    #     bdt = get_bdt_str(rootfile)

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

def extract_scans(rootfiles, correct_minimum=False):
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
                'iToy' : 'itoy',
                #'trackedParam_bsvj_bkgfitua2_npars3_p1' : 'p1',
                #'trackedParam_bsvj_bkgfitua2_npars3_p2' : 'p2',
                #'trackedParam_bsvj_bkgfitua2_npars3_p3' : 'p3'
                }
            # Add the tracked params
            listofbranches = limit.GetListOfBranches()
            for i_branch in range(listofbranches.GetEntries()):
                branch = listofbranches[i_branch].GetName()
                if branch.startswith('trackedParam_'):
                    keys[branch] = branch.replace('trackedParam_','')

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
                #print('*'*50, is_bestfit)
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

    return scans


# def extract_scans(rootfiles, correct_minimum=False):
#     if isinstance(rootfiles, str): rootfiles = [rootfiles]
#     scans = []

#     for rootfile in rootfiles:
#         with bsvj.open_root(rootfile) as tf:
#             limit = tf.Get('limit')

#             tracked_params = {}
#             listofbranches = limit.GetListOfBranches()
#             for i_branch in range(listofbranches.GetEntries()):
#                 branch = listofbranches[i_branch].GetName()
#                 if branch.startswith('trackedParam_'):
#                     tracked_params[branch.replace('trackedParam_', '')] = []

#             all_mus = []
#             all_deltanlls = []
#             all_quantiles = []
#             i_toys = []
#             for _ in limit:
#                 if limit.deltaNLL < 1.e5:
#                     all_mus.append(limit.r)
#                     all_deltanlls.append(limit.deltaNLL)
#                     all_quantiles.append(limit.quantileExpected)
#                     i_toys.append(limit.iToy)
#                     for key, vals in tracked_params.items():
#                         vals.append(getattr(limit, 'trackedParam_'+key))

#         all_mus = np.array(all_mus)
#         all_deltanlls = np.array(all_deltanlls)
#         all_quantiles = np.array(all_quantiles)
#         i_toys = np.array(i_toys)
#         for key, vals in list(tracked_params.items()):
#             tracked_params[key] = np.array(vals)

#         for i_toy in sorted(set(i_toys)):
#             scan = bsvj.AttrDict()
#             sel = (i_toys == i_toy)
#             mus = all_mus[sel]
#             deltanlls = all_deltanlls[sel]
#             quantiles = all_quantiles[sel]
#             order = np.argsort(mus)
#             mus = mus[order]
#             deltanlls = deltanlls[order]
#             quantiles[order]
#             if correct_minimum:
#                 minimum = np.min(deltanlls)
#                 logger.warning('Shifting curve by {0:.4f}'.format(minimum))
#                 deltanlls = deltanlls - minimum
#             scan.mus = mus
#             scan.deltanlls = deltanlls
#             scan.quantiles = quantiles
#             for key, vals in tracked_params.items(): scan[key] = vals
#             scans.append(scan)
#     return scans


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


def apply_ranges(ax):
    xmin = bsvj.read_arg('--xmin', type=float).xmin
    xmax = bsvj.read_arg('--xmax', type=float).xmax
    ymin = bsvj.read_arg('--ymin', type=float).ymin
    ymax = bsvj.read_arg('--ymax', type=float).ymax
    if xmax is not None: ax.set_xlim(right=xmax)
    if xmin is not None: ax.set_xlim(left=xmin)
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

        for rootfile in rootfiles:
            name = name_from_combine_rootfile(rootfile)
            for scan in extract_scans(rootfile, correctminimum):
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
        ax.set_ylim(-10, 25)
        apply_ranges(ax)
        ax.legend(framealpha=0., fontsize=25)


@scripter
def trackedparam():
    param = bsvj.pull_arg('param', type=str).param

    print('*'*15, param)

    rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    clean = bsvj.pull_arg('--clean', action='store_true').clean

    drwstr = '-'
    alpha = 1.
    if clean:
        drwstr += '-'
        alpha = .2

    with quick_ax(outfile=outfile) as ax:
        for rootfile in rootfiles:
            name = name_from_combine_rootfile(rootfile)
            #print('*'*30, name)
            for scan in extract_scans(rootfile):
                #print('*'*45, scan)
                ax.plot(scan.df['mu'], scan.df[param], drwstr, label=name, alpha=alpha, linewidth=2)
                if clean:
                    cln = clean_scan(scan)
                    ax.plot(cln.df['mu'], cln.df[param], label=name, linewidth=2)

        ax.set_xlabel('$\mu$')
        ax.set_ylabel(param)
        apply_ranges(ax)
        ax.legend(framealpha=0., fontsize=25)


@scripter
def mtdist():
    rootfile = bsvj.pull_arg('rootfile', type=str).rootfile
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='muscan.png').outfile

    from scipy.interpolate import make_interp_spline # type:ignore
    
    with bsvj.open_root(rootfile) as f:
        ws = bsvj.get_ws(f)

    mt = ws.var('mt')
    mt_binning = bsvj.binning_from_roorealvar(mt)
    mt_bin_centers = .5*(mt_binning[1:]+mt_binning[:-1])
    mt_bin_widths = mt_binning[1:] - mt_binning[:-1]

    y_bkg_init = bsvj.pdf_values(ws.pdf('shapeBkg_roomultipdf_bsvj'), mt_bin_centers)
    pdf_raw_norm_prefit = np.sum(y_bkg_init)
    logger.warning('y_bkg_prefit: %s', y_bkg_init)
    bkg_norm_init = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()
    y_bkg_init *= bkg_norm_init
    logger.warning('y_bkg_init after norm: %s', y_bkg_init)


    ws.loadSnapshot('MultiDimFit')

    mu = ws.var('r').getVal()
    print(r'$\mu$ is equal to = ', mu)

    data = ws.data('data_obs')
    _, y_data, _ = bsvj.roodataset_values(data)
    errs_data = np.sqrt(y_data)

    bkg = ws.pdf('shapeBkg_roomultipdf_bsvj')
    y_bkg = bsvj.pdf_values(bkg, mt_bin_centers)
    logger.warning('y_bkg_postfit: %s', y_bkg)
    pdf_raw_norm_postfit = np.sum(y_bkg)
    bkg_norm = ws.function('n_exp_final_binbsvj_proc_roomultipdf').getVal()
    y_bkg *= bkg_norm



    sig = ws.embeddedData('shapeSig_sig_bsvj')
    _, y_sig, _ = bsvj.roodataset_values(sig)


    bsvj.logger.info(
        'Found {} entries in data, {} entries in signal (should match with datacard!)'
        .format(y_data.sum(), y_sig.sum())
        )

    fig, (ax, ax2) = plt.subplots(
        2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16)
        )

    ax.plot([], [], ' ', label=name_from_combine_rootfile(rootfile))
    ax2.plot([mt_binning[0], mt_binning[-1]], [0,0], c='gray')

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

    y_sb = y_bkg + mu*y_sig
    ax.step(
        mt_binning[:-1], y_sb, where='post', c='r',
        label=r'$B_{{fit}}+\mu_{{fit}}$S ($\mu_{{fit}}$={0:.1f})'.format(mu)
        )
    ax2.step(
        mt_binning[:-1], (y_sb - y_data) / np.sqrt(y_data), where='post', c='r',
        )

    ax.step(mt_binning[:-1], y_sig, where='post', label=r'S ($\mu$=1)', c='g')

    ax.legend(framealpha=0.0, fontsize=25)
    ax.set_ylabel('$N_{events}$')
    ax.set_xlabel(r'$m_{T}$ (GeV)')
    ax.set_yscale('log')
    ax2.set_ylabel('(pdf - data) / sqrt(data)', fontsize=18)

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

    #     for i in range(1, len(mu)):
    #         if mu[i] == mu[i-1]: continue
    #         keep.append(i)
    #     if len(keep) < len(mu):
    #         logger.warning('Removing {} duplicate mu values'.format(len(mu) - len(keep)))
    #         mu = mu[keep]
    #         dnll = dnll[keep]
    #         quantile = quantile[keep]


    np.testing.assert_array_equal(obs.df['mu'], asimov.df['mu'])


    # I do not understand why not simply q_obs = 2.*dnll_obs?
    # Or is this just to offset a potentially bad best fit?
    # If so, why not just shift the whole dnll array so its minimum is at 0...
    dnll_obs = obs.df['dnll']
    q_obs = []
    for i, mu in enumerate(obs.df['mu']):
        if mu < obs.bestfit.df['mu']:
            dnll_obs_min = np.min(dnll_obs[:i+1])  # Why?
            dnll_obs_constrained = dnll_obs[i] - dnll_obs_min
        else:
            dnll_obs_constrained = dnll_obs[i]
        q_obs.append(2.*max(dnll_obs_constrained, 0.))
    q_obs = np.array(q_obs)
    assert q_obs.shape == (obs.n,)

    q_A = 2. * asimov.df['dnll']
    q_A[q_A < 0.] = 0.  # Set negative values to 0

    # Also this formula I don't fully understand
    s_exp = { q : (1.-norm.cdf(np.sqrt(q_A) - norm.ppf(q))) / q for q in quantiles}

    assert np.all(  ((q_obs >= 0.) & (q_obs <= q_A)) | (q_obs > q_A)  )

    # This is just black magic
    sb = np.where(
        q_obs <= q_A,
        1. - norm.cdf( np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs+q_A) , np.sqrt(q_obs)) )
        )
    b = np.where(
        q_obs <= q_A,
        norm.cdf( np.sqrt(q_A)-np.sqrt(q_obs) ),
        1. - norm.cdf( safe_divide(.5*(q_obs-q_A) , np.sqrt(q_obs)) )
        )
    s = sb / b
    return bsvj.AttrDict(s=s, b=b, sb=sb, q_obs=q_obs, q_A=q_A, obs=obs, asimov=asimov, s_exp=s_exp)


def interpolate_95cl_limit(cls):
    mu = cls.obs.df['mu']
    def interpolate(cl, thing):
        # print('Interpolating')
        # select = ((cl < .20) & (mu>0))
        select = ((cl < .50) & (cl > .001) & (mu>0))
        if select.sum() == 0:
            logger.error('0.01<cl<0.20 & mu>0 yields NO scan points; can\'t interpolate %s', thing)
            return None

        # print('  {} values left'.format(select.sum()))
        order = np.argsort(cl[select])
        # print('  {:14s}  {:14s}'.format('cl', 'mu'))
        # for c, m in zip(cl[select][order], mu[select][order]):
        #     print('  {:+14.7f}  {:+14.7f}'.format(c, m))
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
        # print('Interpolation result: cl=0.05, mu={}'.format(res))
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
        assert mz == get_mz(asimov_rootfile)
        obs, asimov = extract_scans([obs_rootfile, asimov_rootfile])
        if clean:
            obs = clean_scan(obs)
            asimov = clean_scan(asimov)
        cls = get_cls(obs, asimov)
        limit = interpolate_95cl_limit(cls)
        points.append(bsvj.AttrDict(
            mz = mz,
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

        ax.fill_between(
            [p.mz for p in points if p.limit.twosigma_success],
            [p.limit.twosigma_down for p in points if p.limit.twosigma_success],
            [p.limit.twosigma_up for p in points if p.limit.twosigma_success],
            color=cms_yellow
            )
        ax.fill_between(
            [p.mz for p in points if p.limit.onesigma_success],
            [p.limit.onesigma_down for p in points if p.limit.onesigma_success],
            [p.limit.onesigma_up for p in points if p.limit.onesigma_success],
            color=cms_green
            )

        ax.plot(
            [p.mz for p in points if p.limit.expected is not None],
            [p.limit.expected for p in points if p.limit.expected is not None],
            c='black', linestyle='--', label='Exp'
            )
        ax.plot(
            [p.mz for p in points if p.limit.observed is not None],
            [p.limit.observed for p in points if p.limit.observed is not None],
            c='black', linestyle='-', label='Obs'
            )

        ax.set_xlabel(r'$m_{Z\prime}$ (GeV)')
        ax.set_ylabel(r'$\mu$')
        apply_ranges(ax)
        ax.legend(framealpha=0.0)






@scripter
def bkgfit():
    """
    Bkg fit plots
    """
    jsonfile = bsvj.pull_arg('jsonfile', type=str).jsonfile
    bdtcut = bsvj.pull_arg('bdtcut', type=float).bdtcut
    pdftype = bsvj.pull_arg('pdftype', type=str, choices=['main', 'alt', 'ua2']).pdftype
    logscale = bsvj.pull_arg('--log', action='store_true').log
    trigeff = bsvj.pull_arg('--trigeff', type=int, default=None, choices=[2016, 2017, 2018]).trigeff
    fitmethod = bsvj.pull_arg('--fitmethod', type=str, choices=['scipy', 'auto'], default='scipy').fitmethod
    outfile = bsvj.read_arg('-o', '--outfile', type=str, default='test.png').outfile
    mtrange = bsvj.pull_arg('--range', type=float, nargs=2).range

    input = bsvj.InputData(jsonfile)
    if mtrange is None: mtrange = [180., 650.]
    #mtrange = [200., 600.]
    input = input.cut_mt(mtrange[0], mtrange[1])

    bdt_str = '{:.1f}'.format(bdtcut).replace('.', 'p')
    mt = bsvj.get_mt(input.mt[0], input.mt[-1], input.n_bins, name='mt')
    bin_centers = .5*(input.mt_array[:-1]+input.mt_array[1:])
    bin_width = input.mt[1] - input.mt[0]
    bkg_hist = input.bkg_hist(bdtcut)

    if trigeff:
        import requests
        parameters = np.array(requests.get('https://raw.githubusercontent.com/boostedsvj/triggerstudy/main/bkg/bkg_trigeff_fit_2018.txt').json())
        poly = np.poly1d(parameters)
        f_trig_eff = lambda x: np.where(x<1000., 1./(1.+np.exp(-poly(x))), 1.)
        binning = np.array(bkg_hist.binning)
        mt_bin_centers = .5*(binning[:-1]+binning[1:])
        bkg_hist.vals *= f_trig_eff(mt_bin_centers)
        logger.info('Adjusted bkg histogram for trigger eff')

    bkg_th1 = input.bkg_th1('bkg', bdtcut)

    data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_th1, 1.)

    pdfs = bsvj.pdfs_factory(pdftype, mt, bkg_th1, name=pdftype, trigeff=None)

    import scipy
    from scipy import optimize
    for pdf in pdfs:
        ranges=[]
        par_name=[]
        for par in pdf.parameters:
            ranges.append((par.getMin(), par.getMax()))
            par_name.append(par.GetName())
        ranges=tuple(ranges)
        par_name=tuple(par_name)
        par_name = '('+', '.join(map(lambda i: format(i, '14s'), par_name))+')'
    
        if fitmethod == 'auto':
            pdf.res = bsvj.fit(pdf)
        elif fitmethod == 'scipy':
            print('par_name: ', par_name, ' type: ', type(par_name))
            #pdf.res_bf = scipy.optimize.brute(pdf.expression, ranges=ranges, Ns=20, args=par_name, full_output=True, finish=optimize.fmin)
            pdf.res = bsvj.fit_scipy_robust(pdf.expression, pdf.th1, cache=None)
            # Fill in the fitted parameters
            for p, val in zip(pdf.parameters, pdf.res.x):
                p.setVal(val)


    # bsvj.plot_fits(pdfs, [p.res for p in pdfs], data_datahist, 'qp_' + pdftype + '.pdf')

    if fitmethod == 'auto':
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
    winner = bsvj.do_fisher_test(mt, data_datahist, pdfs)
    pdfs[winner].is_winner = True

    bkg_hist.vals = np.array(bkg_hist.vals)
    bkg_hist.shape = bkg_hist.vals / (bkg_hist.vals.sum()*bin_width)
    
    figure, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12,16))

    #ax.plot([], [], ' ', label='{}, score$\\geq${:.1f}'.format(pdftype, bdtcut))
    ax.plot([], [], ' ', label='{}, cutBased'.format(pdftype))
    ax.step(input.mt[:-1], bkg_hist.vals, where='post', label=r'BKG', c='b')
    ax2.plot([input.mt[0], input.mt[-1]], [0.,0.], c='gray')

    def fit_norm(y, y_data):
        from scipy.optimize import minimize # type:ignore
        def fom(norm):
            y_norm = y * norm
            return ((y_norm-y_data)**2 / y_norm).sum()
        res = minimize(fom, 1.)
        return res.x

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

        chi2 = ((y_pdf-bkg_hist.vals)**2 / y_pdf).sum()
        # chi2 /= (len(bin_centers) - pdf.npars)

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
    if logscale: ax.set_yscale('log')
    plt.savefig(outfile, bbox_inches='tight')
    if not(BATCH_MODE) and cmd_exists('imgcat'): os.system('imgcat ' + outfile)


if __name__ == '__main__':
    batch_mode(bsvj.pull_arg('-b', '--batch', action='store_true').batch)
    debug(bsvj.pull_arg('-d', '--debug', action='store_true').debug)
    fontsize = bsvj.read_arg('--fontsize', type=int, default=18).fontsize
    set_mpl_fontsize(legend=fontsize)
    scripter.run()
