"""
Building blocks to create the boosted SVJ analysis datacards
"""

import uuid
from contextlib import contextmanager
from array import array
from math import sqrt
import numpy as np
import itertools, re, logging, os, os.path as osp, copy, subprocess, json
from collections import OrderedDict
from time import strftime

import ROOT # type:ignore
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetPadColor(0)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")


DEFAULT_LOGGING_LEVEL = logging.INFO
def setup_logger(name='boosted'):
    if name in logging.Logger.manager.loggerDict:
        logger = logging.getLogger(name)
        logger.info('Logger %s is already defined', name)
    else:
        fmt = logging.Formatter(
            fmt = '\033[33m%(levelname)s:%(asctime)s:%(module)s:%(lineno)s\033[0m %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
            )
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logger = logging.getLogger(name)
        logger.setLevel(DEFAULT_LOGGING_LEVEL)
        logger.addHandler(handler)
    return logger
logger = setup_logger()
subprocess_logger = setup_logger('subp')
subprocess_logger.handlers[0].formatter._fmt = '\033[34m[%(asctime)s]\033[0m %(message)s'


def debug(flag=True):
    """Sets the logger level to debug (for True) or warning (for False)"""
    logger.setLevel(logging.DEBUG if flag else DEFAULT_LOGGING_LEVEL)


def mpl_fontsizes(small=14, medium=18, large=24):
    import matplotlib.pyplot as plt # type:ignore
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('figure', titlesize=large)   # fontsize of the figure title


def uid():
    return str(uuid.uuid4())


class AttrDict(dict):
    """
    Like a dict, but with access via attributes
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


# _______________________________________________________________________
# JSON Input interface

def hist_to_th1(name, hist):
    """
    Takes a dict-like histogram (keys: binning, vals, errs) and
    returns a ROOT.TH1F
    """
    n_bins = len(hist.binning)-1
    th1 = ROOT.TH1F(name, name, n_bins, array('f', hist.binning))
    ROOT.SetOwnership(th1, False)
    assert len(hist.vals) == n_bins
    assert len(hist.errs) == n_bins
    for i in range(n_bins):
        th1.SetBinContent(i+1, hist.vals[i])
        th1.SetBinError(i+1, hist.errs[i])
    return th1


def hist_cut_left(hist, i_bin_min):
    """
    Returns a new hist with i_bin_min as bin 0, and any preceeding bin
    thrown out
    """
    if i_bin_min == 0: return hist
    histcopy = dict(hist)
    for key in ['binning', 'vals', 'errs']:
        histcopy[key] = hist[key][i_bin_min:]
    return AttrDict(histcopy)


class Histogram:
    def __init__(self, d):
        self.vals = d['vals']
        self.errs = d['errs']
        self.binning = d['binning']
        self.metadata = d['metadata']



def build_histograms_in_dict_tree(d, parent=None, key=None):
    """
    Traverses a dict-of-dicts, and converts everything that looks like a 
    histogram to a Histogram object.
    """
    n_histograms = 0
    if not isinstance(d, dict): return 0
    if 'binning' in d: # It's a histogram / deprecate to 'type' : 'histogram' later
        parent[key] = Histogram(d)
        return 1
    else:
        for key, value in d.items():
            n_histograms += build_histograms_in_dict_tree(value, parent=d, key=key)
    return n_histograms


def iter_histograms(d):
    """
    Traverses a dict-of-dicts, and yields all the Histogram instances.
    """
    if isinstance(d, Histogram):
        yield d
    elif isinstance(d, dict):
        for v in d.values():
            for _ in iter_histograms(v):
                yield _

def ls_inputdata(d, depth=0, key='<root>'):
    """
    Prints a dict of dicts recursively
    """
    if isinstance(d, Histogram):
        print('  '*depth + key + ' (histogram)')
    elif isinstance(d, dict):
        print('  '*depth + key)
        for k, v in sorted(d.items()):
            ls_inputdata(v, depth+1, k)


class InputData(object):
    """
    Interface class for a JSON file that contains histograms
    """
    def __init__(self, jsonfile):
        self.jsonfile = jsonfile
        with open(jsonfile, 'r') as f:
            self.d = json.load(f)
            n_hists = build_histograms_in_dict_tree(self.d)
            logger.info('Read %s histograms from %s', n_hists, jsonfile)

    def copy(self):
        return copy.deepcopy(self)

    def ls(self):
        ls_inputdata(self.d, key=self.jsonfile)

    @property
    def version(self):
        return self.d.get('version', 1)

    def cut_mt(self, min_mt=None, max_mt=None):
        """
        Returns a copy of the InputData object with a limited mt range.
        """
        copy = self.copy()
        i_bin_min = 0
        i_bin_max = len(copy.mt)-1
        if min_mt is not None:
            i_bin_min = np.argmax(copy.mt_array >= min_mt)
        if max_mt is not None:
            i_bin_max = np.argmax(copy.mt_array >= max_mt)
        copy.mt = copy.mt[i_bin_min:i_bin_max]
        logger.info(
            'Cutting histograms {}<mt<{}, i_bin_min={}, i_bin_max={}'
            .format(min_mt, max_mt, i_bin_min, i_bin_max)
            )
        for h in iter_histograms(copy.d):
            h.binning = copy.mt
            h.vals = h.vals[i_bin_min:i_bin_max-1]
            h.errs = h.errs[i_bin_min:i_bin_max-1]
        return copy

    @property
    def mt(self):
        return self.d['mt']

    @property
    def mt_centers(self):
        return .5*(self.mt_array[:-1] + self.mt_array[1:])

    @property
    def mt_binwidth(self, i=0):
        return self.mt[i+1] - self.mt[i]

    @mt.setter
    def mt(self, val):
        self.d['mt'] = val

    @property
    def mt_array(self):
        return np.array(self.mt)

    @property
    def n_bins(self):
        return len(self.mt)-1

    def bkg_hist(self, bdt):
        bdt = float(bdt)
        if self.version == 1:
            # Backward compatibility with the first json implementation
            return self.d['histograms']['{:.1f}/bkg'.format(bdt)]
        else:
            return self.d['histograms']['{:.3f}'.format(bdt)]['bkg']

    def bkg_th1(self, name, bdt):
        return hist_to_th1(name, self.bkg_hist(bdt))

    def sig_hist(self, bdt, mz, rinv=.3, mdark=10):
        bdt = float(bdt)
        if self.version == 1:
            # Backward compatibility with the first json implementation
            return self.d['histograms']['{:.1f}/mz{:.0f}'.format(bdt, mz)]
        else:
            for h in self.d['histograms']['{:.3f}'.format(bdt)].values():
                if 'mz' not in h.metadata: continue
                if (
                    h.metadata['mz'] == mz
                    and h.metadata['rinv'] == rinv
                    and h.metadata['mdark'] == mdark
                    ):
                    return h
            raise Exception(
                'Could not find a histogram for mz={}, rinv={}, mdark={}'
                .format(mz, rinv, mdark)
                )

    def sig_th1(self, name, bdt, mz, rinv=.3, mdark=10):
        return hist_to_th1(name, self.sig_hist(bdt, mz, rinv, mdark))


# _______________________________________________________________________
# Model building code: Bkg fits, fisher testing, etc.


def dump_fits_to_file(filename, results):
    logger.info('Dumping fit results to ' + filename)
    dirname = osp.dirname(osp.abspath(filename))
    if not osp.isdir(dirname): os.makedirs(dirname)
    with open_root(filename, 'RECREATE') as tf:
        for result in results: result.Write()


def dump_ws_to_file(filename, ws):
    logger.info('Dumping ws {} to {}'.format(ws.GetName(), filename))
    dirname = osp.dirname(osp.abspath(filename))
    if not osp.isdir(dirname): os.makedirs(dirname)
    wstatus = ws.writeToFile(filename, True)
    return wstatus


def eval_expression(expression, pars):
    """
    Evaluates a ROOT TFormula expression in python.
    Only a limited amount of keywords are implemented (pow, log, sqrt, exp).
    """
    # Load keywords in local scope
    from numpy import log, sqrt, exp
    def pow(base, exponent):
        return base ** exponent
    # Python variables can't start with '@'; replace with some keyword
    expression = expression.replace('@', 'PARAMETER')
    # Plug parameters in local scope
    par_dict = {'PARAMETER'+str(i) : p for i, p in enumerate(pars)}
    locals().update(par_dict)
    # logger.warning('Evaluating expr:\n%s\nwith parameters:\n%s', expression, par_dict)
    try:
        return eval(expression)
    except NameError:
        logger.error(
            'Missing variables for expression:\n{0}\nAvailable parameters: {1}'
            .format(expression, list(par_dict.keys()))
            )
        raise


def eval_pdf_python(pdf, parameters, mt_array=None):
    if mt_array is None:
        mt = pdf.parameters[0]
        binning = mt.getBinning()
        mt_array = np.array([ binning.binCenter(i) for i in range(binning.numBins()) ])    
    parameters = list(copy.copy(parameters))
    parameters.insert(0, mt_array)
    return eval_expression(pdf.expression, parameters)


def count_parameters(expr):
    """Returns the number of parameters in an expression (i.e. highest @\d"""
    return max(map(int, re.findall(r'@(\d+)', expr))) + 1


def add_normalization(expr):
    """
    Takes an expression string, and basically adds "@NORM*(...)" around it.
    """
    return '@{0}*('.format(count_parameters(expr)) + expr + ')'


def build_rss(expr, th1):
    """
    Builds a residual-sum-of-squares function between a pdf (expression)
    and a histogram.
    """
    # binning, counts = th1_binning_and_values(h)
    hist = th1_to_hist(th1)
    # counts /= counts.sum() # normalize to 1
    bin_centers = [.5*(l+r) for l, r in zip(hist.binning[:-1], hist.binning[1:])]
    mtarray = np.array(bin_centers)
    def rss(parameters):
        # Insert mT array as first parameter
        parameters = list(copy.copy(parameters))
        parameters.insert(0, mtarray)
        y_pdf = eval_expression(expr, parameters)
        # Normalize pdf to counts too before evaluating, so as the compare only shape
        y_pdf = (y_pdf/y_pdf.sum()) * hist.vals.sum()
        return np.sqrt(np.sum((hist.vals-y_pdf)**2))
    return rss


def build_chi2(expr, h):
    """
    Builds a chi2 function between a pdf (expression) and a histogram.
    """
    hist = th1_to_hist(h)
    # Use the bin centers as the mT array
    mt_array = np.array(.5*(hist.binning[:-1]+hist.binning[1:]))
    def chi2(parameters):
        # Construct the parameters of the expression:
        # [ @0 (mt), @1, ... @N (pdf parameters) ]
        parameters = list(copy.copy(parameters))
        parameters.insert(0, mt_array)
        y_pdf = eval_expression(expr, parameters)
        # Normalize pdf to counts too before evaluating, so as the compare only shape
        y_pdf = (y_pdf/y_pdf.sum()) * hist.vals.sum()
        return np.sum((hist.vals-y_pdf)**2 / y_pdf)
    return chi2


def make_fit_hash(expression, th1, init_vals=None, tag=None, **minimize_kwargs):
    """
    Constructs a hash from all the input data of a fit:
    - The expression (as a string)
    - The histogram (binning and values, not errors)
    - The initial values set for the fit
    - The scipy minimizer arguments
    - Any user provided tag
    """
    import hashlib
    m = hashlib.sha256()
    def add_floats_to_hash(floats):
        for number in floats:
            s = '{:.5f}'.format(number)
            m.update(s)
    m.update(expression)
    hist = th1_to_hist(th1)
    add_floats_to_hash(hist.binning)
    add_floats_to_hash(hist.vals)
    if init_vals is not None: add_floats_to_hash(init_vals)
    if 'tol' in minimize_kwargs: m.update('{:.3f}'.format(minimize_kwargs['tol']))
    if 'method' in minimize_kwargs: m.update(minimize_kwargs['method'])
    if tag: m.update(tag)
    return m.hexdigest()


def fit_roofit(pdf, data_hist=None, init_vals=None, init_ranges=None):
    """
    Main bkg fit entry point for fitting pdf to bkg th1 with RooFit
    """
    # Preparation
    if data_hist is None: data_hist = pdf.th1
    if isinstance(data_hist, ROOT.TH1): data_hist = th1_to_datahist(data_hist, pdf.mt)

    logger.info('Fitting pdf {0} to data_hist {1} with RooFit'.format(pdf, data_hist))

    if init_vals is not None:
        if len(init_vals) != len(pdf.parameters):
            raise Exception('Expected {} values; got {}'.format(len(pdf.parameters)-1, len(init_vals)))
        for par, value in zip(pdf.parameters, init_vals):
            old_range = max(abs(par.getMin()), abs(par.getMax()))
            if abs(value) > old_range:
                new_range = 1.1*abs(value)
                logger.info(
                    'Increasing range for {} ({}) from ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})'
                    .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax(), -new_range, new_range)
                    )
                par.setRange(-new_range, new_range)
            elif abs(value) / old_range < .1:
                new_range = 1.1*abs(value)
                logger.info(
                    'Decreasing range for {} ({}) from ({:.3f}, {:.3f}) to ({:.3f}, {:.3f})'
                    .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax(), -new_range, new_range)
                    )
                par.setRange(-new_range, new_range)
            # par.setRange(value-0.01*abs(value), value+0.01*abs(value))
            par.setVal(value)
            logger.info(
                'Setting {0} ({1}) value to {2}, range is {3} to {4}'
                .format(par.GetName(), par.GetTitle(), value, par.getMin(), par.getMax())
                )
    if init_ranges is not None:
        if len(init_ranges) != len(pdf.parameters):
            raise Exception('Expected {} values; got {}'.format(len(pdf.parameters), len(init_ranges)))
        for par, (left, right) in zip(pdf.parameters, init_ranges):
            par.setRange(left, right)
            logger.info(
                'Setting {0} ({1}) range to {2} to {3}'
                .format(par.GetName(), par.GetTitle(), par.getMin(), par.getMax())
                )

    try:
        res = pdf.pdf.fitTo(
            data_hist,
            ROOT.RooFit.Extended(False),
            ROOT.RooFit.Save(1),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Minimizer("Minuit2"),
            ROOT.RooFit.PrintLevel(2 if logger.level <= logging.DEBUG else -1),
            ROOT.RooFit.Range('Full'),
            ROOT.RooFit.PrintEvalErrors(-1)
            )
    except:
        logger.error('Problem fitting pdf {}'.format(pdf.pdf.GetName()))
        raise

    if logger.level <= logging.INFO: res.Print()
    return res


def single_fit_scipy(expression, histogram, init_vals=None, cache=None, **minimize_args):
    """
    Fits a RooFit-style expression (as a string) to a TH1 histogram.

    If cache is a FitCache object, the fit result is stored in the cache.
    """
    fit_hash = make_fit_hash(expression, histogram, init_vals=init_vals, **minimize_args)
    if cache and cache.get(fit_hash):
        logger.info('Returning cached fit')
        return cache.get(fit_hash) # Second call is cheap
    # Do the fit
    n_fit_pars = count_parameters(expression) - 1 # -1 because par 0 is mT
    logger.info('Fitting {0} with {1} parameters'.format(expression, n_fit_pars))
    from scipy.optimize import minimize # type:ignore
    chi2 = build_chi2(expression, histogram)
    if init_vals is None: init_vals = np.ones(n_fit_pars)
    res = minimize(chi2, init_vals, **minimize_args)
    # Save some extra information in the result
    res.x_init = np.array(init_vals)
    res.expression = expression
    res.hash = hash
    # Set approximate uncertainties; see https://stackoverflow.com/a/53489234
    # Assume ftol ~ function value
    # try:
    #     res.dx = np.sqrt(res.fun * np.diagonal(res.hess_inv))
    # except:
    #     logger.error('Failed to set uncertainties; using found function values as proxies')
    #     res.dx = res.x.copy()
    if cache:
        logger.info('Writing fit to cache')
        cache.write(fit_hash, res)
    return res


def fit_scipy_robust(expression, histogram, cache='auto'):
    """
    Main entry point for fitting an expression to a histogram with Scipy
    """
    logger.info('Robust scipy fit of expression %s to %s', expression, histogram)
    fit_hash = make_fit_hash(expression, histogram, tag='robust')

    if cache is 'auto':
        from fit_cache import FitCache # type: ignore
        cache = FitCache()

    if cache and cache.get(fit_hash):
        res = cache.get(fit_hash) # Second call is cheap
        logger.info('Returning cached fit:\n%s', res)
        return res
    
    # Attempt 1: Fit with loose tolerance BFGS, then strict tolerance Nelder-Mead
    res = single_fit_scipy(
        expression, histogram,
        tol=1e-3, method='BFGS',
        cache=cache
        )
    # Refit with output from first fit
    res = single_fit_scipy(
        expression, histogram,
        init_vals=res.x,
        tol=1e-6, method='Nelder-Mead',
        cache=cache
        )

    if res.success:
        # Fit successful, save in the cache and return
        if cache: cache.write(fit_hash, res)
        logger.info('Converged with simple fitting strategy, result:\n%s', res)
        return res

    # The simple fitting scheme failed; Brute force with many different
    # initial values
    npars = count_parameters(expression)-1 # The mT parameter is not a fit parameter
    init_val_variations = [-1., 1.] # All the possible init values a single fit parameter can have
    init_vals = np.array(list(itertools.product(*[init_val_variations for i in range(npars)])))
    logger.info(
        'Fit did not converge with single try; brute forcing it with '
        '%s different variations of initial values with both BFGS and Nelder-Mead.',
        len(init_vals)
        )
    results = []
    for method in ['BFGS', 'Nelder-Mead']:
        for init_val_variation in init_vals:
            result = single_fit_scipy(
                expression, histogram,
                init_vals=init_val_variation,
                tol=1e-3, method=method
                )
            # Check if fit fn val is not NaN or +/- inf
            if not(np.isnan(result.fun) or np.isposinf(result.fun) or np.isneginf(result.fun)):
                results.append(result)
    if len(results) == 0: raise Exception('Not a single fit of the brute force converged!')
    i_min = np.argmin([r.fun for r in results])        
    res = results[i_min]
    logger.info('Best scipy fit from brute force:\n%s', res)
    if cache: cache.write(fit_hash, res)
    return res


def fit(pdf, th1=None, cache='auto'):
    """
    Main bkg fit entry point for
    - first fitting pdf expression to bkg th1 with scipy
    - then using those initial values in RooFit
    """
    if th1 is None: th1 = getattr(pdf, 'th1', None)
    res_scipy = fit_scipy_robust(pdf.expression, th1, cache=cache)
    res_roofit_wscipy = fit_roofit(pdf, th1, init_vals=res_scipy.x)
    return res_roofit_wscipy



def get_mt(mt_min, mt_max, n_bins, name=None):
    """
    Sensible defaults for the mt axis
    """
    if name is None: name = uid()
    mt = ROOT.RooRealVar(name, 'm_{T}', mt_min, mt_max, 'GeV')
    mt.setBins(n_bins)
    # Manually add the boundaries to it as python attributes for easy access
    mt.mt_min = mt_min
    mt.mt_max = mt_max
    return mt


def get_mt_from_th1(histogram, name=None):
    """
    Returns mT from the x axis of a TH1 histogram.
    Min and max are simply the left/right boundary of the first/last bin,
    and bin width is copied.
    """
    mt = get_mt(
        histogram.GetBinLowEdge(1),
        histogram.GetBinLowEdge(histogram.GetNbinsX()+1),
        histogram.GetNbinsX(),
        name = uid() if name is None else name
        )
    object_keeper.add(mt)
    return mt


def rebuild_rpsbp(pdf):
    name = uid()
    def remake_parameter(parameter):
        variable = ROOT.RooRealVar(
            name + '_' + parameter.GetName(), parameter.GetTitle(),
            1., parameter.getMin(), parameter.getMax()
            )
        object_keeper.add(variable)
        return variable
    return build_rpsbp(
        name, pdf.expression, pdf.mt,
        [remake_parameter(p) for p in pdf.parameters], pdf.th1
        )


def trigeff_expression(year=2018, max_fit_range=800.):
    """
    Returns a TFormula-style expression that represents the trigger
    efficiency as a function of MT.

    The formula contains a switch based on `max_fit_range`: evaluating
    beyond `max_fit_range` will return 1. only.
    (This is needed because the fit is unstable under extrapolation)
    """
    import requests
    parameters = np.array(requests.get(
        'https://raw.githubusercontent.com/boostedsvj/triggerstudy/main/bkg/bkg_trigeff_fit_{}.txt'
        .format(year)).json())
    expr = sigmoid(poly1d(parameters))
    return '({0})*(@0<{1}) + (@0>={1})'.format(expr, max_fit_range)

def poly1d(parameters, mt_par='@0'):
    degree = len(parameters)
    return '+'.join([ '{}*pow({},{})'.format(p, mt_par, degree-i-1) for i, p in enumerate(parameters)])

def sigmoid(expr):
    return '1./(1.+exp(-({})))'.format(expr)


def pdf_expression(pdf_type, npars, mt_scale='1000'):
    # Function from Theorists, combo testing, sequence E, 1, 11, 12, 22
    # model NM has N params on 1-x and M params on x. exponents are (p_i + p_{i+1} * log(x))
    if pdf_type == 'main':
        if npars == 2:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2))'
        elif npars == 3:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})))'
        elif npars == 4:
            expression = 'pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})+@4*pow(log(@0/{0}),2)))'
            # Alternatives to 22:
            # 13: pow(1 - @0/{0}, @1+@2*log(@0/{0})) * pow(@0/{0}, -(@3+@4*log(@0/{0})))
        elif npars == 5:
            expression = 'pow(1 - @0/{0}, @1+@2*log(@0/{0})+@3*pow(log(@0/{0}),2)) * pow(@0/{0}, -(@4+@5*log(@0/{0})))'
            # Alternatives to 32:
            # 14: pow(1 - @0/{0}, @1) * pow(@0/{0}, -(@2+@3*log(@0/{0})+@4*pow(log(@0/{0}),2)+@5*pow(log(@0/{0}),3)))
            # 41: pow(1 - @0/{0}, @1+@2*log(@0/{0})+@3*pow(log(@0/{0}),2)+@4*pow(log(@0/{0}),3)) * pow(@0/{0}, -@5)
        else:
            raise Exception('Unavailable npars for main: {0}'.format(npars))
    elif pdf_type == 'alt':
        if npars == 1:
            expression = 'exp(@1*(@0/{0}))'
        elif npars == 2:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2)'
        elif npars == 3:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2*(1+@3*log(@0/{0})))'
        elif npars == 4:
            expression = 'exp(@1*(@0/{0})) * pow(@0/{0},@2*(1+@3*log(@0/{0})*(1+@4*log(@0/{0}))))'
        else:
            raise Exception('Unavailable npars for alt: {0}'.format(npars))
    else:
        raise Exception('Unknown pdf type {0}'.format(pdf_type))
    return expression.format(mt_scale)


def pdf_parameters(pdf_type, npars, prefix=None):
    if prefix is None: prefix = uid()
    if pdf_type == 'main':
        if npars == 2:
            parameters = [
                ROOT.RooRealVar(prefix + "_p1", "p1", 1., -45., 45.),
                ROOT.RooRealVar(prefix + "_p2", "p2", 1., -10., 10.)
                ]
        elif npars == 3:
            parameters = [
                ROOT.RooRealVar(prefix + "_p1", "p1", 1., -45., 45.),
                ROOT.RooRealVar(prefix + "_p2", "p2", 1., -10., 10.),
                ROOT.RooRealVar(prefix + "_p3", "p3", 1., -15, 15),
                ]
        elif npars == 4:
            parameters = [
                ROOT.RooRealVar(prefix + "_p1", "p1", 1., -95., 95.),
                ROOT.RooRealVar(prefix + "_p2", "p2", 1., -25., 20.),
                ROOT.RooRealVar(prefix + "_p3", "p3", 1., -2., 2.),
                ROOT.RooRealVar(prefix + "_p4", "p4", 1., -2., 2.),
                ]
        elif npars == 5:
            parameters = [
                ROOT.RooRealVar(prefix + "_p1", "p1", 1., -15., 15.),
                ROOT.RooRealVar(prefix + "_p2", "p2", 1., -95., 95.),
                ROOT.RooRealVar(prefix + "_p3", "p3", 1., -25., 25.),
                ROOT.RooRealVar(prefix + "_p4", "p4", 1., -5., 5.),
                ROOT.RooRealVar(prefix + "_p5", "p5", 1., -1.5, 1.5),
                ]
    elif pdf_type == 'alt':
        par_lo = -50.
        par_up = 50.
        parameters = [
            ROOT.RooRealVar(prefix + '_p{0}'.format(i+1), '', 1., par_lo, par_up) \
            for i in range(npars)
            ]
    object_keeper.add_multiple(parameters)    
    return parameters


class PDF(object):
    """
    Container object for a RooParametricShapeBinPdf, with more info
    """
    def __init__(self):
        pass

    def __repr__(self):
        return (
            '<RooParametricShapeBinPdf "{}"'
            '\n  pdf_type   = {}'
            '\n  n_pars     = {}'
            '\n  expression = "{}"'
            '\n  pdf        = "{}"'
            '\n  mt         = "{}"'
            '\n  th1        = "{}"'
            '\n  rgp        = "{}"'
            '\n  parameters = \n    {}'
            '\n  >'
            .format(
                self.name, self.pdf_type, self.n_pars, self.expression,
                self.pdf.GetName(), self.mt.GetName(), self.th1.GetName(), self.rgp.GetName(),
                '\n    '.join(['"'+p.GetName()+'"' for p in self.parameters])
                )
            )

    def evaluate(self, x_vals, varname=None):
        """
        Equivalent to the y_values of pdf.createHistogram('...', mt).
        Result is normalized to 1.
        """
        variable = self.mt if varname is None else self.pdf.getVariables()[varname]
        y = []
        for x in x_vals:
            variable.setVal(x)
            y.append(self.pdf.getVal())
        y = np.array(y)
        return y / (y.sum() if y.sum()!=0. else 1.)


def pdf_factory(pdf_type, n_pars, mt, bkg_th1, name=None, mt_scale='1000', trigeff=None):
    """
    Main factory entry point to generate a single RooParametricShapeBinPDF on a TH1.

    If `trigeff` equals 2016, 2017, or 2018, the bkg trigger efficiency as a 
    function of mT_AK15_subl is prefixed to the expression.
    """
    if pdf_type not in {'main', 'alt'}: raise Exception('Unknown pdf_type %s' % pdf_type)
    if name is None: name = uid()
    logger.info(
        'Building name={} pdf_type={} n_pars={} mt.GetName()="{}", bkg_th1.GetName()="{}"'
        .format(name, pdf_type, n_pars, mt.GetName(), bkg_th1.GetName())
        )
    expression = pdf_expression(pdf_type, n_pars, mt_scale)
    if trigeff in [2016, 2017, 2018]:
        logger.info('Adding trigger efficiency formula to expression')
        expression = '({})*({})'.format(trigeff_expression(trigeff), expression)
    parameters = pdf_parameters(pdf_type, n_pars, name)
    logger.info(
        'Expression: {}; Parameter names: {}'
        .format(expression, ', '.join(p.GetName() for p in parameters))
        )
    generic_pdf = ROOT.RooGenericPdf(
        name+'_rgp', name+'_rgp',
        expression, ROOT.RooArgList(mt, *parameters)
        )
    object_keeper.add(generic_pdf)
    parametric_shape_bin_pdf = ROOT.RooParametricShapeBinPdf(
        name+'_rpsbp', name+'_rpsbp',
        generic_pdf, mt, ROOT.RooArgList(*parameters), bkg_th1
        )
    object_keeper.add(parametric_shape_bin_pdf)
    pdf = PDF()
    pdf.name = name
    pdf.pdf = parametric_shape_bin_pdf
    pdf.rgp = generic_pdf
    pdf.expression = expression # Tag it onto the instance
    pdf.parameters = parameters
    pdf.n_pars = n_pars
    pdf.th1 = bkg_th1
    pdf.pdf_type = pdf_type
    pdf.mt = mt
    logger.info('Created {}'.format(pdf))
    return pdf


def pdfs_factory(pdf_type, mt, bkg_th1, name=None, mt_scale='1000', trigeff=None):
    """
    Like pdf_factory, but returns a list for all available n_pars
    """
    if name is None: name = uid()
    all_n_pars = [2, 3, 4, 5] if pdf_type == 'main' else [1, 2, 3, 4]
    return [ pdf_factory(pdf_type, n_pars, mt, bkg_th1, name+'_npars'+str(n_pars), mt_scale, trigeff=trigeff) for n_pars in all_n_pars]


def to_list(rooarglist):
    return [rooarglist.at(i) for i in range(rooarglist.getSize())]


def get_variables(rooabsarg):
    """
    Returns a list of all variables a RooAbsArg depends on
    """
    argset = ROOT.RooArgList(rooabsarg.getVariables())
    return [argset.at(i) for i in range(argset.getSize())]


def set_pdf_to_fitresult(pdf, res):
    """
    Sets the parameters of a pdf to the fit result. 
    """
    def set_par(par, value):
        par.setRange(value-10., value+10.)
        par.setVal(value)
    import scipy # type: ignore
    if isinstance(res, ROOT.RooFitResult):
        vals = []
        for p_fit, p_pdf in zip(to_list(res.floatParsFinal()), pdf.parameters):
            set_par(p_pdf, p_fit.getVal())
            vals.append(p_fit.getVal())
        return vals
    elif isinstance(res, scipy.optimize.optimize.OptimizeResult):
        for val, p_pdf in zip(res.x, pdf.parameters):
            set_par(p_pdf, val)
        return res.x


def plot_pdf_for_various_fitresults(pdf, fit_results, data_obs, outfile='test.pdf', labels=None, title=''):
    """
    Plots the fitted bkg pdfs on top of the data histogram.
    """
    # First find the mT Roo variable in one of the pdfs
    mt = pdf.parameters[0]
    mt_min = mt.getMin()
    mt_max = mt.getMax()

    # Open the frame
    xframe = mt.frame(ROOT.RooFit.Title(title))
    c1 = ROOT.TCanvas(str(uuid.uuid4()), '', 1000, 800)
    c1.cd()

    # Plot the data histogram
    data_obs.plotOn(xframe, ROOT.RooFit.Name("data_obs"))
    norm = data_obs.sumEntries()

    # Plot the pdfs (its parameters already at fit result)
    colors = [ROOT.kPink+6, ROOT.kBlue-4, ROOT.kRed-4, ROOT.kGreen+1]

    py_chi2 = build_chi2(pdf.expression, data_obs.createHistogram(mt.GetName()))

    base_pdf = pdf
    for i, res in enumerate(fit_results):
        pdf = rebuild_rpsbp(base_pdf)
        vals = set_pdf_to_fitresult(pdf, res)
        logger.info(
            'i=%s; Manual chi2=%.5f, chi2_via_frame=%.5f',
            i, py_chi2(vals), get_chi2_viaframe(mt, pdf.pdf, data_obs, len(vals))[1]
            )
        pdf.plotOn(
            xframe,
            ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),
            ROOT.RooFit.LineColor(colors[i]),
            # ROOT.RooFit.FillColor(ROOT.kOrange),
            ROOT.RooFit.FillStyle(1001),
            ROOT.RooFit.DrawOption("L"),
            ROOT.RooFit.Name(pdf.GetName()),
            ROOT.RooFit.Range("Full")
            )
        chi2 = xframe.chiSquare(pdf.GetName(), "data_obs", len(pdf.parameters)-1)
        par_value_str = ', '.join(['p{}={:.3f}'.format(iv, v) for iv, v in enumerate(vals)])
        label = labels[i] if labels else 'fit'+str(i)
        txt = ROOT.TText(
            .13, 0.13+i*.045,
            "{}, chi2={:.4f}, {}".format(label, chi2, par_value_str)
            )
        txt.SetNDC()
        txt.SetTextSize(0.03)
        txt.SetTextColor(colors[i])
        xframe.addObject(txt) 
        txt.Draw()

    xframe.SetMinimum(0.002)
    xframe.Draw()
    c1.SetLogy()
    c1.SaveAs(outfile)
    if outfile.endswith('.pdf'): c1.SaveAs(outfile.replace('.pdf', '.png'))
    del xframe, c1


def plot_fits(pdfs, fit_results, data_obs, outfile='test.pdf'):
    """
    Plots the fitted bkg pdfs on top of the data histogram.
    """
    # First find the mT Roo variable in one of the pdfs
    mT = pdfs[0].mt
    mT_min = mT.getMin()
    mT_max = mT.getMax()

    # Open the frame
    xframe = mT.frame(ROOT.RooFit.Title("extended ML fit example"))
    c1 = ROOT.TCanvas()
    c1.cd()

    # Plot the data histogram
    data_obs.plotOn(xframe, ROOT.RooFit.Name("data_obs"))

    # Set to fitresult
    for pdf, res in zip(pdfs, fit_results): set_pdf_to_fitresult(pdf, res)

    # Plot the pdfs
    colors = [ROOT.kPink+6, ROOT.kBlue-4, ROOT.kRed-4, ROOT.kGreen+1]
    colors.extend(colors)
    colors.extend(colors)
    for pdf, color in zip(pdfs, colors):
        pdf.pdf.plotOn(
            xframe,
            ROOT.RooFit.Name(pdf.pdf.GetName()),
            ROOT.RooFit.LineColor(color),
            ROOT.RooFit.Range("Full")
            )

    # Add the fit result text labels
    for i, fit_result in enumerate(fit_results):
        n_fit_pars = len(fit_result.floatParsFinal())
        chi2 = xframe.chiSquare(pdfs[i].pdf.GetName(), "data_obs", n_fit_pars)

        par_values = [ 'p{}={:.3f}'.format(i, v.getVal()) for i, v in enumerate(pdfs[i].parameters)]
        par_value_str = ', '.join(par_values)
        
        txt = ROOT.TText(
            .12, 0.12+i*.05,
            "model {}, nP {}, chi2: {:.4f}, {}".format(i, n_fit_pars, chi2, par_value_str)
            )
        txt.SetNDC()
        txt.SetTextSize(0.04)
        txt.SetTextColor(colors[i])
        xframe.addObject(txt) 
        txt.Draw()

    xframe.SetMinimum(0.0002)
    xframe.Draw()
    c1.SetLogy()
    c1.SaveAs(outfile)
    c1.SaveAs(outfile.replace('.pdf', '.png'))
    del xframe, c1


def pdf_ploton_frame(frame, pdf, norm):
    pdf.plotOn(
        frame,
        ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.NumEvent),
        ROOT.RooFit.LineColor(ROOT.kBlue),
        ROOT.RooFit.FillColor(ROOT.kOrange),
        ROOT.RooFit.FillStyle(1001),
        ROOT.RooFit.DrawOption("L"),
        ROOT.RooFit.Name(pdf.GetName()),
        ROOT.RooFit.Range("Full")
        )
    pdf.paramOn(
        frame,
        ROOT.RooFit.Label(pdf.GetTitle()),
        ROOT.RooFit.Layout(0.45, 0.95, 0.94),
        ROOT.RooFit.Format("NEAU")
        )


def data_ploton_frame(frame, data, is_data=True):
    data_graph = data.plotOn(
        frame,
        ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson if is_data else ROOT.RooAbsData.SumW2),
        ROOT.RooFit.DrawOption("PE0"),
        ROOT.RooFit.Name(data.GetName())
        )
    return data_graph


def get_chi2_viaframe(mt, pdf, data, n_fit_parameters):
    """
    Get the chi2 value of the fitted pdf on data by plotting it on
    a temporary frame.

    For some reason this is much faster than the non-plotting method.
    """
    logger.debug('Using plotOn residuals')
    frame = mt.frame(ROOT.RooFit.Title(""))
    pdf_ploton_frame(frame, pdf, norm=data.sumEntries())
    data_graph = data_ploton_frame(frame, data)
    roochi2 = frame.chiSquare(pdf.GetName(), data.GetName(), n_fit_parameters)
    # Number of degrees of freedom: data will contain zeros of mt binning
    # is finer than original data binning; don't count those zeros
    dhist = frame.findObject(data.GetName(),ROOT.RooHist.Class())
    n_bins = 0
    for i in range(dhist.GetN()):
        x = ROOT.Double(0.)
        y = ROOT.Double(0.)
        dhist.GetPoint(i,x,y)
        if y!=0: n_bins += 1
    ndf = n_bins - n_fit_parameters
    chi2 = roochi2 * ndf
    roopro = ROOT.TMath.Prob(chi2, ndf)
    return roochi2, chi2, roopro, ndf


def get_rss_viaframe(mt, pdf, data, norm=None, return_n_bins=False):
    """
    Get the Residual Sum of Squares (RSS) of the fitted pdf on data by plotting it on
    a temporary frame.

    if `return_n_bins` is True, also the number of bins that were used to calculate the
    RSS.

    For some reason this is much faster than the non-plotting method.
    """
    logger.info('Calculating RSS using plotOn residuals')
    rss = 0.
    frame = mt.frame(ROOT.RooFit.Title(""))
    pdf_ploton_frame(frame, pdf, norm=(data.sumEntries() if norm is None else norm))
    data_graph = data_ploton_frame(frame, data)

    hist = data_graph.getHist()
    residuals = frame.residHist(data.GetName(), pdf.GetName(), False, True) # this is y_i - f(x_i)
    xmin, xmax = array('d', [0.]), array('d', [0.])
    data.getRange(mt, xmin, xmax)

    n_bins = 0
    for i in xrange(0, hist.GetN()): # type:ignore
        x, y = hist.GetX()[i], hist.GetY()[i]
        res_y = residuals.GetY()[i]
        left  = x - hist.GetErrorXlow(i)
        right = x + hist.GetErrorXhigh(i)
        if left > xmax[0] and right > xmax[0]: continue
        elif y <= 0.: continue
        if logger.level <= logging.DEBUG:
            y_pdf = y - res_y
            logger.debug(
                '{i} ({left:.2f} to {right:.2f}):'
                # '\n  pdf  : {val_pdf:8.3f}'
                '\n  data : {y:8.3f}'
                '\n  residual : {res_y:8.3f}'
                '\n  pdf : {y_pdf:8.3f}'
                .format(**locals())
                )
        rss += res_y**2
        n_bins += 1
    rss = sqrt(rss)
    logger.info('rss_viaframe: {}'.format(rss))
    return (rss, n_bins) if return_n_bins else rss


def do_fisher_test(mt, data, pdfs, a_crit=.07):
    """
    Does a Fisher test. First computes the cl_vals for all combinations
    of pdfs, then picks the winner.
    
    Returns the pdf that won.
    """
    rsss = [ get_rss_viaframe(mt, pdf.pdf, data, return_n_bins=True) for pdf in pdfs ]
    # Compute test values of all combinations beforehand
    cl_vals = {}
    for i, j in itertools.combinations(range(len(pdfs)), 2):
        n1 = pdfs[i].n_pars
        n2 = pdfs[j].n_pars
        rss1, _      = rsss[i]
        rss2, n_bins = rsss[j]
        f = ((rss1-rss2)/(n2-n1)) / (rss2/(n_bins-n2))
        cl = 1.-ROOT.TMath.FDistI(f, n2-n1, n_bins-n2)
        cl_vals[(i,j)] = cl
    # Get the winner index
    get_winner = lambda i, j: i if cl_vals[(i,j)] > a_crit else j
    winner = get_winner(0,1)
    for i in range(2,len(pdfs)):
        winner = get_winner(winner, i)
    if logger.level <= logging.INFO:
        # Print the table
        logger.info(
            'Winner is pdf {} with {} parameters'
            .format(winner, pdfs[winner].n_pars)
            )
        table = [[''] + range(1,len(pdfs))]
        for i in range(len(pdfs)-1):
            table.append(
                [i] + ['{:6.4f}'.format(cl_vals[(i,j)]) for j in range(i+1,len(pdfs))]
                )
        logger.info('alpha values of pdf i vs j:\n' + tabelize(table))
    return winner
    

# _______________________________________________________________________
# For combine

def gen_datacard(
    jsonfile, bdtcut, signal,
    lock=None, injectsignal=False,
    tag=None, mt_min=180., mt_max=720.,
    trigeff=None,
    ):
    mz = signal['mz']
    rinv = signal['rinv']

    input = InputData(jsonfile)
    input = input.cut_mt(mt_min, mt_max)

    bdt_str = bdtcut.replace('.', 'p')
    mt = get_mt(input.mt[0], input.mt[-1], input.n_bins, name='mt')
    bkg_th1 = input.bkg_th1('bkg', bdtcut)

    data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_th1, 1.)

    pdfs_dict = {
        'main' : pdfs_factory('main', mt, bkg_th1, name='bsvj_bkgfitmain', trigeff=trigeff),
        'alt' : pdfs_factory('alt', mt, bkg_th1, name='bsvj_bkgfitalt', trigeff=trigeff),
        }
    winner_pdfs = []

    from fit_cache import FitCache
    cache = FitCache(lock=lock)

    for pdf_type in ['main', 'alt']:
        pdfs = pdfs_dict[pdf_type]
        ress = [ fit(pdf, cache=cache) for pdf in pdfs ]
        i_winner = do_fisher_test(mt, data_datahist, pdfs)
        winner_pdfs.append(pdfs[i_winner])
        # plot_fits(pdfs, ress, data_datahist, pdf_type + '.pdf')

    systs = [
        ['lumi', 'lnN', 1.016, '-'],
        # Place holders
        ['trigger', 'lnN', 1.02, '-'],
        ['pdf', 'lnN', 1.05, '-'],
        ['mcstat', 'lnN', 1.07, '-'],
        ]

    sig_name = 'mz{:.0f}_rinv{:.1f}'.format(mz, rinv)
    sig_th1 = input.sig_th1(sig_name, bdtcut, mz, rinv)
    sig_datahist = ROOT.RooDataHist(sig_name, sig_name, ROOT.RooArgList(mt), sig_th1, 1.)

    # Some checks
    # assert bkg_th1.GetNbinsX() == sig_th1.GetNbinsX()
    # assert bkg_th1.GetBinLowEdge(1) == sig_th1.GetBinLowEdge(1)
    # n = bkg_th1.GetNbinsX()
    # assert bkg_th1.GetBinLowEdge(n+1) == sig_th1.GetBinLowEdge(n+1)
    # assert sig_th1.GetBinLowEdge(n+1) == mt.getMax()
    # x_sig_datahist, y_sig_datahist = roodataset_values(sig_datahist)
    # np.testing.assert_almost_equal(x_sig_datahist, input.mt_centers)
    # np.testing.assert_almost_equal(y_sig_datahist, input.sighist(bdtcut, mz).vals, decimal=3)
    # print('All passed')
    # return

    if injectsignal:
        logger.info('Injecting signal in data_obs')
        data_datahist = ROOT.RooDataHist("data_obs", "Data", ROOT.RooArgList(mt), bkg_th1+sig_th1, 1.)

    outfile = strftime('dc_%b%d{}/dc_mz{}_rinv{:.1f}_bdt{}.txt'.format('_'+tag if tag else '', mz, rinv, bdt_str))
    if injectsignal: outfile = outfile.replace('.txt', '_injectsig.txt')
    compile_datacard_macro(
        winner_pdfs, data_datahist, sig_datahist,
        outfile,
        systs=systs
        )


class Datacard:
    def __init__(self):
        self.shapes = [] # Expects len-4 or len-5 iterables as elements
        self.channels = [] # Expects len-2 iterables as elements, (name, rate)
        self.rates = OrderedDict() # Expects str as key, OrderedDict as value
        self.systs = []

    def __eq__(self, other):
        return (
            self.shapes == other.shapes
            and self.channels == other.channels
            and self.rates == other.rates
            and self.systs == other.systs
            )

    @property
    def syst_names(self):
        return [s[0] for s in self.systs]

    def syst_rgx(self, rgx):
        """
        Returns a list of all systematics that match a pattern
        (Uses Unix file-like pattern matching, e.g. 'bkg_*')
        """
        import fnmatch
        return [s for s in self.syst_names if fnmatch.fnmatch(s, rgx)]


def read_dc(datacard):
    """
    Returns a Datacard object based on the txt stored in the passed path
    """
    with open(datacard, 'r') as f:
        return read_dc_txt(f.read())


def read_dc_txt(txt):
    """
    Returns a Datacard object based on the passed datacard txt.
    """
    lines = txt.split('\n')
    dc = Datacard()

    blocks = []
    block = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        elif line.startswith('#'):
            continue
        elif line.startswith('---------------'):
            blocks.append(block)
            block = []
            continue
        block.append(line)
    blocks.append(block)
    if len(blocks) != 5:
        raise Exception('Found {} blocks, expected 5'.format(len(blocks)))

    # Shapes are easy
    for line in blocks[1]: dc.shapes.append(line.split()[1:])
    # pprint(dc.shapes)

    # Channels
    channel_table = transpose([ l.split() for l in blocks[2] ])[1:]
    for name, rate in channel_table:
        dc.channels.append((name, int(rate)))
    # pprint(dc.channels)

    # Rates
    rate_table = transpose([ l.split() for l in blocks[3] ])[1:]
    for ch, name, index, rate in rate_table:
        dc.rates.setdefault(ch, OrderedDict())
        dc.rates[ch][name] = int(rate)
    # pprint(dc.rates)

    # Systs
    for line in blocks[4]:
        syst = line.split()
        if len(syst) >= 3:
            try:
                syst[2] = float(syst[2])
            except TypeError:
                pass
        dc.systs.append(syst)
    # pprint(dc.systs)

    return dc


def parse_dc(dc):
    '''
    Very basic datacard formatter
    '''
    txt = ''
    line = '\n' + '-'*83

    txt += (
        'imax {} number of channels'
        '\njmax * number of backgrounds'
        '\nkmax * number of nuisance parameters'
        ).format(len(dc.channels))
    txt += line

    # Careful with workspace path: Should be relative path to DC
    shapes = copy.copy(dc.shapes)
    for i in range(len(shapes)):
        shapes[i][2] = osp.basename(shapes[i][2])

    txt += '\n' + tabelize([['shapes']+s for s in shapes])
    txt += line
    txt += '\n' + tabelize(transpose([('bin', 'observation')] + list(dc.channels)))
    txt += line

    # Format the bin/process table
    # Format per column, and only transpose at str format time
    proc_nr_dict = {}
    proc_nr_counter = [0]
    def proc_nr(proc):
        if not proc in proc_nr_dict:
            proc_nr_dict[proc] = proc_nr_counter[0]
            proc_nr_counter[0] += 1
        return proc_nr_dict[proc]
    table = [['bin', 'process', 'process', 'rate']]
    for bin in dc.rates:
        for proc in dc.rates[bin]:
            table.append([bin, proc, proc_nr(proc), int(dc.rates[bin][proc])])
    txt += '\n' + tabelize(transpose(table))

    txt += line
    txt += '\n' + tabelize(dc.systs)
    txt += '\n'
    return txt


def transpose(l):
    '''Transposes a list of lists'''
    return map(list, zip(*l))


def tabelize(data):
    '''
    Formats a list of lists to a single string (space separated).
    Rows need not be of same length.
    '''
    # Ensure data is strings
    data = [ [ str(i) for i in row ] for row in data ]
    # Determine the row with the most columns
    n_columns = max(map(len, data))
    # Determine how wide each column should be (max)
    col_widths = [0 for i in range(n_columns)]
    for row in data:
        for i_col, item in enumerate(row):
            if len(item) > col_widths[i_col]:
                col_widths[i_col] = len(item)
    # Format
    return '\n'.join(
        ' '.join(
            format(item, str(w)) for item, w in zip(row, col_widths)
            )
        for row in data
        )


def make_multipdf(pdfs, name='roomultipdf'):
    cat = ROOT.RooCategory('pdf_index', "Index of Pdf which is active")
    pdf_arglist = ROOT.RooArgList()
    for pdf in pdfs: pdf_arglist.add(pdf.pdf)
    multipdf = ROOT.RooMultiPdf(name, "All Pdfs", cat, pdf_arglist)
    multipdf.cat = cat
    multipdf.pdfs = pdfs
    norm = ROOT.RooRealVar(name+'_norm', "Number of background events", 1.0, 0., 1.e6)
    object_keeper.add_multiple([cat, norm, multipdf])
    return multipdf, norm


def compile_datacard_macro(bkg_pdf, data_obs, sig, outfile='dc_bsvj.txt', systs=None):
    do_syst = systs is not None
    w = ROOT.RooWorkspace("SVJ", "workspace")

    def commit(thing, *args, **kwargs):
        name = thing.GetName() if hasattr(thing, 'GetName') else '?'
        logger.info('Importing {} ({})'.format(name, thing))
        getattr(w, 'import')(thing, *args, **kwargs)

    # Bkg pdf: May be multiple
    is_multipdf = hasattr(bkg_pdf, '__len__')
    if is_multipdf:
        multipdf, norm = make_multipdf(bkg_pdf)
        commit(multipdf.cat)
        commit(norm)
        commit(multipdf)
    else:
        commit(bkg_pdf, ROOT.RooFit.RenameVariable(bkg_pdf.GetName(), 'bkg'))

    commit(data_obs)
    commit(sig, ROOT.RooFit.Rename('sig'))

    wsfile = outfile.replace('.txt', '.root')
    dump_ws_to_file(wsfile, w)

    # Write the dc
    dc = Datacard()
    dc.shapes.append(['roomultipdf' if is_multipdf else 'bkg', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.shapes.append(['sig', 'bsvj', wsfile, 'SVJ:$PROCESS'] + (['SVJ:$PROCESS_$SYSTEMATIC'] if do_syst else []))
    dc.shapes.append(['data_obs', 'bsvj', wsfile, 'SVJ:$PROCESS'])
    dc.channels.append(('bsvj', int(data_obs.sumEntries())))
    dc.rates['bsvj'] = OrderedDict()
    dc.rates['bsvj']['sig'] = sig.sumEntries()
    dc.rates['bsvj']['roomultipdf' if is_multipdf else 'bkg'] = data_obs.sumEntries()
    # Freely floating bkg parameters
    def systs_for_pdf(pdf):
        for par in pdf.parameters:
            dc.systs.append([par.GetName(), 'flatParam'])
    [systs_for_pdf(p) for p in multipdf.pdfs] if is_multipdf else systs_for_pdf(bkg_pdf)
    # Rest of the systematics
    if is_multipdf: dc.systs.append([multipdf.cat.GetName(), 'discrete'])
    if do_syst: dc.systs.extend(systs)
    txt = parse_dc(dc)

    logger.info('txt datacard:\n%s', txt)
    logger.info('Dumping txt to ' + outfile)
    if not osp.isdir(osp.dirname(outfile)): os.makedirs(osp.dirname(outfile))
    with open(outfile, 'w') as f:
        f.write(txt)


def camel_to_snake(name):
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


class CombineCommand(object):

    comma_separated_args = [
        '--freezeParameters',
        '--trackParameters',
        '--trackErrors',
        ]
    comma_separated_arg_map = { camel_to_snake(v.strip('-')) : v for v in comma_separated_args }
    comma_separated_arg_map['redefine_signal_pois'] = '--redefineSignalPOIs'

    def __init__(self, dc=None, method='MultiDimFit', args=None, kwargs=None):
        self.dc = dc
        self.method = method
        self.args = set() if args is None else args
        self.kwargs = OrderedDict()
        if kwargs: self.kwargs.update(kwargs)
        for key in self.comma_separated_arg_map: setattr(self, key, [])
        self.parameters = OrderedDict()
        self.parameter_ranges = OrderedDict()

    def get_name_key(self):
        """
        'name' parameter for the combine CLI can be either '-n' or '--name';
        ensure consistency
        """
        assert ('-n' in self.kwargs) + ('--name' in self.kwargs) < 2
        if '-n' in self.kwargs:
            return '-n'
        else:
            return '--name'
            
    @property
    def name(self):
        return self.kwargs.get(self.get_name_key(), '')

    @name.setter
    def name(self, val):
        self.kwargs[self.get_name_key()] = val

    @property
    def outfile(self):
        return 'higgsCombine{0}.{1}.mH120.root'.format(self.name, self.method)

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self):
        return '<CombineCommand ' + '\n    '.join(self.parse()) + '\n    >'

    @property
    def str(self):
        return ' '.join(self.parse())

    def add_range(self, name, left, right):
        self.parameter_ranges[name] = [left, right]

    def set_parameter(self, name, value, left=None, right=None):
        self.parameters[name] = value
        if left is not None and right is not None: self.add_range(name, left, right)

    def parse(self):
        """
        Returns the command as a list
        """
        command = ['combine']
        command.append('-M ' + self.method)
        if not self.dc: raise Exception('self.dc must be a valid path')
        command.append(self.dc)
        command.extend(list(self.args))
        command.extend([k+' '+str(v) for k, v in self.kwargs.items()])

        for attr, command_str in self.comma_separated_arg_map.items():
            values = getattr(self, attr)
            if not values: continue
            command.append(command_str + ' ' + ','.join(values))
        
        if self.parameters:
            strs = ['{0}={1}'.format(k, v) for k, v in self.parameters.items()]
            command.append('--setParameters ' + ','.join(strs))

        if self.parameter_ranges:
            strs = ['{0}={1},{2}'.format(parname, *ranges) for parname, ranges in self.parameter_ranges.items()]
            command.append('--setParameterRanges ' + ':'.join(strs))

        return command


def likelihood_scan_factory(
    datacard,
    rmin=0., rmax=2., n=40,
    verbosity=0, asimov=False,
    pdf_type='alt'
    ):
    """
    Returns a good CombineCommand template for a likelihood scan 
    """
    dc = read_dc(datacard)
    cmd = CombineCommand(datacard, 'MultiDimFit')

    cmd.redefine_signal_pois.append('r')
    cmd.add_range('r', rmin, rmax)
    cmd.track_parameters.extend(['r'])

    cmd.args.add('--saveWorkspace')
    cmd.args.add('--saveNLL')
    cmd.kwargs['--algo'] = 'grid'
    cmd.kwargs['--points'] = n
    cmd.kwargs['--X-rtd'] = 'REMOVE_CONSTANT_ZERO_POINT=1'
    cmd.kwargs['--alignEdges'] = 1
    cmd.kwargs['-v'] = verbosity

    if asimov:
        cmd.kwargs['-t'] = '-1'
        cmd.args.add('--toysFreq')
        cmd.kwargs['-n'] = 'Asimov'
    else:
        cmd.kwargs['-n'] = 'Observed'

    cmd.freeze_parameters.append('pdf_index')
    if pdf_type == 'alt':
        cmd.set_parameter('pdf_index', 1)
        cmd.freeze_parameters.extend(dc.syst_rgx('bsvj_bkgfitmain_*'))
    elif pdf_type == 'main':
        cmd.set_parameter('pdf_index', 0)
        cmd.freeze_parameters.extend(dc.syst_rgx('bsvj_bkgfitalt_*'))
    else:
        raise Exception('Unknown pdf_type {}'.format(pdf_type))

    return cmd


@contextmanager
def switchdir(other_dir):
    if other_dir:
        try:
            current_dir = os.getcwd()
            logger.info('Changing into %s', other_dir)
            os.chdir(other_dir)
            yield other_dir
        finally:
            logger.info('Changing back into %s', current_dir)
            os.chdir(current_dir)
    else:
        try:
            yield None
        finally:
            pass


def run_command(cmd, chdir=None):
    with switchdir(chdir):
        logger.warning('Issuing command: ' + cmd)
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=True,
            )
        output = []
        for stdout_line in iter(process.stdout.readline, ""):
            subprocess_logger.info(stdout_line.rstrip("\n"))
            output.append(stdout_line)
        process.stdout.close()
        process.wait()
        returncode = process.returncode

        if returncode == 0:
            logger.info("Command exited with status 0 - all good")
        else:
            logger.error("Exit status {0} for command: {1}".format(returncode, cmd))
            raise subprocess.CalledProcessError(cmd, returncode)
        return output


def run_combine_command(cmd, chdir=None):
    if chdir:
        # Fix datacard to be an absolute path first
        cmd = cmd.copy()
        cmd.dc = osp.abspath(cmd.dc)
    logger.info('Running {0}'.format(cmd))
    return run_command(cmd.str, chdir)


# _______________________________________________________________________
# ROOT and RooFit interface

@contextmanager
def open_root(path, mode='READ'):
    '''Context manager that takes care of closing the ROOT file properly'''
    try:
        tfile = ROOT.TFile.Open(path, mode)
        yield tfile
    finally:
        tfile.Close()


class ROOTObjectKeeper:
    """
    Keeps ROOT objects in it so ROOT doesn't garbage collect them.
    """
    def __init__(self):
        self.objects = {}
    
    def add(self, thing):
        try:
            key = thing.GetName()
        except AttributeError:
            key = str(uuid.uuid4())
        if key in self.objects:
            logger.warning('Variable %s (%s) already in object keeper', thing.GetName(), thing.GetTitle())
        self.objects[key] = thing

    def add_multiple(self, things):
        for t in things: self.add(t)


object_keeper = ROOTObjectKeeper()


def get_arrays(rootfile, treename='limit'):
    """Poor man's uproot: make a dict `branch_name : np.array of values`"""
    with open_root(rootfile) as f:
        tree = f.Get(treename)
        branch_names = [ b.GetName() for b in tree.GetListOfBranches() ]
        r = { b : [] for b in branch_names }
        for entry in tree:
            for b in branch_names:
                r[b].append(getattr(entry, b))
    return {k : np.array(v) for k, v in r.items()}


def get_ws(f, wsname=None):
    """
    Functionality to grab the first workspace that's encountered in a rootfile
    """
    if wsname is None:
        # Pick the first one
        for key in f.GetListOfKeys():
            obj = f.Get(key.GetName())
            if isinstance(obj, ROOT.RooWorkspace):
                ws = obj
                break
        else:
            f.ls()
            raise Exception('No workspace found in {0}'.format(f))
    else:
        ws = f.Get(wsname)
    return ws


def th1_to_hist(h):
    n_bins = h.GetNbinsX()
    # GetBinLowEdge of the right overflow bin is the high edge of the actual last bin
    return AttrDict(
        binning = np.array([h.GetBinLowEdge(i) for i in range(1,n_bins+2)]),
        vals = np.array([h.GetBinContent(i) for i in range(1,n_bins+1)]),
        errs = np.array([h.GetBinError(i) for i in range(1,n_bins+1)])
        )


def th1_binning_and_values(h, return_errs=False):
    """
    Returns the binning and values of the histogram.
    Does not include the overflows.
    """
    import inspect
    logger.warning('DEPRECATED: Use th1_to_hist instead; called by {}'.format(inspect.stack()[1][3]))
    n_bins = h.GetNbinsX()
    # GetBinLowEdge of the right overflow bin is the high edge of the actual last bin
    binning = np.array([h.GetBinLowEdge(i) for i in range(1,n_bins+2)])
    values = np.array([h.GetBinContent(i) for i in range(1,n_bins+1)])
    errs = np.array([h.GetBinError(i) for i in range(1,n_bins+1)])
    return (binning, values, errs) if return_errs else (binning, values)


def th1_to_datahist(histogram, mt=None):
    if mt is None: mt = get_mt_from_th1(histogram)
    datahist = ROOT.RooDataHist(uid(), '', ROOT.RooArgList(mt), histogram, 1.)
    datahist.mt = mt
    return datahist


def binning_from_roorealvar(x):
    binning = [x.getMin()]
    for i in range(x.numBins()):
        binning.append(binning[-1] + x.getBinWidth(i))
    return np.array(binning)


def roodataset_values(data, varname='mt'):
    """
    Works on both RooDataHist and RooDataSet!
    """
    x = []
    y = []
    for i in range(data.numEntries()):
        s = data.get(i)
        x.append(s[varname].getVal())
        y.append(data.weight())
    return np.array(x), np.array(y)


def pdf_values(pdf, x_vals, varname='mt'):
    """
    Equivalent to the y_values of pdf.createHistogram('...', mt)
    """
    variable = pdf.getVariables()[varname]
    y = []
    for x in x_vals:
        variable.setVal(x)
        # logger.info('{}: set {} to {}'.format(pdf.GetName(), variable.GetName(), x))
        # logger.info('  got pdf = {}'.format(pdf.getVal()))
        y.append(pdf.getVal())
    y = np.array(y)
    return y / (y.sum() if y.sum()!=0. else 1.)



# # DEPRECATED: Fit cache

# FIT_CACHE_FILE = 'fit_cache.pickle'

# def _read_fit_cache():
#     import pickle
#     if osp.isfile(FIT_CACHE_FILE):
#         logger.info('Reading cached fits from %s', FIT_CACHE_FILE)
#         with open(FIT_CACHE_FILE, 'rb') as f:
#             return pickle.load(f)
#     else:
#         return {}

# def _write_fit_cache(cache_dict):
#     import pickle
#     with open(FIT_CACHE_FILE, 'wb') as f:
#         pickle.dump(cache_dict, f)

# def _get_from_fit_cache(fit_hash):
#     return _read_fit_cache().get(fit_hash, None)

# def _add_one_to_fit_cache(key, result):
#     d = _read_fit_cache()
#     d[key] = result
#     _write_fit_cache(d)


# # DEPRECATED: Old fitting functions

# def _fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=None, hash=None, **minimize_kwargs):
#     """
#     The actual entry point to the scipy fit
#     """
#     n_fit_pars = count_parameters(expression) - 1 # -1 because par 0 is mT
#     logger.info('Fitting {0} with {1} parameters'.format(expression, n_fit_pars))
#     from scipy.optimize import minimize
#     chi2 = build_chi2(expression, histogram)
#     if init_vals is None: init_vals = np.ones(n_fit_pars)
#     res = minimize(chi2, init_vals, **minimize_kwargs)
#     # Save some extra information in the result
#     res.x_init = np.array(init_vals)
#     res.expression = expression
#     res.hash = hash
#     # Set approximate uncertainties; see https://stackoverflow.com/a/53489234
#     # Assume ftol ~ function value
#     # try:
#     #     res.dx = np.sqrt(res.fun * np.diagonal(res.hess_inv))
#     # except:
#     #     logger.error('Failed to set uncertainties; using found function values as proxies')
#     #     res.dx = res.x.copy()
#     return res


# def brute_force_init_vals(npars, values):
#     return np.array(list(itertools.product(*[values for i in range(npars)])))


# def fit_expr_to_histogram_robust(expression, histogram):
#     """
#     Heuristic around `fit_pdf_expression_to_histogram_python`.
#     First attempts a single fit, and only goes for the bruteforce if the single fit
#     did not converge properly.
#     """
#     fit_hash = make_fit_hash(expression, histogram)
#     res = _get_from_fit_cache(fit_hash)
#     if res is not None:
#         logger.warning('Returning cached robust fit')
#         return res
#     res = fit_pdf_expression_to_histogram_python(
#         expression, histogram,
#         tol=1e-3, method='BFGS'
#         )
#     # Refit with output from first fit
#     res = fit_pdf_expression_to_histogram_python(
#         expression, histogram,
#         init_vals=res.x,
#         tol=1e-6, method='Nelder-Mead'
#         )
#     if not res.success:
#         logger.info('Fit did not converge with single try; brute forcing it')
#         results = []
#         for method in ['BFGS', 'Nelder-Mead']:
#             results = fit_pdf_expression_to_histogram_python(
#                 expression, histogram,
#                 init_vals=brute_force_init_vals(count_parameters(expression)-1, [-1., 1.]),
#                 tol=1e-3, method=method
#                 )
#         results = [ r for r in results if not(np.isnan(r.fun) or np.isposinf(r.fun) or np.isneginf(r.fun)) ]
#         if len(results) == 0: raise Exception('Not a single fit of the brute force converged!')
#         i_min = np.argmin([r.fun for r in results])        
#         res = results[i_min]
#     logger.info('Final scipy fit result:\n%s', res)
#     _add_one_to_fit_cache(fit_hash, res)
#     return res


# def fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=None, cache=True, cache_dict=None, **minimize_kwargs):
#     """
#     Fits a a background pdf expression to a TH1 histogram with scipy.optimize.minimize.
#     Assumes @0 in the expression is mT, and the histogram is binned in mT.

#     If `cache` is True, it tries to save the fit result to a file. If init_vals is
#     a 2-dimensional array, the fit is repeated for each initial value.

#     If `cache_dict` is given, no read/write to the cache file is performed.
#     """
#     _write_cache = False
#     if cache and cache_dict is None:
#         # If cache is enabled, and no cache dict was specified, treat this as one fit result
#         # and do the read/write IO of the cache file
#         cache_dict = _read_fit_cache()
#         _write_cache = True
    
#     if init_vals is not None:
#         init_vals = np.array(init_vals)
#         if len(init_vals.shape) > 1:
#             logger.info('Will run fit for %s different initial values', init_vals.shape[0])
#             results = [
#                 fit_pdf_expression_to_histogram_python(
#                     expression, histogram, init_vals=x_init,
#                     cache=cache, cache_dict=cache_dict, **minimize_kwargs
#                     ) for x_init in init_vals
#                 ]
#             _write_fit_cache(cache_dict)
#             return results

#     fit_hash = make_fit_hash(expression, histogram, init_vals, **minimize_kwargs)

#     if cache and fit_hash in cache_dict:
#         # Nothing new to save, so return immediately
#         logger.warning('Returning cached fit')
#         return cache_dict[fit_hash]
#     else:
#         res = _fit_pdf_expression_to_histogram_python(expression, histogram, init_vals=init_vals, hash=fit_hash, **minimize_kwargs)
#         cache_dict[fit_hash] = res
#         if _write_cache: _write_fit_cache(cache_dict)
#         return res
