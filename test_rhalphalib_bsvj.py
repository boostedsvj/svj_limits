from __future__ import print_function, division
import sys
import os
from contextlib import contextmanager
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import ROOT

rl.util.install_roofit_helpers()
#rl.ParametericSample.PreferRooParametricHist = False


def expo_sample(norm, scale, obs, loc=0):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning, loc=loc) * norm
    return (np.diff(cdf), obs.binning, obs.name)


def gaus_sample(norm, loc, scale, obs):
    cdf = scipy.stats.norm.cdf(loc=loc, scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)


def obs_to_TH1(channel):
    sumw = channel.getObservation()
    return rl.util._to_TH1((sumw,sumw), channel.observable.binning, channel.observable.name)


@contextmanager
def quick_ax(outfile, figsize=(12,12)):
    import matplotlib.pyplot as plt
    try:
        fig = plt.figure(figsize=figsize)
        ax = fig.gca()
        yield ax
    finally:
        plt.savefig(outfile, bbox_inches='tight')


def set_mpl_fontsize(small=22, medium=28, large=32, legend=None):
    import matplotlib.pyplot as plt
    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', titlesize=small)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=medium if legend is None else legend)    # legend fontsize
    plt.rc('figure', titlesize=large)  # fontsize of the figure title
set_mpl_fontsize()


def get_color_cycle():
    import matplotlib.pyplot as plt
    from itertools import cycle
    colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    return colors


def plot_tf(mtbins, fail_bkg, pass_bkg, fitresult):
    mtpts = mtbins[:-1] + 0.5 * np.diff(mtbins)
    mtscaled = (mtpts - min(mtbins)) / (max(mtbins) - min(mtbins))

    bkg_eff = pass_bkg.Integral()/fail_bkg.Integral()
    tfs_th1 = pass_bkg
    tfs_th1.Divide(fail_bkg)
    tfs = rl.util._to_numpy(tfs_th1, read_sumw2=True)
    errs = np.sqrt(tfs[3])

    tf_name = "tf_mc"
    npar = len([f for f in fitresult.floatParsFinal() if tf_name in f.GetName()])-1
    tf_mc = rl.BasisPoly(tf_name, (npar,), ["mt"])
    tf_mc.update_from_roofit(fitresult)
    tf_mc_vals, tf_mc_band = tf_mc(mtscaled, nominal=True, errorband=True)

    outfile = "tf_test_bsvj.png"
    with quick_ax(outfile=outfile) as ax:
        colors = get_color_cycle()
        pcolor = next(colors)
        ax.errorbar(mtpts, tfs[0], yerr=errs, label="MC", color=pcolor)
        pcolor = next(colors)
        ax.plot(mtpts, bkg_eff * tf_mc_vals, label="fit", color=pcolor)
        ax.fill_between(mtpts, bkg_eff * tf_mc_band[0], bkg_eff * tf_mc_band[1], alpha=0.2, color=pcolor)
        ax.legend(fontsize=18, framealpha=0.0)
        ax.set_xlabel(r'$m_{\mathrm{T}}$ [GeV]')
        ax.set_ylabel(f'TF')


def plot_hist(th1, ax, **kwargs):
    arrs = rl.util._to_numpy(th1, read_sumw2=True)
    def get_kwargs(orig, keys):
        return {k : orig.get(k, None) for k in keys}
    step_keys = ['where','label','color']
    ax.step(arrs[1][:-1], arrs[0], **get_kwargs(kwargs, step_keys))
    fill_keys = ['where','alpha','color']
    fill_dict = get_kwargs(kwargs, fill_keys)
    fill_dict['step'] = fill_dict.pop('where')
    errs = np.sqrt(arrs[3])
    arr_dn = arrs[0]-errs
    arr_up = arrs[0]+errs
    ax.fill_between(arrs[1][:-1], arr_dn, arr_up, **fill_dict)


def plot_sr_cr(fail_bkg, pass_bkg):
    outfile = "srcr_test_bsvj.png"
    ranges = []
    with quick_ax(outfile=outfile) as ax:
        colors = get_color_cycle()
        pcolor = next(colors)
        plot_hist(fail_bkg, ax, where='post', label="fail", alpha=0.2, color=pcolor)
        pcolor = next(colors)
        plot_hist(pass_bkg, ax, where='post', label="pass", alpha=0.2, color=pcolor)
        ax.legend(fontsize=18, framealpha=0.0)
        ax.set_xlabel(r'$m_{\mathrm{T}}$ [GeV]')
        ax.set_ylabel(f'Number of events')
        ax.set_yscale('log')


def test_rhalphabet(tmpdir):
    throwPoisson = False
    tfFromMC = True

    jec = rl.NuisanceParameter("CMS_jec", "lnN")
    massScale = rl.NuisanceParameter("CMS_massScale", "shape")
    lumi = rl.NuisanceParameter("CMS_lumi", "lnN")

    mtmin = 180.
    mtmax = 650.
    mtstep = 10.
    mtbins = np.arange(mtmin, mtmax+mtstep, mtstep)
    nmt = len(mtbins) - 1
    mt = rl.Observable("mt", mtbins)

    mtpts = mtbins[:-1] + 0.5 * np.diff(mtbins)
    mtscaled = (mtpts - mtmin) / (mtmax - mtmin)

    msig = 350
    template_info = {
        "yields": {
#           "svj": {"pass": 35100, "fail": 7760},
#            "bkg": {"pass": 379780, "fail": 315180},
            "svj": {"pass": 35100, "fail": 17680},
            "bkg": {"pass": 379780, "fail": 2278568},
        },
        "locs": {
            "svj": {"pass": msig, "fail": msig},
            "bkg": {"pass": mtmin, "fail": mtmin},
        },
        "scales": {
            "svj": {"pass": 0.25*msig, "fail": msig},
#            "bkg": {"pass": 93.8, "fail": 95.9},
            "bkg": {"pass": 93.8, "fail": 90.7},
        },
    }

    if tfFromMC:
        # Build bkg MC pass+fail model and fit to polynomial
        bkgmodel = rl.Model("bkgmodel")
        for region in ["pass", "fail"]:
            ch = rl.Channel(region)
            bkgmodel.addChannel(ch)
            # mock template
            templates = {
                "bkg": expo_sample(
                    norm=template_info["yields"]["bkg"][region],
                    loc=template_info["locs"]["bkg"][region],
                    scale=template_info["scales"]["bkg"][region],
                    obs=mt,
                ),
            }
            ch.setObservation(templates["bkg"])

        bkgeff = template_info["yields"]["bkg"]["pass"] / template_info["yields"]["bkg"]["fail"]
        tf_mc = rl.BernsteinPoly("tf_mc", (2,), ["mt"], limits=(0, 10))
        tf_mc_params = bkgeff * tf_mc(mtscaled)
        failCh = bkgmodel["fail"]
        passCh = bkgmodel["pass"]
        bkgparams = np.array([rl.IndependentParameter("bkgparam_mtbin%d" % (i), 0) for i in range(mt.nbins)])
        initial_bkg = failCh.getObservation().astype(float)  # was integer, and numpy complained about subtracting float from it
        scaledparams = initial_bkg * (1 + 1.0 / np.maximum(1.0, np.sqrt(initial_bkg))) ** (bkgparams)
        fail_bkg = rl.ParametericSample("fail_bkg", rl.Sample.BACKGROUND, mt, scaledparams)
        failCh.addSample(fail_bkg)
        pass_bkg = rl.TransferFactorSample("pass_bkg", rl.Sample.BACKGROUND, tf_mc_params, fail_bkg)
        passCh.addSample(pass_bkg)

        bkgfit_ws = ROOT.RooWorkspace("bkgfit_ws")
        simpdf, obs = bkgmodel.renderRoofit(bkgfit_ws)
        bkgfit = simpdf.fitTo(
            obs,
            ROOT.RooFit.Extended(True),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Save(),
            ROOT.RooFit.Minimizer("Minuit2", "migrad"),
            ROOT.RooFit.PrintLevel(-1),
        )
        bkgfit_ws.add(bkgfit)
        if "pytest" not in sys.modules:
            bkgfit_ws.writeToFile(os.path.join(str(tmpdir), "svjModel_bkgfit.root"))
        if bkgfit.status() != 0:
            raise RuntimeError("Could not fit bkg")

        param_names = [p.name for p in tf_mc.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_mc.name + "_deco", bkgfit, param_names)
        tf_mc.parameters = decoVector.correlated_params.reshape(tf_mc.parameters.shape)
        tf_mc_params_final = tf_mc(mtscaled)

        fail_th1 = obs_to_TH1(bkgmodel["fail"])
        pass_th1 = obs_to_TH1(bkgmodel["pass"])
        plot_sr_cr(fail_th1, pass_th1)
        plot_tf(mtbins, fail_th1, pass_th1, bkgfit)

    # build actual fit model now
    model = rl.Model("svjModel")

    for region in ["pass", "fail"]:
        ch = rl.Channel(region)
        model.addChannel(ch)

        templates = {
            "svj": gaus_sample(
                norm=template_info["yields"]["svj"][region],
                loc=template_info["locs"]["svj"][region],
                scale=template_info["scales"]["svj"][region],
                obs=mt,
            ),
        }
        for sName in templates:
            # some mock expectations
            templ = templates[sName]
            stype = rl.Sample.SIGNAL if sName == "svj" else rl.Sample.BACKGROUND
            sample = rl.TemplateSample(ch.name + "_" + sName, stype, templ)

            # mock systematics
            jecup_ratio = np.random.normal(loc=1, scale=0.05, size=mt.nbins)
            massUp = np.linspace(0.9, 1.1, mt.nbins)
            massDn = np.linspace(1.2, 0.8, mt.nbins)

            # for jec we set lnN prior, shape will automatically be converted to norm systematic
            sample.setParamEffect(jec, jecup_ratio)
            sample.setParamEffect(massScale, massUp, massDn)
            sample.setParamEffect(lumi, 1.027)

            ch.addSample(sample)

        # make up a data_obs, with possibly different yield values, excluding signal
        templates = {
            "bkg": expo_sample(
                norm=template_info["yields"]["bkg"][region],
                loc=template_info["locs"]["bkg"][region],
                scale=template_info["scales"]["bkg"][region],
                obs=mt,
            ),
        }
        yields = sum(tpl[0] for tpl in templates.values())
        if throwPoisson:
            yields = np.random.poisson(yields)
        data_obs = (yields, mt.binning, mt.name)
        ch.setObservation(data_obs)

    tf_data = rl.BernsteinPoly("tf_data", (2,), ["mt"], limits=(0, 10))
    tf_data_params = tf_data(mtscaled)
    if tfFromMC:
        tf_params = bkgeff * tf_mc_params_final * tf_data_params
    else:
        tf_params = bkgeff * tf_data_params

    failCh = model["fail"]
    passCh = model["pass"]

    bkgparams = np.array([rl.IndependentParameter("bkgparam_mtbin%d" % (i), 0) for i in range(mt.nbins)])
    initial_bkg = failCh.getObservation().astype(float)  # was integer, and numpy complained about subtracting float from it
    scaledparams = initial_bkg * (1 + 1.0 / np.maximum(1.0, np.sqrt(initial_bkg))) ** (bkgparams)
    fail_bkg = rl.ParametericSample("fail_bkg", rl.Sample.BACKGROUND, mt, scaledparams)
    failCh.addSample(fail_bkg)
    pass_bkg = rl.TransferFactorSample("pass_bkg", rl.Sample.BACKGROUND, tf_params, fail_bkg)
    passCh.addSample(pass_bkg)

    with open(os.path.join(str(tmpdir), "svjModel.pkl"), "wb") as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), "svjModel"))


if __name__ == "__main__":
    if not os.path.exists("tmp"):
        os.mkdir("tmp")
    test_rhalphabet("tmp")
