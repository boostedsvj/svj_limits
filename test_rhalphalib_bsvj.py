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

#basis = "Bernstein"
basis = "Chebyshev"

def expo_sample(norm, scale, obs, loc=0, vals=None):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning, loc=loc) * norm
    if vals is None: vals = np.diff(cdf)
    return (vals, obs.binning, obs.name)


def gaus_sample(norm, loc, scale, obs):
    cdf = scipy.stats.norm.cdf(loc=loc, scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)


def obs_to_TH1(channel):
    sumw = get_obs(channel)
    return rl.util._to_TH1((sumw,sumw), channel.observable.binning, channel.observable.name)


def get_obs(channel):
    obs = channel.getObservation()
    if isinstance(obs,tuple): return obs[0]
    else: return obs


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
    tf_mc = rl.BasisPoly(tf_name, (npar,), ["mt"], basis=basis)
    tf_mc.update_from_roofit(fitresult)
    tf_mc_vals, tf_mc_band = tf_mc(mtscaled, nominal=True, errorband=True)

    outfile = "tmp/tf_test_bsvj.png"
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
    outfile = "tmp/srcr_test_bsvj.png"
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
    mtscaled = (mtpts - min(mtpts)) / (max(mtpts) - min(mtpts))

    msig = 350
    template_info = {
        "yields": {
            "svj": {"pass": 22640, "fail": 9619},
            "bkg": {"pass": 277440.1886037716, "fail": 472348.698199576},
        },
        "locs": {
            "svj": {"pass": msig, "fail": msig},
            "bkg": {"pass": mtmin, "fail": mtmin},
        },
        "scales": {
            "svj": {"pass": 0.25*msig, "fail": msig},
            "bkg": {"pass": 94.2, "fail": 92.1},
        },
        "vals": {
            "bkg": {
                "pass": np.array([27694.437526301597, 25199.442008055747, 23933.818656462245, 20964.187123346026, 18319.932810074883, 15612.968794859946, 14525.593370080926, 13237.611986682517, 13336.220533201704, 11083.544259571005, 10079.225086054299, 8634.173412132077, 7552.1853848984465, 7260.119516679086, 6812.459799276665, 6074.417684029322, 4919.891813983209, 4549.854085942265, 4125.581354004797, 3671.634960473515, 3425.505970105529, 2922.416475141421, 2692.8409902649, 2447.8570429496467, 2036.1059530228376, 1866.5779191209003, 1809.868225327693, 1558.2388484934345, 1330.1253253752366, 1254.3396296054125, 1104.8163758879527, 865.0637027751654, 799.2259566113353, 810.760053393431, 724.6756397057325, 611.2978768311441, 544.1322553129867, 484.32290973793715, 427.09226111136377, 363.7518107574433, 348.393721383065, 295.7854479914531, 273.7027771472931, 245.82586858421564, 222.1324712857604, 199.30874756351113, 188.72418217454106]),
                "fail": np.array([45511.04403647699, 41079.93012525796, 40727.236316832365, 34559.57936263963, 31946.14810642379, 28878.48818682891, 25820.406854955247, 23166.404159538797, 21632.704532302916, 18499.283225926803, 16938.807095905882, 15601.445229024393, 13386.164538178826, 12101.94353529031, 11723.452738028369, 9968.97439501679, 8947.14102957386, 7981.194391448749, 7301.207673351048, 6364.0393990497105, 5687.775319714332, 5151.827819460072, 4674.958705232013, 3984.1076764897443, 3728.472815330606, 3182.867656349903, 2830.3143796238583, 2567.943558348343, 2264.352760080714, 2026.9894705256447, 1752.9333959259093, 1568.080979383085, 1460.6756752429064, 1305.7658144892193, 1126.8351017078385, 985.7203802638687, 880.427430575015, 799.109357524896, 746.6917918575928, 657.3273142203689, 549.4357224991545, 517.6786587073002, 439.8946387998294, 393.87680231104605, 340.0927645144984, 312.42521196906455, 276.5220663778018]),
            }
        },
        "errs": {
            "bkg": {
                "pass": np.array([631.2181869262554, 655.3758529791091, 714.9869986512043, 656.0741138556303, 470.6310343954024, 443.53216728071857, 418.14581996309556, 403.4256742144957, 611.4472373112185, 416.50324158329187, 407.12306678579563, 249.2919698140988, 241.84160992128878, 332.1914872956116, 370.5177986259088, 272.1201523464285, 98.65530740347248, 95.52598345660434, 87.71910472290624, 82.29656320295523, 79.96778582590201, 67.45872197791849, 66.38003996801268, 61.298545662483974, 53.513718464011816, 50.89395650703222, 50.42907588910356, 43.98447878249719, 39.19306572081908, 41.395170835186676, 35.27519355597898, 24.847341552686213, 24.275459648772895, 29.754391060999662, 28.429830002600543, 23.616795802231117, 20.80947479489803, 19.240678632959074, 18.777020827685543, 15.417700513572225, 17.093567724910024, 13.32405839813672, 14.244323650828605, 10.992186582122415, 9.85599957773463, 10.141288632267049, 11.503355948560422]),
                "fail": np.array([760.8567947606394, 725.2335724844102, 924.7935764792179, 615.3961553050276, 651.647096708545, 627.7979761059515, 548.7652264416174, 529.1127690388306, 544.9601627629468, 449.411747729661, 346.5823467995972, 413.69490761356076, 237.13327597049351, 266.50569509295866, 366.3737398320233, 213.4492822619158, 284.3828445671478, 211.0307284748164, 158.51657895508623, 102.98171996532403, 95.28305459795044, 88.33895564078259, 82.69350973396938, 70.89950581870676, 70.72349674157685, 61.66976308413847, 54.257129945943895, 53.14789531717416, 47.37930523453283, 43.21735002661446, 37.11752082685987, 31.82212342437327, 32.0847603536062, 32.57502514593756, 24.295399399432448, 23.453917154281786, 21.8665497521834, 19.715356365692028, 22.03313264407499, 20.28984755550465, 14.455323142397619, 14.625956340233568, 9.570181657352864, 7.93838311909737, 5.379828626575783, 9.233103638155825, 7.16338492393075]),
            }
        },
    }

    # different initializations
    def make_init(npar, basis):
        inits = None
        if basis=='Bernstein':
            inits = np.ones(npar+1)
        elif basis=='Chebyshev':
            inits = np.zeros(npar+1)
            inits[0] = 1
        return inits

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
                    vals=template_info["vals"]["bkg"][region],
                ),
            }
            ch.setObservation(templates["bkg"])
            ch._observation = (ch._observation, template_info["errs"]["bkg"][region]**2)

        bkgeff = template_info["yields"]["bkg"]["pass"] / template_info["yields"]["bkg"]["fail"]
        npar_mc = 2
        tf_mc = rl.BasisPoly("tf_mc", (npar_mc,), ["mt"], basis=basis, init_params=make_init(npar_mc,basis))
        tf_mc_params = bkgeff * tf_mc(mtscaled)
        failCh = bkgmodel["fail"]
        passCh = bkgmodel["pass"]
        bkgparams = np.array([rl.IndependentParameter("bkgparam_mtbin%d" % (i), 0) for i in range(mt.nbins)])
        initial_bkg = get_obs(failCh).astype(float)  # was integer, and numpy complained about subtracting float from it
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
            ROOT.RooFit.SumW2Error(False),
            ROOT.RooFit.Strategy(0),
            ROOT.RooFit.Save(),
            ROOT.RooFit.Minimizer("Minuit2", "migrad"),
            ROOT.RooFit.PrintLevel(-1),
        )
        bkgfit_ws.add(bkgfit)
        if "pytest" not in sys.modules:
            bkgfit_ws.writeToFile(os.path.join(str(tmpdir), "svjModel_bkgfit.root"))
        if bkgfit.status() != 0:
#            raise RuntimeError("Could not fit bkg")
            print(f"Could not fit bkg: status {bkgfit.status()}")

        bkgmodel.readRooFitResult(bkgfit)

        # save MC TF details for later use in plotting
        paramfile = os.path.join(str(tmpdir), "svjModel_mctf")
        np.save(paramfile, [par.value for par in tf_mc.parameters.flatten()])
        param_names = [p.name for p in tf_mc.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_mc.name + "_deco", bkgfit, param_names)
        decofile = os.path.join(str(tmpdir), "svjModel_deco")
        np.save(decofile, decoVector._transform)
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
                vals=template_info["vals"]["bkg"][region],
            ),
        }
        yields = sum(tpl[0] for tpl in templates.values())
        if throwPoisson:
            yields = np.random.poisson(yields)
        data_obs = (yields, mt.binning, mt.name)
        ch.setObservation(data_obs)

    npar_data = 2
    tf_data = rl.BasisPoly("tf_data", (npar_data,), ["mt"], basis=basis, init_params=make_init(npar_data,basis))
    tf_data_params = tf_data(mtscaled)
    if tfFromMC:
        tf_params = bkgeff * tf_mc_params_final * tf_data_params
    else:
        tf_params = bkgeff * tf_data_params

    failCh = model["fail"]
    passCh = model["pass"]

    bkgparams = np.array([rl.IndependentParameter("bkgparam_mtbin%d" % (i), 0) for i in range(mt.nbins)])
    initial_bkg = get_obs(failCh).astype(float)  # was integer, and numpy complained about subtracting float from it
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
