from __future__ import print_function, division
import sys
import os
from contextlib import contextmanager
import rhalphalib as rl
import numpy as np
import scipy.stats
import pickle
import ROOT
import argparse
import subprocess
import shlex

rl.util.install_roofit_helpers()

basis = "Bernstein"

def predef_sample(obs, data):
    return (data["vals"], obs.binning, obs.name)


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


def arg_name(args):
    return f"sig-{args.sig}_bkg-{args.bkg}_obs-{args.obs}_nmc-{args.nmc}_ndata-{args.ndata}"


def plot_tf(mtbins, fail_bkg, pass_bkg, fitresult, args):
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

    outfile = f"{args.dir}/tf_mc_bsvj.png"
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


def plot_sr_cr(fail_bkg, pass_bkg, args):
    outfile = f"{args.dir}/srcr_mc_bsvj.png"
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

def run_and_log(cmd, dir, fname):
    print(f"Running: {cmd}")
    with open(os.path.join(dir,fname), 'w') as logfile:
        subprocess.run(shlex.split(cmd), cwd=dir, stdout=logfile, stderr=logfile, text=True)


def get_bias(fname, test):
    # undo change from install_roofit_helpers
    ROOT.TH1.AddDirectory(True)

    file = ROOT.TFile.Open(fname)
    tree = file.Get("tree_fit_sb")
    hist = ROOT.TH1F(test, test, 50, -5, 5)
    n = tree.Draw(f"(r-0)/(0.5*(rLoErr+rHiErr))>>{test}", "fit_status==0 || fit_status==1", "goff")
    return hist.GetMean(), hist.GetRMS()


def test_rhalphabet(args, sig_data, bkg_data, obs_data):
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

    msig = 250

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
            template = predef_sample(mt, bkg_data[region])
            ch.setObservation(template)
            ch._observation = (ch._observation, bkg_data[region]["errs"]**2)

        bkgeff = np.sum(bkg_data["pass"]["vals"]) / np.sum(bkg_data["fail"]["vals"])
        npar_mc = args.nmc-1
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
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Strategy(0),
            ROOT.RooFit.Save(),
            ROOT.RooFit.Minimizer("Minuit2", "migrad"),
            ROOT.RooFit.PrintLevel(10 if args.verbose else -1),
        )
        bkgfit_ws.add(bkgfit)
        if "pytest" not in sys.modules:
            bkgfit_ws.writeToFile(os.path.join(str(args.dir), "svjModel_bkgfit.root"))
        if bkgfit.status() != 0:
#            raise RuntimeError("Could not fit bkg")
            print(f"Could not fit bkg: status {bkgfit.status()}")

        bkgmodel.readRooFitResult(bkgfit)

        # save MC TF details for later use in plotting
        paramfile = os.path.join(str(args.dir), "svjModel_mctf")
        np.save(paramfile, [par.value for par in tf_mc.parameters.flatten()])
        param_names = [p.name for p in tf_mc.parameters.reshape(-1)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_mc.name + "_deco", bkgfit, param_names)
        #print(decoVector.correlated_str)
        decofile = os.path.join(str(args.dir), "svjModel_deco")
        np.save(decofile, decoVector._transform)
        tf_mc.parameters = decoVector.correlated_params.reshape(tf_mc.parameters.shape)
        tf_mc_params_final = tf_mc(mtscaled)

        fail_th1 = obs_to_TH1(bkgmodel["fail"])
        pass_th1 = obs_to_TH1(bkgmodel["pass"])
        plot_sr_cr(fail_th1, pass_th1, args)
        plot_tf(mtbins, fail_th1, pass_th1, bkgfit, args)

    # build actual fit model now
    model = rl.Model("svjModel")

    for region in ["pass", "fail"]:
        ch = rl.Channel(region)
        model.addChannel(ch)

        template_sig = predef_sample(mt, sig_data[region])
        # some mock expectations
        sample = rl.TemplateSample(ch.name + "_sig", rl.Sample.SIGNAL, template_sig)

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
        template_obs = predef_sample(mt, obs_data[region])
        yields = template_obs[0]
        if throwPoisson:
            yields = np.random.poisson(yields)
        data_obs = (yields, mt.binning, mt.name)
        ch.setObservation(data_obs)

    npar_data = args.ndata-1
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

    with open(os.path.join(str(args.dir), "svjModel.pkl"), "wb") as fout:
        pickle.dump(model, fout)

    svjModelPath = os.path.join(str(args.dir), "svjModel")
    model.renderCombine(svjModelPath)
    if args.mode!="dryrun":
        subprocess.run(["bash", "build.sh"], cwd=svjModelPath)

    # build simple model for bias tests
    model2 = rl.Model("simpleModel")
    for region in ["pass", "fail"]:
        # copy sig and obs from full model, use bkg histograms directly
        ch = rl.Channel(region)
        model2.addChannel(ch)
        ch.addSample(model[region][ch.name+"_sig"])
        template_bkg = predef_sample(mt, bkg_data[region])
        bkg_sample = rl.TemplateSample(ch.name+"_bkg", rl.Sample.BACKGROUND, template_bkg)
        ch.addSample(bkg_sample)
        ch._observable = model[region]._observable
        ch._observation = model[region]._observation

    with open(os.path.join(str(args.dir), "simpleModel.pkl"), "wb") as fout:
        pickle.dump(model2, fout)

    simpleModelPath = os.path.join(str(args.dir), "simpleModel")
    model2.renderCombine(simpleModelPath)
    if args.mode!="dryrun":
        subprocess.run(["bash", "build.sh"], cwd=simpleModelPath)

    if args.mode=="bias":
        toy_cmd = "combine -M GenerateOnly -d model_combined.txt --saveWorkspace --toysFrequentist --bypassFrequentistFit --saveToys --name _{model} -t 300 -s 995 -v 0 --expectSignal 0.0 --setParameterRanges r=-3.0,5.0"
        fit_cmd = "combine -M FitDiagnostics -d model_combined.txt --noPreFitValue --savePredictionsPerToy --toysFrequentist --saveToys --name _{test} -t 300 -s 995 -v 0 --expectSignal 0.0 --X-rtd MINIMIZER_MaxCalls=100000 --setParameterRanges r=-5.0,5.0 --cminDefaultMinimizerStrategy 0 --toysFile {toyfile}"
        # self test
        # toys generated from rhalphabet background model
        # fit to rhalphabet background model
        toy_cmd_self = toy_cmd.format(model="svjModel")
        fit_cmd_self = fit_cmd.format(test="self", toyfile="higgsCombine_svjModel.GenerateOnly.mH120.995.root")
        run_and_log(toy_cmd_self, svjModelPath, "log_gen_self.log")
        run_and_log(fit_cmd_self, svjModelPath, "log_fit_self.log")

        # bias test
        # toys generated from MC background directly
        # fit to rhalphabet background model
        toy_cmd_bias = toy_cmd.format(model="simpleModel")
        fit_cmd_bias = fit_cmd.format(test="bias", toyfile="../simpleModel/higgsCombine_simpleModel.GenerateOnly.mH120.995.root")
        run_and_log(toy_cmd_bias, simpleModelPath, "log_gen_bias.log")
        run_and_log(fit_cmd_bias, svjModelPath, "log_fit_bias.log")

        # print results after Combine spew
        result_self = get_bias(os.path.join(svjModelPath, "higgsCombine_self.FitDiagnostics.mH120.995.root"), "self")
        print(f"self: mean = {result_self[0]:.4f}, stdev = {result_self[1]:.4f}")
        result_bias = get_bias(os.path.join(svjModelPath, "higgsCombine_bias.FitDiagnostics.mH120.995.root"), "bias")
        print(f"bias: mean = {result_bias[0]:.4f}, stdev = {result_bias[1]:.4f}")


if __name__ == "__main__":

    bkg = {"MC": {}, "smooth": {}, "fake": {}}
    # bias_study_20260807/regular/merged_20260807/bkg_sel-antiloosertcutbased_ddt=0.1_mt.json
    bkg["MC"]["fail"] = {
        "vals": np.array([288141.38752316753,269407.2147902497,252804.19023079262,228657.95430641156,207807.24484041147,189163.40492648818,171856.5262428308,155678.14664101973,143247.43561054114,129054.5512668509,115582.62581490912,104652.99883200787,93970.61376103945,86686.57287738658,77620.86720058508,70401.98366039619,63997.7711375393,56727.678735362366,51980.663534175605,46085.701999424025,41613.9747704491,36362.91602727026,33012.11931016669,29537.487635223195,27187.24608475156,23085.75159085542,20146.485513418913,18266.643162924796,16053.329891137779,14176.836655143648,12165.380147825927,10263.944158747792,9242.990741532296,8119.557855028659,6646.001165270805,6063.807853285223,5234.572751887143,4690.165772747248,3918.750500243157,3532.448755275458,2919.263520631939,2579.059108994901,2259.828122854233,1955.999055672437,1746.689328867942,1564.3625855781138,1349.5834227576852]),
        "errs": np.array([1785.8842585186796,1750.4743222823008,1803.8056712474242,1572.8702757858546,1460.5531168788195,1419.5922447802168,1322.1779910512064,1298.8283758571965,1249.0801548559305,1229.7054923597752,992.4759033393677,924.843058706853,857.1960045931104,909.3037251059645,799.2760195370241,717.5973359160615,731.369398405049,630.4123099454363,685.1233505619404,582.6654105554143,601.3579141913402,400.55635515701834,448.65631318095194,409.99449456731423,596.5985260477077,407.896412847357,295.6416948806689,257.20061183013615,240.4491680635607,310.56554602017,234.69599038724894,204.20440283467306,198.84035522535567,121.76292368125605,102.39189887913903,99.1422652449659,89.9930217876848,83.9716365970292,72.91890292612705,70.62662280836528,57.30201815381456,54.04455232002697,46.90556860645281,41.63977607152207,37.50013593908469,37.715413376740926,32.94947597965393]),
    }
    # bias_study_20260807/regular/smooth_20260807/bkg_sel-antiloosertcutbased_ddt\=0.1_mt_wide_smooth.json
    bkg["smooth"]["fail"] = {
        "vals": np.array([278533.5915702801,260835.41368958273,242812.12776779843,224756.59332237436,206992.9373437459,189783.5488390035,173339.11692435908,157806.16037825486,143350.84839372264,130121.09187367672,117990.40121171028,106822.61604658907,96553.57711713735,87250.70391867163,78921.82812464926,71366.29512886133,64394.138526139206,57961.25441656184,52109.7155533606,46828.87965910406,42062.56625394551,37752.5943228945,33791.337940645055,30132.156791880378,26786.290867916952,23755.574473867164,21028.768035285124,18578.393283335492,16383.10369431029,14411.030256753487,12640.88227773623,11056.665465480757,9643.018295732918,8388.821668587032,7291.66973614883,6341.129118250843,5514.1124613713855,4797.251653945878,4172.296032688438,3628.9239729054716,3150.2446533437587,2727.33526208324,2354.506871402281,2028.063096027727,1745.1680376016793,1503.0668915494891,1298.2743670115033]),
        "errs": np.array([500.03481927210055,501.98199552067746,502.8857250002012,500.8027719601805,499.8633946561505,500.9037485486765,503.9183251185113,506.20012550464753,505.8469680835726,502.6259140048519,504.08194345647604,507.81598226362837,511.1547516204213,506.87921980659917,497.65341325952573,499.69582427437416,512.1663058623436,524.3018846017985,525.729706874233,520.1410987488307,516.7474567328184,517.7334113257367,526.0124789367786,531.5913871225704,529.4468253436515,526.1291451389746,531.0899846057671,544.5520757754875,555.5276066671088,560.6105060656739,557.9149734703898,549.7752692849721,540.5330747834541,532.5401084985299,522.4149569986264,512.6382197555492,509.5544138558283,508.26057431084206,511.23805414798414,516.154969974688,533.5679192687056,563.717981100247,599.8361886559888,633.916620807778,659.2602675674167,671.2463284479219,667.0112818438622]),
    }
    bkg["fake"]["fail"] = bkg["smooth"]["fail"]
    # bias_study_20260807/regular/merged_20260807/bkg_sel-rtcutbased_ddt\=0.1_mt.json
    bkg["MC"]["pass"] = {
        "vals": np.array([25553.111139234214,23687.431316580565,22394.075393573148,19678.231506056152,18464.726244880003,16175.779406012618,14744.273975052754,13577.069436682621,12390.181059993338,10823.352646310814,10075.220489688916,8995.956856843783,8175.03228028724,7373.295597264776,6591.50761317322,6074.291638907511,5570.651509478223,5071.291944093304,4572.885039018001,3937.6858026357368,3682.9572260770947,3279.7719049989246,2938.291728928685,2622.1095529422164,2297.6027755350806,1963.6125096897595,1899.9059703480452,1625.91901233932,1442.4441337427124,1236.9617746286094,1091.3808683226816,883.5099325496703,771.0283229243942,746.955849446822,680.254356124904,509.7412402983755,461.78124475991353,392.4589970605448,360.95183352800086,315.2484003910795,278.59602298680693,221.43371677957475,189.82085445942357,178.85419704625383,150.79990029893816,129.83915504766628,116.63618205301464]),
        "errs": np.array([248.97917164577217,286.8291190837316,383.37883508799274,178.3559429264958,240.7308330635058,187.13969564861387,186.9052851966286,255.10286714071108,319.19935291632146,179.64762700572615,248.75222230196755,120.06999618806404,113.3113907549284,107.853520322161,98.52754152132454,95.5738007890168,224.17432298288907,195.9456378644362,82.02561678262938,71.90615032361518,68.74881528928142,64.78389989542175,62.121948110119526,55.29595256835904,52.759542374678546,44.001001367087305,47.40956521215322,42.2530651353153,37.00797362582923,33.62037998810922,32.08136120206618,23.913101741536618,20.47308983667378,23.761541271646394,24.501591974631893,11.222154605190706,12.38888340904769,9.9825602771704,11.574380719198992,10.254682290127256,11.782538983705583,3.9613561284718606,3.5698982481918358,3.6457147375497545,2.8625580810684905,2.7064597997493833,2.586115504387575]),
    }
    # bias_study_20260807/regular/smooth_20260807/bkg_sel-rtcutbased_ddt\=0.1_mt_wide_smooth.json
    bkg["smooth"]["pass"] = {
        "vals": np.array([24588.539063977292,22943.48801011798,21316.919217062677,19706.197563608574,18113.24580589489,16554.587711790944,15035.754195225394,13601.918354837706,12324.253830531074,11168.709215113151,10107.723201824718,9138.124129524378,8261.148731480242,7475.004874816148,6773.554546975516,6139.563402449903,5561.016509350468,5029.9436510073765,4540.76364998223,4091.4952204205133,3682.6045620283794,3307.6555984040924,2959.550414382849,2639.7986022106734,2351.0856744817083,2087.313177047756,1846.9303363925128,1628.472038609442,1431.3588636760776,1259.0497453647044,1110.3619113342352,972.6100326056346,847.436654981838,737.2840422175971,641.3320779271967,557.6210800155154,484.4690099991519,420.5024947831607,364.6164553820584,315.851465381391,273.48313843393794,236.50973378832288,204.30062117809933,176.38606018538943,152.37318613134778,131.88489708293,114.54433959147093]),
        "errs": np.array([50.90633358647332,51.32034640716677,51.684246480950065,51.40942780209957,51.13111057809641,51.858462075287,52.9184351564796,52.47052492570073,52.47453754762201,53.68856544254438,54.40859517212769,54.255783807939835,53.50979964514558,52.70330217110621,52.28004161635356,52.51185467780744,52.643849649043766,52.23476197272861,51.654777128203726,50.98064814389369,49.96927482670214,49.344912656560645,49.29513652465025,48.43745964593918,47.55384369162157,48.875761030448174,49.969665517922614,50.939509712880295,52.02079174369495,51.91273482699839,52.88639796554576,54.32536982520547,53.96216986543881,52.943015190711094,54.08629685836392,56.66076083901921,58.425359202812906,58.89056099963148,58.51208943675633,57.70176817482859,56.51046533400582,55.94510372011074,56.37800321992705,57.817737033910376,59.911485429734476,62.000998307906876,63.2373597361721]),
    }
    # bias_study_20260807/scaled_bkg/smooth_20260807/bkg_sel-rtcutbased_ddt\=0.1_mt_wide_smooth.json
    # = bkg["smooth"]["fail"] * 0.087
    bkg["fake"]["pass"] = {
        "vals": np.array([24241.053934615196,22700.764010761803,21132.179612571817,19560.788593476766,18014.800046004817,16517.04994493686,15085.874772421199,13734.026202647745,12475.96610482174,11324.567311102672,10268.821305951033,9296.877918470638,8403.153304710035,7593.515050022207,6868.644752693407,6211.079243904691,5604.285559543912,5044.425293583687,4535.160079346261,4075.5637089019297,3660.7467395290473,3285.64561995395,2940.893559458875,2622.431405276118,2331.23738496592,2067.471139916278,1830.154478847609,1616.8959408480998,1425.8377168484585,1254.206215264298,1100.1484860342484,962.2725301270805,839.24141889493,730.0874460765939,634.6012283505321,551.8747383225113,479.89866078293176,417.50955576172,363.11904998381993,315.82884224441244,274.1689076653815,237.36268510096681,204.9150615453983,176.504336932754,151.8837002243389,130.8133980536132,112.99010211022716]),
        "errs": np.array([50.90633358647332,51.32034640716677,51.684246480950065,51.40942780209957,51.13111057809641,51.858462075287,52.9184351564796,52.47052492570073,52.47453754762201,53.68856544254438,54.40859517212769,54.255783807939835,53.50979964514558,52.70330217110621,52.28004161635356,52.51185467780744,52.643849649043766,52.23476197272861,51.654777128203726,50.98064814389369,49.96927482670214,49.344912656560645,49.29513652465025,48.43745964593918,47.55384369162157,48.875761030448174,49.969665517922614,50.939509712880295,52.02079174369495,51.91273482699839,52.88639796554576,54.32536982520547,53.96216986543881,52.943015190711094,54.08629685836392,56.66076083901921,58.425359202812906,58.89056099963148,58.51208943675633,57.70176817482859,56.51046533400582,55.94510372011074,56.37800321992705,57.817737033910376,59.911485429734476,62.000998307906876,63.2373597361721]),
    }

    # m = 250
    sig = {"smooth": {}, "fakeWidth20": {}}
    # bias_study_20260807/regular/smooth_20260807/SVJ_s-channel_mMed-250_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-antiloosertcutbased_ddt\=0.1_mt_smooth.json
    sig["smooth"]["fail"] = {
        "vals": np.array([791.5806142852589,876.2432205269495,951.9978721449696,1014.6449677229543,1061.2299110305476,1089.2317995851201,1090.1937586064666,1060.7999428305232,1007.7575861649733,942.6491538803242,875.4904137338225,807.79337756228,739.9525685542284,674.6930536984694,613.3255874379154,556.2361996890331,506.24674707215246,464.29973616035943,426.50448970669436,390.94103095451055,357.8849747887239,329.71517304407433,307.54956991054985,289.174944596326,271.95644494005535,254.91546870900208,237.08599953322798,218.80089833932684,201.7686558313289,186.5321134506327,172.72477170710275,160.1539985906047,149.66228297314154,140.93887151514838,132.72068399802168,124.95033362329792,117.88590301915008,112.53190134225801,107.68569490179011,101.41522760464647,93.31036214076852,82.94098418661504,70.77165410086165,58.79971436399397,48.151167635631396,39.220425731867806,32.159937468316215]),
        "errs": np.array([9.20690229813298,9.047963341189218,8.758305532623424,8.437215309874977,8.188195839135329,8.101909971567343,8.116359163519512,8.138236022997273,8.171380543095198,8.231270981597481,8.318993944586783,8.355170676085118,8.314715118943782,8.325603240271635,8.401773703794582,8.46081800960943,8.454429293777315,8.478664934007101,8.573642840613612,8.613518339353002,8.555140047006095,8.453927265206206,8.416128992857844,8.473272389228265,8.529772567975662,8.499202493091754,8.451380155467941,8.451487219142875,8.502411603842702,8.603000634057038,8.647390032683028,8.61432052709624,8.540224528192804,8.501769018108394,8.482032182868867,8.462103846437572,8.484937621964223,8.312134041327344,8.306984409141009,8.472884605501763,8.644176055529044,8.743440696330728,8.96511059612392,9.695609855695748,10.700100079089664,11.56295220056244,12.06266602557268]),
    }
    # bias_study_20260807/scaled_bkg_sigWidth20/smooth_20260807/SVJ_s-channel_mMed-250_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-antiloosertcutbased_ddt\=0.1_mt_smooth.json
    sig["fakeWidth20"]["fail"] = {
        "vals": np.array([20.642873980704852,92.51494269969096,322.90887871185214,877.7573372505136,1858.212297541234,3063.674140432788,3933.8354647646697,3933.8354647646697,3063.674140432788,1858.212297541234,877.7573372505136,322.90887871185214,92.51494269969096,20.642873980704852,3.5871936157776707,0.48547386401583986,0.05116856881232149,0.004200171900546281,0.00026850800630243577,1.3368226467098212e-05,5.18342391376793e-07,1.5652583936446494e-08,3.6811349100719933e-10,6.74223377135916e-12,9.617279958382127e-14,1.0683832976497416e-15,9.243326651370472e-18,6.22810450721684e-20,3.2682092993644875e-22,1.3356424419565449e-24,4.2510571152625686e-27,1.0537317077590104e-29,2.034180733717669e-32,3.058267040671018e-35,3.5808625528590293e-38,3.265323983081654e-41,2.3189494640579317e-44,1.2825747037331374e-47,5.524594863641213e-51,1.853295111053801e-54,4.841892069681497e-58,9.851718823436665e-62,1.5611164194295638e-65,1.9265707147830112e-69,1.8516589421089724e-73,1.386000601447804e-77,8.079643899218167e-82]),
        "errs": np.array([9.20690229813298,9.047963341189218,8.758305532623424,8.437215309874977,8.188195839135329,8.101909971567343,8.116359163519512,8.138236022997273,8.171380543095198,8.231270981597481,8.318993944586783,8.355170676085118,8.314715118943782,8.325603240271635,8.401773703794582,8.46081800960943,8.454429293777315,8.478664934007101,8.573642840613612,8.613518339353002,8.555140047006095,8.453927265206206,8.416128992857844,8.473272389228265,8.529772567975662,8.499202493091754,8.451380155467941,8.451487219142875,8.502411603842702,8.603000634057038,8.647390032683028,8.61432052709624,8.540224528192804,8.501769018108394,8.482032182868867,8.462103846437572,8.484937621964223,8.312134041327344,8.306984409141009,8.472884605501763,8.644176055529044,8.743440696330728,8.96511059612392,9.695609855695748,10.700100079089664,11.56295220056244,12.06266602557268]),
    }
    # bias_study_20260807/regular/smooth_20260807/SVJ_s-channel_mMed-250_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-rtcutbased_ddt\=0.1_mt_smooth.json
    sig["smooth"]["pass"] = {
        "vals": np.array([1559.939446631649,1716.0007565000124,1843.8621887540342,1933.430949057905,1979.0872500127582,1975.0039105829037,1913.1649116303797,1796.9885220222227,1635.0801371346815,1446.0251479031178,1253.6152083162315,1072.4310646845977,910.678438524635,770.8099832617697,651.4156809576932,551.1708749110683,468.2791557587155,401.1311166840083,350.0719940325394,311.62283680397496,280.1411858134862,253.00740646110944,229.6181249880338,208.97350382323597,189.70363414597588,172.03221289733293,154.87515596247533,137.7435176818577,120.75497158422665,103.94591519359994,87.23241743969307,71.71400462745565,59.649928408961216,50.70069126383767,43.29786837937861,36.81821838777725,31.139965188808915,26.44159702293264,23.009285545015196,20.17533328736744,17.759356144508605,15.718937500997605,14.083855358493215,12.668021106102264,11.350968156648328,10.214791633184369,9.417607827302515]),
        "errs": np.array([14.289221782834932,14.0186691449662,13.570496413091785,13.101161299453269,12.740394591328187,12.587792335396786,12.648144956514077,12.910650821367934,13.09548141526676,13.226095237639163,13.475568146034778,13.657621369517214,13.693645473970344,13.684319651182534,13.72163051349427,13.712020842578147,13.632467121484531,13.561777615728056,13.469700808962818,13.502235446470827,13.617234183451826,13.588620901895338,13.304983513157199,13.216068903921748,13.36059699839922,13.425741732478244,13.40404456261155,13.534264646090849,14.121918755190842,14.618597532040289,14.853949690350067,14.876740924596655,15.046364956364435,15.464883958928633,15.681739507249228,15.725914429453203,15.64748161701256,14.966990215792046,14.49017700908483,15.100053188313845,15.649687070403177,15.607930195525507,15.381792926679962,15.96514604029225,16.88224996828362,17.819433496467912,19.0249140129509]),
    }
    # bias_study_20260807/scaled_bkg_sigWidth20/smooth_20260807/SVJ_s-channel_mMed-250_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-rtcutbased_ddt\=0.1_mt_smooth.json
    sig["fakeWidth20"]["pass"] = {
        "vals": np.array([27.328703942578233,122.47875376595765,427.49285564571323,1162.0460612978084,2460.0515310721526,4055.9392862970763,5207.929132147784,5207.929132147784,4055.9392862970763,2460.0515310721526,1162.0460612978084,427.49285564571323,122.47875376595765,27.328703942578233,4.749016653491547,0.642709513895669,0.06774108438313628,0.0055605268183711535,0.00035547258668100636,1.7697937975988706e-05,6.862235252739118e-07,2.072215490611138e-08,4.873383726707798e-10,8.925913650461365e-12,1.2732132016171032e-13,1.4144121049208803e-15,1.223706242337974e-17,8.245267803323306e-20,4.326719450411692e-22,1.7682313472189446e-24,5.627892775714875e-27,1.3950151467848122e-29,2.693012760211045e-32,4.048780930889998e-35,4.74063508102702e-38,4.322899635663612e-41,3.0700126067850834e-44,1.6979759889696363e-47,7.313904913252968e-51,2.4535417624289342e-54,6.41008781131574e-58,1.3042501121834232e-61,2.0667320105895254e-65,2.5505499252652257e-69,2.451375669823076e-73,1.8348995462844495e-77,1.0696485203057576e-81]),
        "errs": np.array([14.289221782834932,14.0186691449662,13.570496413091785,13.101161299453269,12.740394591328187,12.587792335396786,12.648144956514077,12.910650821367934,13.09548141526676,13.226095237639163,13.475568146034778,13.657621369517214,13.693645473970344,13.684319651182534,13.72163051349427,13.712020842578147,13.632467121484531,13.561777615728056,13.469700808962818,13.502235446470827,13.617234183451826,13.588620901895338,13.304983513157199,13.216068903921748,13.36059699839922,13.425741732478244,13.40404456261155,13.534264646090849,14.121918755190842,14.618597532040289,14.853949690350067,14.876740924596655,15.046364956364435,15.464883958928633,15.681739507249228,15.725914429453203,15.64748161701256,14.966990215792046,14.49017700908483,15.100053188313845,15.649687070403177,15.607930195525507,15.381792926679962,15.96514604029225,16.88224996828362,17.819433496467912,19.0249140129509]),
    }

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--sig", type=str, required=True, choices=list(sig.keys()), help="signal type")
    parser.add_argument("--bkg", type=str, required=True, choices=list(bkg.keys()), help="background type")
    parser.add_argument("--obs", type=str, default=None, choices=list(bkg.keys()), help="observation type (None -> reuse background)")
    parser.add_argument("--nmc", type=int, default=1, help="number of params for MC TF")
    parser.add_argument("--ndata", type=int, default=1, help="number of params for data TF")
    parser.add_argument("--verbose", default=False, action="store_true", help="verbose fit printouts")
    parser.add_argument("--mode", type=str, default="dryrun", choices=["dryrun","build","bias"], help="mode of operation (determines what combine commands run after creating model & cards)")
    args = parser.parse_args()

    if args.obs is None: args.obs = args.bkg

    args.dir = "tmp_"+arg_name(args)
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)
    test_rhalphabet(args, sig[args.sig], bkg[args.bkg], bkg[args.obs])
