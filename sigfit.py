import json, os
import numpy as np
from time import strftime
from contextlib import contextmanager
from boosted_fits import Decoder, cut_histograms, get_mt, th1_to_hist, do_fisher_test
import ROOT

def get_f_chi2(f_th1, f_vals):
    f_np = th1_to_hist(f_th1)
    chi2 = np.sum((f_np['vals']-f_vals)**2/(f_np['errs']**2))
    return(chi2)

mt_min = 180
mt_max = 650
sigfile = "smooth_20241101/SVJ_s-channel_mMed-350_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.json"
with open(sigfile, 'r') as f:
    sig = json.load(f, cls=Decoder)
    sig = cut_histograms(sig,mt_min,mt_max)

mt = sig['central'].binning
n_bins = len(mt)-1
sig_name = 'sig'
mtname = 'mt'
rebin = 1

mtvar = get_mt(mt[0], mt[-1], int(n_bins/rebin), name=mtname)
sig_th1 = sig['central'].th1(sig_name, rebin, mtname)
sig_datahist = ROOT.RooDataHist(sig_name, sig_name, ROOT.RooArgList(mtvar), sig_th1, 1.)

import rhalphalib as rl
rl.util.install_roofit_helpers()

mt_array = np.array(mt)
mtpts = mt_array[:-1] + 0.5*np.diff(mt_array)
mtscaled = (mtpts - min(mtpts))/(max(mtpts) - min(mtpts))
rl_mt = rl.Observable(mtname, mt_array)

npar_vals = range(20)
basis_mc = 'Bernstein'
sigmodels = []
sigyield = sig_th1.Integral()
for npar_mc in npar_vals:
    # use rhalphalib model to create workspace
    sigmodel = rl.Model("sigmodel")
    ch = rl.Channel("default")
    sigmodel.addChannel(ch)
    ch.setObservation(sig_th1, read_sumw2=True)

    # signal model from polynomial function
    f_mc = rl.BasisPoly("f_mc", (npar_mc,), ["mt"], basis=basis_mc)
    f_mc_params = sigyield * f_mc(mtscaled)
    ch_sig = rl.ParametericSample("default_sig", rl.Sample.SIGNAL, rl_mt, f_mc_params)
    ch.addSample(ch_sig)

    # do the fit
    sigfit_ws = ROOT.RooWorkspace("sigfit_ws")
    simpdf, obs = sigmodel.renderRoofit(sigfit_ws)
    sigfit = simpdf.fitTo(
        obs,
        ROOT.RooFit.Extended(True),
        ROOT.RooFit.SumW2Error(True),
        ROOT.RooFit.Strategy(0),
        ROOT.RooFit.Save(),
        ROOT.RooFit.Minimizer("Minuit2", "migrad"),
        ROOT.RooFit.PrintLevel(-1),
    )
    sigfit_ws.add(sigfit)

    # actual goodness of fit: sigfit.minNll()
    # but this is not a saturated gof, so can't compare nll from different fits
    # therefore, just use chi2 instead (approximation)

    # treat fit result as unweighted histo
    f_mc.update_from_roofit(sigfit)
    f_mc_vals = f_mc(mtscaled, nominal=True)
    chi2 = get_f_chi2(sig_th1, sigyield * f_mc_vals)

    # save for later
    sigmodels.append({
        "model": sigmodel,
        "f": f_mc,
        "ws": sigfit_ws,
        "fit": sigfit,
        "results": (len(f_mc.parameters), chi2),
    })

results = [sm["results"] for sm in sigmodels]
i_winner = do_fisher_test(results, n_bins)
print(f'chose n_pars={results[i_winner][0]} for f_mc')

outdir = strftime(f'dc_%Y%m%d_{sig["central"].metadata["selection"]}')
if not os.path.isdir(outdir): os.makedirs(outdir)

'''
# save winner
wsfile = f'{outdir}/dc_{os.path.basename(sigfile).replace(".json","")}.root'
sigmodel = sigmodels[i_winner]['model']
f_mc = sigmodels[i_winner]['f']
sigfit = sigmodels[i_winner]['fit']
sigfit_ws = sigmodels[i_winner]['ws']
print(f'MC fit status: {sigfit.status()}')
sigfitfile = wsfile.replace("/dc_", "/sigfit_")
sigfitf = ROOT.TFile.Open(sigfitfile, "RECREATE")
sigfitf.cd()
sigfit.Write("sigfit")
sigfit_ws.Write()
sigfitf.Close()
'''

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

def plot_f(mtbins, sig_th1, fitresults, outdir):
    mtpts = mtbins[:-1] + 0.5 * np.diff(mtbins)
    mtscaled = (mtpts - min(mtbins)) / (max(mtbins) - min(mtbins))

    sigyield = sig_th1.Integral()
    sig_np = th1_to_hist(sig_th1)
    vals = sig_np['vals']
    errs = sig_np['errs']

    outfile = f"{outdir}/f_sigfit.png"
    with quick_ax(outfile=outfile) as ax:
        colors = get_color_cycle()
        pcolor = next(colors)
        ax.errorbar(mtpts, vals, yerr=errs, label="MC", color=pcolor)

        for fitresult, chi2 in fitresults:
            f_name = "f_mc"
            npar = len([f for f in fitresult.floatParsFinal() if f_name in f.GetName()])-1
            f_mc = rl.BasisPoly(f_name, (npar,), ["mt"])
            f_mc.update_from_roofit(fitresult)
            f_mc_vals, f_mc_band = f_mc(mtscaled, nominal=True, errorband=True)

            pcolor = next(colors)
            ax.plot(mtpts, sigyield * f_mc_vals, label=f"fit (n={npar+1}, $\chi^2$={chi2:.1f})", color=pcolor)
            ax.fill_between(mtpts, sigyield * f_mc_band[0], sigyield * f_mc_band[1], alpha=0.2, color=pcolor)
        ax.legend(fontsize=18, framealpha=0.0)
        ax.set_xlabel(r'$m_{\mathrm{T}}$ [GeV]')
        ax.set_ylabel(f'Distribution')
        ax.set_ylim(0., max(vals)*1.1)

def plot_chi2(results, outdir):
    outfile = f"{outdir}/chi2_sigfit.png"
    with quick_ax(outfile=outfile) as ax:
        nvals = [r[0] for r in results]
        chi2vals = [r[1] for r in results]
        ax.plot(nvals, chi2vals, marker='o')
        ax.set_yscale('log')
        ax.set_xlabel('n')
        ax.set_ylabel('$\chi^2$')
        ax.set_xticks([0,5,10,15,20])
        ax.minorticks_on()

npars_to_plot = [6, 10, 14, 16]
sigfits = [(sm["fit"],sm["results"][1]) for sm in sigmodels if any(x==sm["results"][0] for x in npars_to_plot)]
plot_f(mt_array, sig_th1, sigfits, outdir)
plot_chi2(results, outdir)
