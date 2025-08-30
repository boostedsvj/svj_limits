#==============================================================================
# plot_bias_or_self_study.py ----------------------------------------------------------
#------------------------------------------------------------------------------
# Author(s): Brendan Regnery, Sara Nabili, Kevin Pedro ------------------------
#------------------------------------------------------------------------------
# Plot the results of the bias test, it's a rather simple script made from ----
#     an example given by Sara. It is kept here for smooth running in later ---
#     reruns of the analysis. -------------------------------------------------
#------------------------------------------------------------------------------

import ROOT
import matplotlib
matplotlib.use('Agg') # prevents opening displays (fast), must use before pyplot
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

from run_rhalpha import handle_signals_explim, get_signame, get_rinj

ROOT.gROOT.SetBatch(True)  # Prevents ROOT from opening canvases
ROOT.TH1.AddDirectory(True)

def main():
    imgs = ["pdf", "png"]

    parser = argparse.ArgumentParser(description="Configure base directory and parameters for bias test")
    parser.add_argument('--base-dir', type=str, required=True,
                        help=("Directory path containing the test subdirectories, e.g. '/path/to/test/',"
                              " it expects two subdirectories: rinj0 and rinj1"))
    parser.add_argument('--suff', type=str, default="",
                        help="Suffix for rinj0 and rinj1 dir names")
    parser.add_argument('--sel', type=str, required=True,
                        help="Selection type")
    parser.add_argument('--test', type=str, required=True, choices=['bias', 'self'],
                        help="Test type: 'bias' or 'self'")
    parser.add_argument("--signals", dest="signals", type=str, default="",
                        help="text file w/ list of signal parameters")
    parser.add_argument('--inj-types', type=int, choices=[0,1], nargs='*', default=[0,1],
                        help="which signal injection type(s) to plot")
    parser.add_argument("--explim", type=str, default="explim_{sel}.txt",
                        help="generated file with expected limit values per signal from Asimov scan")
    parser.add_argument("--rinj", type=float, default=0,
                        help="toy signal injection strength; if -x, strength = expected limit * x")
    parser.add_argument('-s', '--seed', type=int, default=1001,
                        help="toy seed")
    args = parser.parse_args()
    args, signals, _ = handle_signals_explim(args)

    # Switch between bias and self test
    test = args.test
    test_title = None
    if test == 'bias' : test_title = 'Bias'
    elif test == 'self' : test_title = 'Self'

    # Directories for r_inj = 0 and r_inj = inj
    results = {}
    for inj_type in args.inj_types:
        base = f"{args.base_dir}/rinj{inj_type}{args.suff}"
        if not os.path.exists(base):
            print(f"inj_type {inj_type} not found at {base}, skipping")
            continue

        rinj = 0 if inj_type==0 else args.rinj

        results[inj_type] = {"mean": [], "emean": [], "sigma": [], "esigma": []}

        for signal in signals:
            rinj_signal = get_rinj(rinj,signal)
            path = f"{base}/higgsCombineObserveddc_{get_signame(signal)}_sel-{args.sel}_mt_smooth.FitDiagnostics.mH120.{args.seed}.root"
            afile = ROOT.TFile(path)
            atree = afile.Get("tree_fit_sb")

            # Create histograms for r/rErr
            histo = ROOT.TH1F(f"histo_rinj{inj_type}_{signal.mMed}", f"Histogram of (r-rinj)/(0.5*(rLoErr+rHiErr)) for MZ = {signal.mMed} GeV (r_inj={rinj_signal})", 50, -10, 10)

            # Draw r/rErr into the histogram
            atree.Draw(f"(r - {rinj_signal})/(0.5*(rLoErr+rHiErr))>>{histo.GetName()}", "fit_status==0 || fit_status==1")

            # Perform Gaussian fit
            fit = histo.Fit("gaus", "S", "Q")

            # Store means and sigmas for both fits
            results[inj_type]["mean"].append(fit.Parameter(1))
            results[inj_type]["emean"].append(fit.ParError(1))
            results[inj_type]["sigma"].append(fit.Parameter(2))
            results[inj_type]["esigma"].append(fit.ParError(2))

            # Save the histograms as PDFs and PNGs
            canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
            histo.Draw()
            for img in imgs:
                canvas.SaveAs(f"{args.base_dir}/{test}_rinj{inj_type}_MZ_{signal.mMed}_GeV.{img}")

    # keep consistent colors
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Plot mean for both r_inj = 0 and r_inj = 1
    inj_label = ["0", "Expected"]
    xvals = [signal.mMed_val() for signal in signals]
    for inj_type in results:
        plt.errorbar(xvals, results[inj_type]["mean"], yerr=results[inj_type]["emean"], color=color_cycle[inj_type], marker='.', markersize=8, linestyle='--', label=f'$r_{{inj}}=${inj_label[inj_type]}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('Mean', fontsize=12)
    plt.title(f'{test_title} Test: Gaussian $\mu$ of $(r-r_{{inj}})/\\langle\\varepsilon_{{r}}\\rangle$', fontsize=15)
    plt.ylim(-3.0, 3.0)
    plt.fill_between(plt.xlim(), 0.5, -0.5, color='#f88379', alpha=0.25, label='$\pm 0.5$')
    for img in imgs:
        plt.savefig(f'{args.base_dir}/{test}_pull_mean.{img}')
    plt.close()

    # Plot sigma for both r_inj = 0 and r_inj = 1
    for inj_type in results:
        plt.errorbar(xvals, results[inj_type]["sigma"], yerr=results[inj_type]["esigma"], color=color_cycle[inj_type], marker='.', markersize=8, linestyle='--', label=f'$r_{{inj}}=${inj_label[inj_type]}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('$\sigma$', fontsize=12)
    plt.title(f'{test_title} Test: Gaussian $\sigma$ of $(r-r_{{inj}})/\\langle\\varepsilon_{{r}}\\rangle$', fontsize=15)
    plt.ylim(0.5, 2.5)
    for img in imgs:
        plt.savefig(f'{args.base_dir}/{test}_pull_stdev.{img}')
    plt.close()

main()

