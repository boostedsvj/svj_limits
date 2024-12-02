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

ROOT.gROOT.SetBatch(True)  # Prevents ROOT from opening canvases

def main():

    parser = argparse.ArgumentParser(description="Configure base directory and parameters for bias test")
    parser.add_argument('--base_dir', type=str, required=True,
                        help=("Directory path contianing the test subdirectories, e.g. '/path/to/test/'," 
                              " it expects two subdirectories: siginj0 and siginj1"))
    parser.add_argument('--sel', type=str, required=True, choices=['cutbased', 'bdt=0.65', 'bdt=0.67'],
                        help="Selection type: 'cutbased' or 'bdt=0.67' or add additional bdt values as needed")
    parser.add_argument('--test', type=str, required=True, choices=['bias', 'self'],
                        help="Test type: 'bias' or 'self'")
    parser.add_argument('--mz', nargs='+', type=int, default=[200, 250, 300, 350, 400, 450, 500, 550],
                        help="List of mass points. Default is [200, 250, 300, 350, 400, 450, 500, 550]")
    parser.add_argument('--inj_value', type=float, default=None,
                        help="Set all injection values to a specific number for testing (e.g., 0.3)")   
 
    args = parser.parse_args()

    # Switch between bias and self test
    test = args.test
    test_title = None
    if test == 'bias' : test_title = 'Bias'
    elif test == 'self' : test_title = 'Self'

    # injected signal strength
    # Needs to be added for cutbased
    if args.inj_value is not None:
        # If the user specifies --inj_value, override all injection values
        inj = {mz: args.inj_value for mz in args.mz}
    elif args.sel == 'bdt=0.67':
        # Expected limit values for bdt=0.67
        inj_values = [0.271, 0.136, 0.165, 0.189, 0.210, 0.249, 0.265, 0.397]
        inj = {mz: inj_val for mz, inj_val in zip(args.mz, inj_values)}
    else:
        # Default case with warning
        print("WARNING: using 0.2 as injected signal value, please correct if unintended")
        inj = {mz: 0.2 for mz in args.mz}

    # Directories for r_inj = 0 and r_inj = inj
    path_rinj0 = [(f"{args.base_dir}/siginj0/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-{mz}_mDark-10_"
                   f"rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-{args.sel}_smooth.root") for mz in args.mz]
    path_rinj1 = [(f"{args.base_dir}/siginj1/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-{mz}_mDark-10_"
                   f"rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-{args.sel}_smooth.root") for mz in args.mz]

    mean_rinj0, emean_rinj0, sigma_rinj0, esigma_rinj0 = [], [], [], []
    mean_rinj1, emean_rinj1, sigma_rinj1, esigma_rinj1 = [], [], [], []

    for i, (mz, f_rinj0, f_rinj1) in enumerate(zip(args.mz, path_rinj0, path_rinj1)):
        afile_rinj0 = ROOT.TFile(f_rinj0)
        atree_rinj0 = afile_rinj0.Get("tree_fit_sb")
        
        afile_rinj1 = ROOT.TFile(f_rinj1)
        atree_rinj1 = afile_rinj1.Get("tree_fit_sb")
        
        # Create histograms for r/rErr for both r_inj=0 and r_inj=1
        histo_rinj0 = ROOT.TH1F(f"histo_rinj0_{mz}", f"Histogram of (r-rinj)/(0.5*(rHiErr+rLoErr)) for MZ = {mz} GeV (r_inj=0)", 50, -10, 10)
        histo_rinj1 = ROOT.TH1F(f"histo_rinj1_{mz}", f"Histogram of (r-rinj)/(0.5*(rHiErr+rLoErr)) for MZ = {mz} GeV (r_inj={inj[mz]})", 50, -10, 10)
        
        # Draw r/rErr for both r_inj=0 and r_inj=1 into their respective histograms
        atree_rinj0.Draw("(r)/rErr>>%s" % histo_rinj0.GetName(), "fit_status==0 || fit_status==1")
        atree_rinj1.Draw(f"(r - {inj[mz]})/rErr>>%s" % histo_rinj1.GetName(), "fit_status==0 || fit_status==1")
        
        # Perform Gaussian fits for both histograms
        fit_rinj0 = histo_rinj0.Fit("gaus", "S", "Q")
        fit_rinj1 = histo_rinj1.Fit("gaus", "S", "Q")
        
        # Store means and sigmas for both fits
        mean_rinj0.append(fit_rinj0.Parameter(1))
        emean_rinj0.append(fit_rinj0.ParError(1))
        sigma_rinj0.append(fit_rinj0.Parameter(2))
        esigma_rinj0.append(fit_rinj0.ParError(2))
        
        mean_rinj1.append(fit_rinj1.Parameter(1))
        emean_rinj1.append(fit_rinj1.ParError(1))
        sigma_rinj1.append(fit_rinj1.Parameter(2))
        esigma_rinj1.append(fit_rinj1.ParError(2))

        # Save the histograms as PDFs and PNGs
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        histo_rinj0.Draw()
        canvas.SaveAs(f"{test}_siginj_rinj0_MZ_{mz}_GeV.pdf")
        canvas.SaveAs(f"{test}_siginj_rinj0_MZ_{mz}_GeV.png")
        
        histo_rinj1.Draw()
        canvas.SaveAs(f"{test}_siginj_rinj1_MZ_{mz}_GeV.pdf")
        canvas.SaveAs(f"{test}_siginj_rinj1_MZ_{mz}_GeV.png")

    # Plot mean for both r_inj = 0 and r_inj = 1
    plt.errorbar(args.mz, mean_rinj0, yerr=emean_rinj0, marker='.', markersize=8, linestyle='--', label='$r_{inj}=0$')
    plt.errorbar(args.mz, mean_rinj1, yerr=emean_rinj1, marker='.', markersize=8, linestyle='--', label=f'$r_{{inj}}=Expected$')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('Mean', fontsize=12)
    plt.title(f'{test_title} Test: Mean of Gaussian Fit to (r-$r_{{inj}})$/rErr', fontsize=15)
    plt.ylim(-1.5, 1.5)
    plt.fill_between(plt.xlim(), 0.5, -0.5, color='#f88379', alpha=0.25, label='$\pm 0.5$')
    plt.savefig(f'{test}_pull_mean.png')
    plt.savefig(f'{test}_pull_mean.pdf')
    plt.close()

    # Plot sigma for both r_inj = 0 and r_inj = 1
    plt.errorbar(args.mz, sigma_rinj0, yerr=esigma_rinj0, marker='.', markersize=8, linestyle='--', label='$r_{inj}=0$')
    plt.errorbar(args.mz, sigma_rinj1, yerr=esigma_rinj1, marker='.', markersize=8, linestyle='--', label=f'$r_{{inj}}=Expected$')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('$\sigma$', fontsize=12)
    plt.title(f'{test_title} Test: $\sigma$ of Gaussian Fit to (r-$r_{{inj}})$/rErr', fontsize=15)
    plt.ylim(0.5, 2.5)
    plt.savefig(f'{test}_pull_stdev.png')
    plt.savefig(f'{test}_pull_stdev.pdf')
    plt.close()

main()

