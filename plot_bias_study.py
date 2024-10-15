import ROOT
import matplotlib
matplotlib.use('Agg') # prevents opening displays (fast), must use before pyplot
import matplotlib.pyplot as plt
import numpy as np

ROOT.gROOT.SetBatch(True)  # Prevents ROOT from opening canvases

def main():
    mz = [200, 250, 300, 350, 400, 450, 500, 550]
    
    # Directories for r_inj = 0 and r_inj = 1
    idir_rinj0 = '/uscms/home/bregnery/work/CMSSW_11_3_4/src/boosted/svj_limits/bias_test/siginj0/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    idir_rinj1 = '/uscms/home/bregnery/work/CMSSW_11_3_4/src/boosted/svj_limits/bias_test/siginj1/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    
    path_rinj0 = [idir_rinj0 + str(i) + "_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt=0.65_smooth.root" for i in mz]
    path_rinj1 = [idir_rinj1 + str(i) + "_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-bdt=0.65_smooth.root" for i in mz]

    mean_rinj0, emean_rinj0, sigma_rinj0, esigma_rinj0 = [], [], [], []
    mean_rinj1, emean_rinj1, sigma_rinj1, esigma_rinj1 = [], [], [], []

    for i, (f_rinj0, f_rinj1) in enumerate(zip(path_rinj0, path_rinj1)):
        afile_rinj0 = ROOT.TFile(f_rinj0)
        atree_rinj0 = afile_rinj0.Get("tree_fit_sb")
        
        afile_rinj1 = ROOT.TFile(f_rinj1)
        atree_rinj1 = afile_rinj1.Get("tree_fit_sb")
        
        # Create histograms for r/rErr for both r_inj=0 and r_inj=1
        histo_rinj0 = ROOT.TH1F(f"histo_rinj0_{mz[i]}", "Histogram of (r-rinj)/(0.5*(rHiErr+rLoErr)) for MZ = %d GeV (r_inj=0)" % mz[i], 50, -10, 10)
        histo_rinj1 = ROOT.TH1F(f"histo_rinj1_{mz[i]}", "Histogram of (r-rinj)/(0.5*(rHiErr+rLoErr)) for MZ = %d GeV (r_inj=1)" % mz[i], 50, -10, 10)
        
        # Draw r/rErr for both r_inj=0 and r_inj=1 into their respective histograms
        atree_rinj0.Draw("(r)/rErr>>%s" % histo_rinj0.GetName(), "fit_status==0 || fit_status==1")
        atree_rinj1.Draw("(r - 1)/rErr>>%s" % histo_rinj1.GetName(), "fit_status==0 || fit_status==1")
        
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
        canvas.SaveAs(f"bias_siginj_rinj0_MZ_{mz[i]}_GeV.pdf")
        canvas.SaveAs(f"bias_siginj_rinj0_MZ_{mz[i]}_GeV.png")
        
        histo_rinj1.Draw()
        canvas.SaveAs(f"bias_siginj_rinj1_MZ_{mz[i]}_GeV.pdf")
        canvas.SaveAs(f"bias_siginj_rinj1_MZ_{mz[i]}_GeV.png")

    # Plot mean for both r_inj = 0 and r_inj = 1
    plt.errorbar(mz, mean_rinj0, yerr=emean_rinj0, marker='.', markersize=8, linestyle='--', label='$r_{inj}=0$')
    plt.errorbar(mz, mean_rinj1, yerr=emean_rinj1, marker='.', markersize=8, linestyle='--', label='$r_{inj}=1$')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('Mean', fontsize=12)
    plt.title('Bias Test: Mean of Gaussian Fit to (r-$r_{inj})$/rErr', fontsize=15)
    plt.ylim(-1.5, 1.5)
    plt.fill_between(plt.xlim(), 0.5, -0.5, color='#f88379', alpha=0.25, label='$\pm 0.5$')
    plt.savefig('bias_pull_mean.png')
    plt.savefig('bias_pull_mean.pdf')
    plt.close()

    # Plot sigma for both r_inj = 0 and r_inj = 1
    plt.errorbar(mz, sigma_rinj0, yerr=esigma_rinj0, marker='.', markersize=8, linestyle='--', label='$r_{inj}=0$')
    plt.errorbar(mz, sigma_rinj1, yerr=esigma_rinj1, marker='.', markersize=8, linestyle='--', label='$r_{inj}=1$')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)', fontsize=12)
    plt.ylabel('$\sigma$', fontsize=12)
    plt.title('Bias Test: $\sigma$ of Gaussian Fit to (r-$r_{inj})$/rErr', fontsize=15)
    plt.ylim(0.5, 2.5)
    plt.savefig('bias_pull_stdev.png')
    plt.savefig('bias_pull_stdev.pdf')
    plt.close()

main()

