import ROOT
import matplotlib.pyplot as plt
import numpy as np

def main():
    mz=[200,250,300,350,400,450,500,550]
    #idir = '/home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits/toyfits_Sep20/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    #idir = '/home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits/toyfits_Sep24/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    #idir = '/home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits/bias_genaltfitua2/toyfits_Oct28_mz'
    idir = '/home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits/bias_genmainfitua2_siginj0p2/hadd/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    #idir = '/home/snabili/data/svj_limits/CMSSW_11_3_4/src/boosted/svj_limits/selftest_siginj0p2/hadd/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
    path = [idir+str(i)+"_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.root" for i in mz]
    #path = [idir + str(i) + "_genaltfitua2/fitDiagnosticsObserveddc_SVJ_s-channel_mMed-"+str(i)+"_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.root" for i in mz]
    mean=[]
    emean=[]
    sigma=[]
    esigma=[]
    lab = "$r_{inj}$=0.2"
    i=0
    for f in path:
        afile = ROOT.TFile(f)
        atree = afile.Get("tree_fit_sb")
        histo = ROOT.TH1F("histo", "histo", 50, -50, 50)
        atree.Draw("(r-0.2)/rErr>> %s"%histo.GetName(), "fit_status==0 || fit_status==1")
        fit = histo.Fit("gaus","S","Q")
        #fit = histo.Fit("gaus",-5,5)
        mean.append(fit.Parameter(1))
        emean.append(fit.ParError(1))
        sigma.append(fit.Parameter(2))
        esigma.append(fit.ParError(2))
    plt.errorbar(mz,mean,yerr=emean, marker='.',markersize=8,linestyle='--',label=f'{lab:s}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
    plt.ylabel('mean',fontsize=12)
    plt.title('Bias Test: Mean of Gaussian Fit to r/rErr',fontsize=15)
    #plt.ylim(-1.5,1.5)
    plt.fill_between(plt.xlim(), 0.5, -0.5,color='#f88379',alpha=0.25, label='$\pm 0.5$')
    plt.savefig('bias_pull_mean_siginj0p2_genmainfitua2_rpm2.pdf')
    plt.savefig('bias_pull_mean_siginj0p2_genmainfitua2_rpm2.png')
    plt.close()

    plt.errorbar(mz,sigma,yerr=esigma, marker='.',markersize=8,linestyle='--',label=f'{lab:s}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
    plt.ylabel('$\sigma$',fontsize=12)
    plt.title('Bias Test: $\sigma$ of Gaussian Fit to (r-$r_{inj})$/rErr',fontsize=15)
    plt.ylim(1.5,5.)
    plt.savefig('bias_pull_stdev_siginj0p2_genmainfitua2_rpm2.pdf')
    plt.savefig('bias_pull_stdev_siginj0p2_genmainfitua2_rpm2.png')
    plt.close()


main()

