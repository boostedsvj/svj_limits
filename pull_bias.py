import ROOT
import matplotlib.pyplot as plt
import numpy as np
import os
import glob



def main():
    mz = np.linspace(200,550,8)
    #dir = 'badcondor_BiasSelftests/bias_genmainfitua2mod_siginjExp'
    #dir = 'self_genua2modfitua2mod_siginj0p0_4par'
    dir = 'bias_genmainfitua2mod_siginj0p0'
    idir = os.path.join(os.getcwd(),dir)
    fitDia_files = glob.glob(idir+'/hadd/fit*')
    Genera_files = glob.glob(idir + '/toysgen*/higg*smooth0*.root')
    print('*'*50, [np.round(ROOT.TFile(g).Get("w").var('r').getVal(),2) for g in Genera_files],'*'*50)
    mean=[]
    emean=[]
    sigma=[]
    esigma=[]
    #lab = "$r_{inj}$=$r_{exp}$"
    lab = "$r_{exp}=0$"
    i=0
    for f,g in zip(fitDia_files,Genera_files):
        exp_sig = np.round(ROOT.TFile(g).Get("w").var('r').getVal(),2)
        print(exp_sig)
        afile = ROOT.TFile(f)
        atree = afile.Get("tree_fit_sb")
        histo = ROOT.TH1F("histo", "histo", 30, -10, 10)
        atree.Draw("(r-"+str(exp_sig)+")/rErr>> %s"%histo.GetName(), "fit_status==0 || fit_status==1")
        fit = histo.Fit("gaus","S","Q")
        mean.append(fit.Parameter(1))
        emean.append(fit.ParError(1))
        sigma.append(fit.Parameter(2))
        esigma.append(fit.ParError(2))
        i+=1
    plt.errorbar(mz,mean,yerr=emean, marker='.',markersize=8,linestyle='--',label=f'{lab:s}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
    plt.ylabel('mean',fontsize=12)
    plt.title('Bias Test: Mean of Gaussian Fit to (r-$r_{exp}$)/rErr',fontsize=15)
    #plt.title('Self Test: Mean of Gaussian Fit to r/rErr',fontsize=15)
    #plt.ylim(-1.5,1.5)
    plt.fill_between(plt.xlim(), 0.5, -0.5,color='#f88379',alpha=0.25, label='$\pm 0.5$')
    plt.savefig('biastest/pull_plots/mean_siginjExp_ua2mod_4pars.png')
    #plt.savefig('selftest/pull_plots/mean_expmod3.png')
    plt.close()

    plt.errorbar(mz,sigma,yerr=esigma, marker='.',markersize=8,linestyle='--',label=f'{lab:s}')
    plt.legend()
    plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
    plt.ylabel('$\sigma$',fontsize=12)
    plt.title('Bias Test: $\sigma$ of Gaussian Fit to (r-$r_{exp}$)/rErr',fontsize=15)
    #plt.title('Self Test: $\sigma$ of Gaussian Fit to r/rErr',fontsize=15)
    #plt.ylim(1.5,5.)
    plt.savefig('biastest/pull_plots/stdev_siginjExp_ua2mod_4pars.png')
    #plt.savefig('selftest/pull_plots/stdev_siginj0p0_ua2mod_3pars.png')
    plt.close()


main()

