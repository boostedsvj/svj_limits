import ROOT as rt
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

# to run: python pull_bias_compare.py -nosig toyfits_Sep24/ -sig toyfits_Oct09/ 

def main():
	argParser = argparse.ArgumentParser()
	argParser.add_argument("-nosig",  "--nosig_toydir", help="nosiginj directory")
	argParser.add_argument("-sig",    "--sig_toydir",   help="siginj   directory")
	args = argParser.parse_args()

	mz = np.arange(200,551,50)
	nosig_idir = os.getcwd() + '/' + args.nosig_toydir + 'fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
	sig_idir   = os.getcwd() + '/' + args.sig_toydir   + 'fitDiagnosticsObserveddc_SVJ_s-channel_mMed-'
	nosig_files = [nosig_idir + str(i) + '_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.root' for i in mz]
	sig_files   = [sig_idir   + str(i) + '_mDark-10_rinv-0p3_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.root' for i in mz]
	fitvalues = {}
	#sig_fitvalues   = {}
	nsinj = np.zeros(8)
	sinj  = [0.4410, 0.3895,0.3481,0.2918,0.3741,0.2,0.6364,0.7172]
	SINJ = {0:nsinj, 1:sinj}
	l=0
	for s in (nosig_files, sig_files):
		fitvalues[l] = {}
		mean, sigma, emean, esigma = [],[],[],[]
		i=0
		for f in s:
			afile = rt.TFile(f)
			atree = afile.Get("tree_fit_sb")
			histo = rt.TH1F("histo", "histo", 50, -10, 10)
			sig = SINJ[l][i]
			atree.Draw("(r-%f)/rErr>> %s"%(SINJ[l][i],histo.GetName()), "fit_status==0 || fit_status==1")
			#atree.Draw("(r-"+str(sinj)+")/rErr>> %s"%histo.GetName(), "fit_status==0 || fit_status==1")
			fit = histo.Fit("gaus","S","Q")
			mean.append(fit.Parameter(1))
			emean.append(fit.ParError(1))
			sigma.append(fit.Parameter(2))
			esigma.append(fit.ParError(2))
			i+=1
		print(mean,emean,sigma,esigma)
		fitvalues[l]['mean']   = mean
		fitvalues[l]['emean']  = emean
		fitvalues[l]['sigma']  = sigma
		fitvalues[l]['esigma'] = esigma
                #print(fitvalues)
		print(fitvalues)#[l][j] for l in [0,1] for j in ['mean','emean','sigma','esigma']])
		l+=1
	plt.errorbar(mz,fitvalues[0]['mean'],yerr=fitvalues[0]['emean'], marker='.',markersize=8,linestyle='--',label=r'$r_{inj}$=0')
	plt.errorbar(mz,fitvalues[1]['mean'],yerr=fitvalues[1]['emean'], marker='.',markersize=8,linestyle='--',label='sig inj')
	plt.legend()
	plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
	plt.ylabel('mean',fontsize=12)
	plt.title(r'Bias Test: Mean of Gaussian Fit to (r-$r_{inj}$)/rErr',fontsize=15)
	plt.ylim(-1.5,1.5)
	plt.fill_between(plt.xlim(), 0.5, -0.5,color='#f88379',alpha=0.25, label='$\pm 0.5$')
	plt.savefig('bias_pull_mean_siginj_comparison.pdf')
	plt.close()

	plt.errorbar(mz,fitvalues[0]['sigma'],yerr=fitvalues[0]['esigma'], marker='.',markersize=8,linestyle='--',label=r'$r_{inj}$=0')
	plt.errorbar(mz,fitvalues[1]['sigma'],yerr=fitvalues[1]['esigma'], marker='.',markersize=8,linestyle='--',label='sig inj')
	plt.legend()
	plt.xlabel('$M_{Z^{\prime}}$ (GeV)',fontsize=12)
	plt.ylabel('$\sigma$',fontsize=12)
	plt.title(r'Bias Test: $\sigma$ of Gaussian Fit to (r-$r_{inj}$)/rErr',fontsize=15)
	plt.ylim(0.5,3.5)
	plt.savefig('bias_pull_stdev_siginj_comparison.pdf')
	plt.close()


main()
