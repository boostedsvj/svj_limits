import os
import boosted_fits as bsvj
import quick_plot as tools
import ROOT as r
import numpy as np


rootfiles = bsvj.pull_arg('rootfiles', type=str, nargs='+').rootfiles
clean  = bsvj.pull_arg('--clean', action='store_true').clean
mz     = bsvj.pull_arg('--mz',type=float).mz
rinv   = bsvj.pull_arg('--rinv',type=float).rinv


#make a new ttree (limit) in a new root file, write cls into it 
file = r.TFile.Open("Limit_mz{}_rinv{}.root".format(mz,rinv),"RECREATE")
tree = r.TTree("limit","limit")

base_qtys = ["observed"]
keys = ["expected", "twosigma_down","onesigma_down","onesigma_up","twosigma_up"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in ["mZprime","mDark","rinv","xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")

qobj = r.quantile_t()

for qty in qtys:
  tree.Branch(qty, r.addressof(qobj,qty), '{}/D'.format(qty))


for observed, asimov in zip(*tools.organize_rootfiles(rootfiles)):
  obs, asimov = tools.extract_scans([observed, asimov], correct_minimum=True)
  if clean:
    obs    = tools.clean_scan(obs)
    asimov = tools.clean_scan(asimov)
  cls = tools.get_cls(obs, asimov)
limit = tools.interpolate_95cl_limit(cls)

observed = rootfiles[1]
f = r.TFile(observed,'read')
t = f.Get("limit")
N = t.GetEntries()
for i in range(N):
  t.GetEntry(i)
  qobj.trackedParam_mZprime = t.trackedParam_mZprime
  qobj.trackedParam_mDark = t.trackedParam_mDark
  qobj.trackedParam_rinv = t.trackedParam_rinv
  qobj.trackedParam_xsec = t.trackedParam_xsec
  qobj.observed = limit.observed
  qobj.expected = limit.expected
  qobj.twosigma_down = limit.twosigma_up
  qobj.onesigma_down = limit.onesigma_up
  qobj.onesigma_up = limit.onesigma_down
  qobj.twosigma_up = limit.twosigma_down
  tree.Fill()
  
file.Write()
file.Close()
