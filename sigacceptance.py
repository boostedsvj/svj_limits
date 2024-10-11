import os
import boosted_fits as bsvj
import quick_plot as tools
import ROOT as r
import numpy as np


rootfiles = bsvj.pull_arg('rootfiles', type=str).rootfiles
clean  = bsvj.pull_arg('--clean', action='store_true').clean
mz     = bsvj.pull_arg('--mz',type=float).mz
rinv   = bsvj.pull_arg('--rinv',type=float).rinv


#make a new ttree (limit) in a new root file, write cls into it 
file = r.TFile.Open("SigAcc/SigAccEff_mz{}_rinv{}.root".format(mz,rinv),"RECREATE")
tree = r.TTree("limit","limit")

base_qtys = ["quantileExpected"]
keys = ["limit"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in ["mZprime","mDark","rinv","xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")

qobj = r.quantile_t()

for qty in qtys:
  tree.Branch(qty, r.addressof(qobj,qty), '{}/D'.format(qty))
with bsvj.open_root(rootfiles) as f:
  ws = bsvj.get_ws(f)
sig = ws.data('sig')
_, y_sig, _ = bsvj.roodataset_values(sig)
lumi=21071

print(sum(y_sig)/(lumi * tools.get_xsec(mz)))
qobj.limit = sum(y_sig)/(lumi * tools.get_xsec(mz))
qobj.quantileExpected = 0.5
qobj.trackedParam_mZprime = mz
qobj.trackedParam_mDark = 10
qobj.trackedParam_rinv = rinv
qobj.trackedParam_xsec = tools.get_xsec(mz)
tree.Fill()

file.Write()
file.Close()
