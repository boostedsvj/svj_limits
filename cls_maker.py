import os
import boosted_fits as bsvj
import quick_plot as tools
import ROOT as r
import numpy as np

limits = tools.LimitObj()
outfile = bsvj.read_arg('-o', '--outfile', type=str, default='limits.root').outfile

#make a new ttree (limit) in a new root file, write cls into it
file = r.TFile.Open(outfile,"RECREATE")
tree = r.TTree("limit","limit")

base_qtys = ["quantileExpected"]
keys = ["limit"]
qtys = base_qtys + keys + ["trackedParam_{}".format(q) for q in ["mZprime","mDark","rinv","xsec"]]
r.gROOT.ProcessLine("struct quantile_t { "+" ".join(["Double_t {};".format(qty) for qty in qtys])+" };")

qobj = r.quantile_t()

for qty in qtys:
    tree.Branch(qty, r.addressof(qobj,qty), '{}/D'.format(qty))

for key,result in limits.results.items():
    limit = result['limit']
    limit_map = {limit.expected: 0.5, limit.observed:-1, limit.twosigma_down:0.975, limit.onesigma_down:0.84, limit.onesigma_up: 0.16, limit.twosigma_up:0.025}
    for lim,quant in limit_map.items():
        qobj.limit = lim
        qobj.quantileExpected = quant
        qobj.trackedParam_mZprime = key[0]
        qobj.trackedParam_mDark = key[1]
        qobj.trackedParam_rinv = key[2]
        qobj.trackedParam_xsec = bsvj.get_xs(mz)
        tree.Fill()
file.Write()
file.Close()
