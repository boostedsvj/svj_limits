import ROOT as rt
import argparse
import numpy as np

# code to retrive bkg fitting function parameter values
# to run:        python parvalues.py -mz 300 -md 10 -r 0p3 -d 20241010

argParser = argparse.ArgumentParser()
argParser.add_argument("-mz", "--mMed", help="mediator mass")
argParser.add_argument("-md", "--mDark", help="darak hadron mass")
argParser.add_argument("-r", "--rinv", help="rinv")
argParser.add_argument("-d", "--date", help="date for DC folder")

args = argParser.parse_args()
print("args=%s" % args)

#file_ = rt.TFile('dc_'+str(args.date)+'_cutbased/dc_svj_mMed-350_mDark-10_rinv-0p3_nosys.root')
file_ = rt.TFile('dc_'+str(args.date)+'_cutbased/dc_SVJ_s-channel_mMed-'+str(args.mMed)+'_mDark-'+str(args.mDark)+'_rinv-'+str(args.rinv)+'_alpha-peak_MADPT300_13TeV-madgraphMLM-pythia8_sel-cutbased_smooth.root','read')

svj=file_.Get('SVJ')

varsname = svj.allVars()

iter = varsname.createIterator()
var = iter.Next()

pars = []
while var:
   pars.append(var.GetName())
   var = iter.Next()


for p in pars:
   print(svj.var(p))
