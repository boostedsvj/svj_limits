#!/usr/bin/env python3

import argparse
import ROOT as r

def fprint(*args):
    import sys
    print(*args)
    sys.stdout.flush()

def obj_print(obj, parg):
    obj.Print(parg)
    if hasattr(obj,'numEntries'):
        for i in range(obj.numEntries()):
            obj.get(i)
            w = obj.weight()
            err = r.TMath.Sqrt(obj.weightSquared())
            fprint(f"  fSumw[{i}]={w}, error={err}")

def main(filename, workspace, verbose):
    file = r.TFile.Open(filename)
    ws = file.Get(workspace)
    parg = "V" if verbose else ""

    ws.Print(parg)

    alls = ['allCats', 'allData', 'allFunctions', 'allGenericObjects', 'allPdfs', 'allVars']
    for all in alls:
        all_fn = getattr(ws,all)
        fprint('\n'+'*'*10+f" {all} "+'*'*10+'\n')
        for obj in sorted(all_fn(), key=lambda obj: obj.GetName()):
            obj_print(obj, parg)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, required=True, help="input filename")
    parser.add_argument("-w", "--workspace", type=str, required=True, help="input workspace name")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="enable verbose printouts")
    args = parser.parse_args()

    main(args.file, args.workspace, args.verbose)
