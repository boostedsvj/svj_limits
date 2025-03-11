import numbers
import warnings
import numpy as np
import argparse
import rhalphalib as rl

# from hessian.py
def get_objs(args):
    import ROOT

    fws, wsname = args.workspace.split(':')
    fin = ROOT.TFile.Open(fws)
    w = fin.Get(wsname)

    ffit, fitname = args.fit.split(':')
    fin2 = ROOT.TFile.Open(ffit)
    fit = fin2.Get(fitname)

    return w, fit

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Get decorrelated parameter expressions from covariance matrix",
    )
    subparsers = parser.add_subparsers()

    p_hessian = subparsers.add_parser('hessian')
    p_hessian.add_argument("-w", "--workspace", metavar="ROOTFILE:WORKSPACE", help="Workspace to load", required=True)
    p_hessian.add_argument("-f", "--fit", metavar="ROOTFILE:FIT_NAME", help="Fit result to load", required=True)
    p_hessian.add_argument("-m", "--model", help="Model to load", default="ModelConfig")
    p_hessian.add_argument("-p", "--prefix", help="Prefix for parameter names", default="bsvj")
    p_hessian.set_defaults(mode='hessian')

    p_cov = subparsers.add_parser('cov')
    p_cov.add_argument("-t", "--tree", metavar="ROOTFILE:TREE_NAME", help="Tree to load", required=True)
    p_cov.add_argument("-p", "--parameters", nargs='+', help="Parameters to use", required=True)
    p_cov.set_defaults(mode='cov')

    args = parser.parse_args()

    if args.mode=='hessian':
        w, fit = get_objs(args)

        param_names = [p.GetName() for p in fit.floatParsFinal() if p.GetName().startswith(args.prefix)]
        decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult("@", fit, param_names)
    elif args.mode=='cov':
        import uproot as up
        fname, tname = args.tree.split(':')
        f = up.open(fname)
        pnames = args.parameters
        params = f[tname].arrays(pnames)

        # compute covariance matrix from scan in limit tree
        X = np.asarray([params[p] for p in pnames])
        means = np.mean(X, axis=1)
        cov = np.cov(X)
        decoVector = rl.DecorrelatedNuisanceVector("@", means, cov)

    # output
    print(decoVector.correlated_str.replace("{","").replace("}",""))

