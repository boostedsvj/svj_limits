import argparse
import os
import json
import ROOT as r
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot as up
from matplotlib.offsetbox import AnchoredText
import rhalphalib as rl

verbose = False

def get_fit_val(fitDiag, val, fittype='fit_s', substitute=1.):
    if val in fitDiag.Get(fittype).floatParsFinal().contentsString().split(','):
        return fitDiag.Get(fittype).floatParsFinal().find(val).getVal()
    else:
        return substitute


def prepareTF(
    rl_poly,
    xaxis=(180, 650),
    xstep=10
    ):
    fmtbins = np.arange(xaxis[0], xaxis[1]+xstep, xstep)
    fmtpts = fmtbins[:-1] + 0.5 * np.diff(fmtbins)
    fmtscaled = (fmtpts - min(fmtpts)) / (max(fmtpts) - min(fmtpts))

    TF = rl_poly(fmtscaled, nominal=True)
    if verbose: print(TF)
    return TF, fmtpts

def base_plot(TF, fmtpts, ax=None):
    if ax is None:
        ax = plt.gca()

    plt.plot(fmtpts, TF)

    ax.set_xlim(round(np.min(fmtpts)), round(np.max(fmtpts)))
    ax.set_xlabel(r'$m_{T}$', ha='right', x=1)
    ax.set_ylabel(r'TF', ha='right', y=1)
    return ax


def singleTF(tf_poly):
    TF, fmtpts = prepareTF(tf_poly)
    fig, ax = plt.subplots()
    base_plot(TF, fmtpts, ax=ax)

    _deg_str = "$({})$".format(",".join([n +"=" + str(o) for n, o in zip(tf_poly.dim_names, tf_poly.order)]))
    _label =  "{} {}".format(tf_poly.basis, _deg_str)

    at = AnchoredText(_label, loc='lower right', frameon=False, prop=dict(size=20))
    ax.add_artist(at);
    return ax


def combinedTF(tf1, tf2):
    TF1, fmtpts = prepareTF(tf1)
    TF2, _ = prepareTF(tf2)
    fig, ax = plt.subplots()
    base_plot(TF1 * TF2, fmtpts, ax=ax)

    _deg_str = "$({})$".format(",".join([n +"=" + str(o) for n, o in zip(tf1.dim_names, tf1.order)]))
    _label =  "{} {}".format(tf1.basis, _deg_str)
    _deg_str = "$({})$".format(",".join([n +"=" + str(o) for n, o in zip(tf2.dim_names, tf2.order)]))
    _label += "\nx\n{} {}".format(tf2.basis, _deg_str)

    at = AnchoredText(_label, loc='lower right', frameon=False, prop=dict(size=20, multialignment='center'))
    ax.add_artist(at);
    return ax


def plotTF_ratio(in_ratio, region):
    fig, ax = plt.subplots()
    mtpts = np.arange(180,650+10,10)
    hep.histplot(in_ratio, mtpts)
    ax.tick_params(axis='y', which='minor', left=False, right=False)

    ax.set_title('{} QCD Ratio'.format(region), pad=15, fontsize=26)
    ax.set_xlabel(r'$m_{T}$', ha='right', x=1)
    ax.set_ylabel(r'(Pass QCD) / (Fail QCD * eff)', ha='right', y=1)
    return ax


if __name__ == '__main__':
    plt.style.use([hep.cms.style.ROOT])
    plt.switch_backend('agg')
    np.seterr(divide='ignore', invalid='ignore')

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", default='', help="Model/Fit dir")
    parser.add_argument("-f",
                        "--fit",
                        default='fitDiagnostics.root',
                        dest='fit',
                        help="fitDiagnostics file")
    parser.add_argument("-o",
                        "--output-folder",
                        default='plots',
                        dest='output_folder',
                        help="Folder to store plots - will be created ? doesn't exist.")
    parser.add_argument("--label",
                        default=None,
                        choices={"Preliminary", "Private Work", "Supplementary"},
                        type=str,
                        help="year label")
    parser.add_argument("--MC",
                        action='store_true',
                        dest='isMC',
                        help="Use 'simulation' label")
    args = parser.parse_args()

    if args.output_folder.split("/")[0] != args.dir:
        args.output_folder = os.path.join(args.dir, args.output_folder)
    os.makedirs(args.output_folder, exist_ok=True)

    # Get fitDiagnostics File
    rf = r.TFile.Open(os.path.join(args.dir, args.fit))

    # Define TFs
    basis1, basis2 = ("Bernstein","Bernstein")
    degs = (2,)
    degsMC = (2,)
    tf_MC = rl.BasisPoly("tf_mc", degsMC, ['mt'], basis=basis1)
    tf_res = rl.BasisPoly("tf_data", degs, ['mt'], basis=basis2)

    # Set to values from fit
    tf_res.update_from_roofit(rf.Get('fit_s'))
    pars_final = rf.Get('fit_s').floatParsFinal()
    par_names = pars_final.contentsString().split(',')

    ax = singleTF(tf_res)
    ax.set_title("Residual (Data/MC) TF", x=0, ha='left', fontsize='small')
    hep.cms.label(llabel=args.label, ax=ax, loc=2, data=not args.isMC)
    ax.figure.savefig('{}/TF_data.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_data.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")

    if verbose: print('MC_nuis','\n'.join([f'{p} = {pars_final.find(p).getVal()}' for p in par_names if 'tf_mc' in p]))
    MC_nuis = [round(pars_final.find(p).getVal(), 3) for p in par_names if 'tf_mc' in p]
    _vect = np.load(os.path.join(os.path.dirname(args.dir), 'svjModel_deco.npy'))
    if verbose: print('_vect',_vect)
    _MCTF_nominal = np.load(os.path.join(os.path.dirname(args.dir), 'svjModel_mctf.npy'))
    if verbose: print('_MCTF_nominal',_MCTF_nominal)
    _values = _vect.dot(np.array(MC_nuis)) + _MCTF_nominal
    tf_MC.set_parvalues(_values)

    if verbose: print('singleTF MC')
    ax = singleTF(tf_MC)
    ax.set_title("Tagger Response TF", x=0, ha='left', fontsize='small')
    hep.cms.label(llabel=args.label, ax=ax, loc=2, data=not args.isMC)
    ax.figure.savefig('{}/TF_MC.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_MC.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")

    if verbose: print('combinedTF')
    ax = combinedTF(tf_MC, tf_res)
    ax.set_title("Effective Transfer Factor", x=0, ha='left', fontsize='small')
    hep.cms.label(llabel=args.label, ax=ax, loc=2, data=not args.isMC)
    ax.figure.savefig('{}/TF_eff.png'.format(args.output_folder), dpi=300, bbox_inches="tight")
    ax.figure.savefig('{}/TF_eff.pdf'.format(args.output_folder), transparent=True, bbox_inches="tight")

    f = up.open(os.path.join(args.dir, args.fit))
    region = 'prefit'
    fail_qcd = f['shapes_{}/{}/bkg'.format(region, 'fail')].values
    pass_qcd = f['shapes_{}/{}/bkg'.format(region, 'pass')].values

    fail_qcd = np.array(fail_qcd)
    pass_qcd = np.array(pass_qcd)
    q = np.sum(pass_qcd) / np.sum(fail_qcd)
    in_data_rat = (pass_qcd / (fail_qcd * q))

    ax = plotTF_ratio(in_data_rat, region=region)
    ax.figure.savefig('{}/{}{}.png'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", dpi=300)
    ax.figure.savefig('{}/{}{}.pdf'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", transparent=True)

    region = 'postfit'
    fail_qcd = f['shapes_{}/{}/bkg'.format('fit_s', 'fail')].values
    pass_qcd = f['shapes_{}/{}/bkg'.format('fit_s', 'pass')].values

    fail_qcd = np.array(fail_qcd)
    pass_qcd = np.array(pass_qcd)
    q = np.sum(pass_qcd) / np.sum(fail_qcd)
    in_data_rat = (pass_qcd / (fail_qcd * q))

    ax = plotTF_ratio(in_data_rat, region=region)
    ax.figure.savefig('{}/{}{}.png'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", dpi=300)
    ax.figure.savefig('{}/{}{}.pdf'.format(args.output_folder, "TF_ratio_", region),
                      bbox_inches="tight", transparent=True)
