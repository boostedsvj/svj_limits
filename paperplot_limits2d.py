import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse
import ROOT
import glob
import mplhep as hep

matplotlib.style.use(hep.style.CMS)


def construct_limits_table(base_dir: str):
    """
    Constructing all limits as a dictionary of numpy arrays
    """
    return_dict = {
        "limit": [],
        "quantile": [],
        "mZprime": [],
        "mDark": [],
        "rinv": [],
        "xsec": [],
    }
    for file in glob.glob(base_dir + "/*.root"):
        try:
            f = ROOT.TFile.Open(file, "READ")
            ltree = f.Get("limit")
            for entry in ltree:
                return_dict["limit"].append(entry.limit)
                return_dict["quantile"].append(entry.quantileExpected)
                return_dict["mZprime"].append(entry.trackedParam_mZprime)
                return_dict["mDark"].append(entry.trackedParam_mDark)
                return_dict["rinv"].append(entry.trackedParam_rinv)
                return_dict["xsec"].append(entry.trackedParam_xsec)
        except:
            pass
    return {k: np.array(v) for k, v in return_dict.items()}


def make_interpolator(table, mDark: float, quantile=None, observed=None, xsec=False):
    """Returning the 2 interpolation to allow for arbitrary mesh grid"""
    assert not (quantile is None and observed is None), (
        "Must specifiy either quantile or observed"
    )
    assert not (quantile is not None and observed is not None), (
        "Cannot specify both quantile and observed"
    )

    x = table["mZprime"]
    y = table["rinv"]
    z = table["limit"]
    if xsec == True:
        z = z * table["xsec"]
    if observed is True:
        quantile = -1

    filter = (table["quantile"] == quantile) & (table["mDark"] == mDark)
    (x, y, z) = x[filter], y[filter], z[filter]
    return matplotlib.tri.LinearTriInterpolator(matplotlib.tri.Triangulation(x, y), z)


def _make_mesh(n_entries):
    """Slightly shifting the lower boundary for mMed to avoid axis tick clashing"""
    return np.meshgrid(np.linspace(201, 550, n_entries), np.linspace(0, 1, n_entries))


def plot_2d_color(ax, interp, n_entries=50):
    """Plotting the 2d clolor map"""
    x, y = _make_mesh(n_entries)
    im = ax.pcolormesh(
        *(x, y, interp(x, y)),
        norm=matplotlib.colors.LogNorm(vmin=3e-1, vmax=3e2),
        linewidth=0.0,
        edgecolors="None",
    )
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("95% CL upper limit on cross section [fb]", va="top")
    return ax


def plot_limit_contour(ax, interp, n_entries=50, color="red", **kwargs):
    x, y = _make_mesh(n_entries)
    ax.contour(*(x, y, interp(x, y)), levels=np.array([1]), colors=color, **kwargs)
    return ax.plot([], [], color=color, **kwargs)


def plot_limit_band(ax, interp1, interp2, n_entries=2000, color="red", **kwargs):
    x, y = _make_mesh(n_entries)
    z = (interp1(x, y) < 1) ^ (interp2(x, y) < 1)
    ax.contourf(
        *(x, y, z),
        levels=np.array([0.5, 2]),
        cmap=matplotlib.colors.ListedColormap([color]),
        alpha=0.2,
    )
    return ax.fill(np.NaN, np.NaN, color=color, alpha=0.2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Script for plotting limits ready for publication")
    parser.add_argument(
        "--base_dir",
        "-d",
        type=str,
        default="Limits",
        help="Directory containig the root files that start the limit information files here should be generated by the cls_maker.py script",
    )
    parser.add_argument(
        "--mDark",
        "-m",
        type=float,
        default=10,
        choices=[10],
        help="Dark meson mass point to use for plotting",
    )
    parser.add_argument(
        "--label", type=str, default="Preliminary", help="Plot label used for display"
    )
    parser.add_argument(
        "--observed",
        "-o",
        type=str,
        default="False",
        choices=["True", "False", "Dummy"],
        help="Whether or not to include the observed limit",
    )
    args = parser.parse_args()

    table = construct_limits_table(args.base_dir)

    fig = plt.figure(constrained_layout=True, figsize=(11, 11))
    spec = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    ax = fig.add_subplot(spec[0, 0])
    ax.set_xlabel("", horizontalalignment="right", x=1.0)
    ax.set_ylabel("", horizontalalignment="right", y=1.0)

    legend_entries = {}

    exp_central_fb = make_interpolator(table, mDark=args.mDark, quantile=0.5, xsec=True)
    plot_2d_color(ax, exp_central_fb)

    exp_central = make_interpolator(table, mDark=args.mDark, quantile=0.5)
    exp_up = make_interpolator(table, mDark=args.mDark, quantile=0.16)
    exp_lo = make_interpolator(table, mDark=args.mDark, quantile=0.84)
    p1 = plot_limit_contour(ax, exp_central, color="red")
    p2 = plot_limit_band(ax, exp_up, exp_lo, color="red")
    legend_entries[r"Exp. limit $\pm 1\sigma$"] = (p1[0], p2[0])

    if args.observed != "False":
        obs = make_interpolator(
            table, mDark=args.mDark, quantile=-1 if args.observed == "True" else 0.16
        )
        po = plot_limit_contour(ax, obs, color="black")
        legend_entries[
            "Obs. limit (dummy)" if args.observed == "Dummy" else "Obs. limit"
        ] = po[0]

    hep.cms.text(text=args.label, ax=ax, loc=0)
    hep.cms.lumitext(text="138 $fb^{-1} (13 TeV)$")

    ax.legend(
        list(legend_entries.values()),
        list(legend_entries.keys()),
        title="$m_{Dark}$ = 10 GeV",
        loc="upper right",
        frameon=True,
    )

    ax.set_xlabel("$m_{Z'}$ [GeV]")
    ax.set_ylabel("$r_{inv}$")
    fig.savefig(f"limits2d_{args.label.replace(' ', '-')}_mDark-{args.mDark}.pdf")
