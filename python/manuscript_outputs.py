"""Reproduce manuscript figures, interval summaries, and pointwise envelopes."""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cse import Bank, Estimator, empirical, features, generate

ROOT = Path(__file__).resolve().parents[1]
REPORTED = ROOT / "results" / "reported"


def interval_summary():
    data = pd.read_csv(ROOT / "results" / "inference_replicates.csv")
    columns = [f"{kind}_{par}_{stat}" for kind in ["basic", "percentile"]
               for par in ["phi", "v"] for stat in ["coverage", "length"]]
    columns += ["global_reject", "real_reject", "imaginary_reject", "boundary"]
    out = data.groupby(["n", "phi", "v"], sort=True)[columns].mean().reset_index()
    out.to_csv(ROOT / "results" / "inference_table.csv", index=False)
    return out


def reference_envelopes(draws=499, seed=20261104):
    point = pd.read_csv(REPORTED / "SP500_point_estimate.csv").iloc[0]
    theta = point[["phi", "v"]].to_numpy(dtype=float)
    n = 2262
    spec = features(n, diagnostic=True)
    ix = np.flatnonzero((spec.pairs[:, 0] == 0) &
                        (spec.pairs[:, 1] == len(spec.tau) - 1))[0]
    simulations = generate(theta, n, np.random.default_rng(seed), draws)
    values = empirical(simulations, spec)
    curves = values[:, ix * len(spec.lam):(ix + 1) * len(spec.lam)]
    output = ROOT / "results" / "reference_curves"
    output.mkdir(parents=True, exist_ok=True)
    for part, component in [("real", curves.real), ("imaginary", curves.imag)]:
        lower, upper = np.quantile(component, [.025, .975], axis=0, method="linear")
        path = REPORTED / f"SP500_{part}_curve.csv"
        frame = pd.read_csv(path)
        if not np.allclose(frame["lambda"], spec.lam, atol=1e-12, rtol=0):
            raise ValueError("The curve endpoints do not match the evaluation grid.")
        frame["lower"], frame["upper"] = lower, upper
        frame.to_csv(output / path.name, index=False, float_format="%.17g")
    pd.DataFrame({
        "draw": np.repeat(np.arange(1, draws + 1), len(spec.lam)),
        "lambda": np.tile(spec.lam, draws),
        "real": curves.real.ravel(), "imaginary": curves.imag.ravel()
    }).to_csv(ROOT / "results" / "SP500_reference_draws.csv", index=False)
    settings = {"n": n, "phi": float(theta[0]), "v": float(theta[1]),
                "draws": draws, "seed": seed, "stationary_initialization": True,
                "quantile_pair": [.1, .9], "endpoint_count": len(spec.lam),
                "probabilities": [.025, .975], "quantile_method": "linear",
                "target": "empirical centered integrated spectrum",
                "parameter_refitting": False, "simultaneous": False}
    (ROOT / "results" / "reference_envelope_settings.json").write_text(
        json.dumps(settings, indent=2) + "\n", encoding="utf-8")


def style():
    plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10,
                         "mathtext.fontset": "dejavuserif", "axes.titlesize": 10,
                         "axes.labelsize": 10, "legend.fontsize": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})


def identification_figure(directory):
    spec = features(1024)
    bank = Bank(spec.tau, upper_phi=.98)
    est = Estimator(spec, bank, upper_phi=.98)
    phi = np.linspace(.005, .98, 121)
    variance = np.linspace(.01, 3, 121)
    targets = [est.model(np.array(theta))[0] for theta in [(.6, 1), (.9, 1)]]
    surfaces = [np.empty((len(variance), len(phi))) for _ in targets]
    for row, v in enumerate(variance):
        for col, p in enumerate(phi):
            model = est.model(np.array([p, v]))[0]
            for out, target in zip(surfaces, targets):
                out[row, col] = np.dot(spec.weights, (model-target)**2)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.05), layout="constrained")
    for i, (ax, values, theta) in enumerate(zip(axs, surfaces, [(.6, 1), (.9, 1)])):
        im = ax.pcolormesh(phi, variance, np.log10(np.maximum(values, 1e-15)),
                           shading="auto", cmap="viridis", vmin=-10, vmax=-4,
                           rasterized=True)
        ax.plot(theta[0], theta[1], marker="x", color="white", ms=7, mew=1.6)
        ax.set(xlabel=r"Persistence $\phi$", ylabel=r"Latent variance $v$")
        ax.set_title(r"$\boldsymbol{\theta}_0=(%.1f,1)^\top$" % theta[0])
        ax.text(0, 1.03, f"({chr(97+i)})", transform=ax.transAxes, va="bottom")
    fig.colorbar(im, ax=axs, label=r"$\log_{10}Q^{\mathrm{grid}}(\boldsymbol{\theta})$",
                 ticks=np.arange(-10, -3))
    fig.savefig(directory / "Fig1.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)


def baseline_figures(directory):
    data = pd.read_csv(REPORTED / "simulation_baseline.csv")
    methods = ["CSE", "TD2", "TD16", "QML"]
    colors = ["#174e73", "#bb512e", "#43835b", "#8c5b99"]
    styles = ["-", "--", "-.", ":"]
    for number, parameter in [(2, "phi"), (3, "v")]:
        fig, axs = plt.subplots(3, 3, figsize=(7.6, 8.1), sharex=True)
        for i, phi in enumerate([.3, .6, .9]):
            for j, variance in enumerate([.25, 1., 2.]):
                ax = axs[i, j]
                part = data[(data.phi == phi) & (data.v == variance)]
                for method, color, ls in zip(methods, colors, styles):
                    name = "LSV_QML" if method == "QML" and "LSV_QML" in part.method.values else method
                    z = part[part.method == name].sort_values("n")
                    ax.plot(z.n, z[f"rmse_{parameter}"], label=method,
                            color=color, ls=ls, marker="o", ms=3, lw=1.3)
                ax.set_title(r"$\phi=%.1f,\ v=%g$" % (phi, variance))
                ax.text(0, 1.04, f"({chr(97+3*i+j)})", transform=ax.transAxes)
                ax.set_xscale("log", base=2)
                ax.set_xticks([256, 512, 1024, 2048], labels=["256", "512", "1024", "2048"])
                ax.set_ylim(bottom=0)
                ax.grid(alpha=.18)
                if j == 0:
                    ax.set_ylabel(r"RMSE($\phi$)" if parameter == "phi" else r"RMSE($v$)")
                if i == 2:
                    ax.set_xlabel("Sample size")
        handles, labels = axs[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False)
        fig.tight_layout(rect=(0, .045, 1, 1))
        fig.savefig(directory / f"Fig{number}.pdf", bbox_inches="tight")
        plt.close(fig)


def application_figure(directory, filename="Fig4.pdf"):
    fig, axs = plt.subplots(2, 2, figsize=(8, 5.45), sharex=True)
    for row, (app, title) in enumerate([("EURUSD", "EURUSD"), ("SP500", "S&P 500")]):
        for col, part in enumerate(["real", "imaginary"]):
            directory_curves = ROOT / "results" / "reference_curves" if app == "SP500" else REPORTED
            d = pd.read_csv(directory_curves / f"{app}_{part}_curve.csv")
            ax = axs[row, col]
            x = d["lambda"] / np.pi
            if app == "SP500":
                ax.fill_between(x, 1000*d.lower, 1000*d.upper, color="#a9c6da", alpha=.55,
                                label="95% pointwise reference envelope")
            ax.plot(x, 1000*d.observed, color="#174e73", lw=1.6, label="Empirical spectrum")
            ax.plot(x, 1000*d.model, color="#c35732", ls="--", lw=1.5, label="Fitted SV spectrum")
            ax.axhline(0, color=".65", lw=.5, zorder=0)
            ax.set_title(f"{title}: {part} component")
            label = r"$10^3\operatorname{Re}\mathcal{K}(\lambda;0.1,0.9)$" if col == 0 else r"$10^3\operatorname{Im}\mathcal{K}(\lambda;0.1,0.9)$"
            ax.set_ylabel(label)
            ax.text(0, 1.03, f"({chr(97+2*row+col)})", transform=ax.transAxes, va="bottom")
            ax.set_xlim(0, 1)
            ax.grid(alpha=.15)
            if row == 1:
                ax.set_xlabel(r"Integration endpoint $\lambda/\pi$")
            ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axs[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(.5, .003), fontsize=8)
    fig.tight_layout(rect=(0, .065, 1, 1))
    fig.savefig(directory / filename, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("task", choices=["intervals", "envelopes", "figures", "all"])
    parser.add_argument("--output", type=Path, default=ROOT / "figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.task in ["intervals", "all"]:
        interval_summary()
    if args.task in ["envelopes", "figures", "all"]:
        reference_envelopes()
    if args.task in ["figures", "all"]:
        style()
        identification_figure(args.output)
        baseline_figures(args.output)
        application_figure(args.output)


if __name__ == "__main__":
    main()
