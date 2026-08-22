#!/usr/bin/env python3
"""
Benchmark: GLUE-anchored metacell pairing vs SCENIC+'s stock label-random pairing.

Reproduces `pairing_sensitivity.csv` and `pairing_benchmark.png` (both written
next to this script).

DESIGN
------
Ground truth is a SHARED latent geometry -- cell-type centroids plus N_PROG
independent within-cell-type programs -- which is what a GLUE co-embedding
approximates. Peak i and gene i are both driven by program i, so (peak i, gene i)
is a true link and (peak i, gene j) is a genuine decoy. Modality-specific noise
added to the latent coordinates stands in for integration error, and sweeping it
maps the operating envelope.

Two pairing methods are compared on identical data:
  * stock  -- verbatim reimplementation of generate_pseudocells_for_numpy:
              uniform random draw within a label, independently per modality,
              with random.seed(x) inside the loop (the upstream behaviour).
  * GLUE   -- build_metacells_for_group from 02_pair/glue_metacells.py.

Usage:  python benchmark_pairing.py [--quick]
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "02_pair"))
from glue_metacells import aggregate, build_metacells_for_group, l2_normalize  # noqa: E402

CT = ["HSC", "MPP", "GMP", "MEP", "Mono", "Ery"]
N_PROG, N_NOISE, DIM = 40, 200, 30
SEED = 7

_rng = np.random.default_rng(SEED)
BASE = {c: _rng.normal(0, 1, DIM) * 2.0 for c in CT}
PROG = _rng.normal(0, 1, (N_PROG, DIM))
PROG /= np.linalg.norm(PROG, axis=1, keepdims=True)


def build(noise: float, n_rna=5000, n_atac=4000, seed=11):
    """Two unpaired modalities sharing one latent geometry, plus per-modality noise."""
    r = np.random.default_rng(seed)

    def sample(n):
        lab = r.choice(CT, n)
        act = r.uniform(0, 1, (n, N_PROG))          # per-cell program activities
        Z = np.vstack([BASE[l] for l in lab]) + act @ PROG
        Z += r.normal(0, noise, (n, DIM))           # modality-specific error
        return lab, act, Z.astype(np.float32)

    def mat(act):
        M = r.normal(0, 1, (len(act), N_PROG + N_NOISE))
        M[:, :N_PROG] += 3.0 * act                  # feature i driven by program i
        return M.astype(np.float32)

    lr, ar, Zr = sample(n_rna)
    la, aa, Za = sample(n_atac)
    rna = ad.AnnData(
        X=mat(ar),
        obs=pd.DataFrame({"celltype": pd.Categorical(lr)},
                         index=[f"r{i}" for i in range(n_rna)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(N_PROG + N_NOISE)]))
    rna.obsm["X_glue"] = Zr
    atac = ad.AnnData(
        X=mat(aa),
        obs=pd.DataFrame({"celltype": pd.Categorical(la)},
                         index=[f"a{i}" for i in range(n_atac)]),
        var=pd.DataFrame(index=[f"chr1:{1000 + i * 1000}-{1501 + i * 1000}"
                                for i in range(N_PROG + N_NOISE)]))
    atac.obsm["X_glue"] = Za
    return rna, atac


def stock_pair(rna, atac, key="celltype", k=25):
    """Verbatim SCENIC+ behaviour: label-random, independent per modality."""
    R, A = [], []
    for g in sorted(set(rna.obs[key].astype(str)) & set(atac.obs[key].astype(str))):
        ri = np.flatnonzero(rna.obs[key].astype(str) == g)
        ai = np.flatnonzero(atac.obs[key].astype(str) == g)
        for x in range(max(1, round(min(len(ri), len(ai)) / k) * 2)):
            random.seed(x)                          # seed inside the loop, as upstream
            R.append(rna.X[np.sort(random.sample(list(ri), k))].mean(0))
            A.append(atac.X[np.sort(random.sample(list(ai), k))].mean(0))
    return np.vstack(R), np.vstack(A)


def glue_pair(rna, atac, key="celltype", k=25):
    Zr = l2_normalize(rna.obsm["X_glue"])
    Za = l2_normalize(atac.obsm["X_glue"])
    R, A = [], []
    for g in sorted(set(rna.obs[key].astype(str)) & set(atac.obs[key].astype(str))):
        ri = np.flatnonzero(rna.obs[key].astype(str) == g)
        ai = np.flatnonzero(atac.obs[key].astype(str) == g)
        n_mc = max(1, round(min(len(ri), len(ai)) / k) * 2)
        lr, la, _ = build_metacells_for_group(
            Zr[ri], Za[ai], n_mc, k, k, np.random.default_rng(0), "both")
        R.append(aggregate(rna.X, [ri[x] for x in lr]))
        A.append(aggregate(atac.X, [ai[x] for x in la]))
    return np.vstack(R), np.vstack(A)


def evaluate(E, Acc):
    """Separation of true (peak i, gene i) links from decoy (peak i, gene j) pairs."""
    true_rho = np.array([spearmanr(Acc[:, i], E[:, i]).statistic for i in range(N_PROG)])
    off = [(i, j) for i in range(N_PROG) for j in range(N_PROG) if i != j][:2000]
    decoy = np.array([spearmanr(Acc[:, i], E[:, j]).statistic for i, j in off])
    auroc = (mannwhitneyu(true_rho, decoy, alternative="greater").statistic
             / (len(true_rho) * len(decoy)))
    return dict(n_metacells=int(E.shape[0]),
                median_rho_true=round(float(np.median(true_rho)), 3),
                median_rho_decoy=round(float(np.median(decoy)), 3),
                frac_true_gt_0p3=round(float((true_rho > 0.3).mean()), 3),
                auroc=round(float(auroc), 3))


def run_sweep(noises, ks):
    rows = []
    for noise in noises:
        rna, atac = build(noise)
        for k in ks:
            E, A = glue_pair(rna, atac, k=k)
            rows.append({**evaluate(E, A), "noise": noise, "k": k, "pair": "GLUE"})
        E, A = stock_pair(rna, atac, k=25)
        rows.append({**evaluate(E, A), "noise": noise, "k": 25, "pair": "stock"})
        print(f"  noise={noise}: done", flush=True)
    return pd.DataFrame(rows)


def make_figure(sw, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    META_GREY = "#8a8a8a"
    C_GLUE = {10: "#9ecae1", 25: "#3b7dd8", 50: "#08306b"}
    C_STOCK = "#c44e52"
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6})

    g, st = sw[sw.pair == "GLUE"], sw[sw.pair == "stock"]
    xmax = sw.noise.max() * 1.10
    XT = [t for t in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0) if t <= xmax]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.1))

    H, L = [], []
    for k in sorted(g.k.unique()):
        d = g[g.k == k].sort_values("noise")
        ln, = axA.plot(d.noise, d.auroc, "-o", ms=3.6, lw=1.5,
                       color=C_GLUE.get(k, "#3b7dd8"), zorder=3)
        H.append(ln); L.append(f"GLUE-anchored, k={k}")
    d = st.sort_values("noise")
    ln, = axA.plot(d.noise, d.auroc, "-s", ms=3.6, lw=1.5, color=C_STOCK, zorder=3)
    H.append(ln); L.append("SCENIC+ stock, k=25")
    axA.axhline(0.5, ls=":", lw=1.0, color=META_GREY, zorder=1)
    axA.set_ylabel("AUROC, true vs. decoy peak\u2013gene links")
    axA.set_title("Recovery of true peak\u2013gene links", loc="left")
    axA.set_xlim(0, xmax); axA.set_ylim(0.40, 1.10); axA.set_xticks(XT)
    axA.text(xmax * 0.98, 0.487, "chance", fontsize=6, color=META_GREY,
             va="top", ha="right")
    axA.legend(H, L, frameon=False, fontsize=6, loc="upper right",
               bbox_to_anchor=(1.005, 0.62), handlelength=1.6, labelspacing=0.35)

    for k in sorted(g.k.unique()):
        d = g[g.k == k].sort_values("noise")
        axB.plot(d.noise, d.median_rho_true, "-o", ms=3.6, lw=1.5,
                 color=C_GLUE.get(k, "#3b7dd8"), zorder=3)
    d = st.sort_values("noise")
    axB.plot(d.noise, d.median_rho_true, "-s", ms=3.6, lw=1.5, color=C_STOCK, zorder=3)
    axB.axhline(0.0, ls=":", lw=1.0, color=META_GREY, zorder=1)
    axB.set_ylabel("Median Spearman \u03c1 on true links")
    axB.set_title("Effect size reaching the region\u2013gene model", loc="left")
    axB.set_xlim(0, xmax); axB.set_ylim(-0.055, 0.66); axB.set_xticks(XT)
    axB.text(xmax * 0.60, 0.255,
             "stock pairing stays at \u03c1\u22480\nat every integration quality",
             fontsize=6, color=C_STOCK, ha="center", va="bottom")

    for ax, lab in ((axA, "a"), (axB, "b")):
        ax.set_xlabel("Modality-specific noise in shared embedding\n"
                      "(integration quality, worse \u2192)")
        ax.text(-0.16, 1.06, lab, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="left")
        ax.text(0.02, 0.98, "\u2191 higher = better", transform=ax.transAxes,
                fontsize=6, color=META_GREY, va="top", ha="left")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"wrote {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="3 noise levels, k=25 only")
    args = ap.parse_args()
    noises = [0.30, 1.00, 2.50] if args.quick else [0.15, 0.30, 0.60, 1.00, 1.50, 2.50]
    ks = [25] if args.quick else [10, 25, 50]

    print("running sweep ...")
    sw = run_sweep(noises, ks)
    csv = HERE / "pairing_sensitivity.csv"
    sw.to_csv(csv, index=False)
    print(f"wrote {csv}")
    make_figure(sw, HERE / "pairing_benchmark.png")

    head = sw[(sw.noise == noises[0])]
    print("\nAt best integration quality tested (noise="
          f"{noises[0]}):")
    for _, r in head.iterrows():
        print(f"  {r.pair:6s} k={int(r.k):2d}  rho_true={r.median_rho_true:.3f}  "
              f"rho_decoy={r.median_rho_decoy:.3f}  AUROC={r.auroc:.3f}")


if __name__ == "__main__":
    main()
