#!/usr/bin/env python3
"""
Does OVERSAMPLING (more metacells at fixed k) buy statistical power, or
pseudoreplication?

Reproduces `oversample_sweep.csv`, `oversample_null.csv` and
`oversampling_tradeoff.png` (all written next to this script).

WHY THIS EXISTS
---------------
Metacell SIZE (k) and metacell COUNT are independent knobs, because anchors may
overlap: you can hold k=50 -- the size with the largest effect size in
`benchmark_pairing.py` -- and still emit many more metacells. Whether the extra
metacells carry information or merely re-weight the same cells is an empirical
question, so it is measured here rather than asserted.

DESIGN
------
Shared latent geometry, as in benchmark_pairing.py: cell-type centroids plus
N_PROG independent within-cell-type programs. Peak i and gene i are both driven
by program i, so (peak i, gene i) is a true link and (peak i, gene j) is a
genuine decoy. Modality-specific noise added to the latent coordinates stands
in for integration error.

TWO SWEEPS, and they answer different questions:

  * oversample_sweep.csv -- RANKING. AUROC separating true links from decoys, and
    median rho on true links, over noise x k x oversample. The first version of
    this sweep used only the easy regime (low noise) and saturated at AUROC
    ~1.0, which is uninformative; the noise levels here are the harder regimes
    where the question has teeth.

  * oversample_null.csv -- CALIBRATION. Under a PURE NULL (no true links at
    all), how often does a naive p-value that treats n_metacells as the sample
    size clear 0.05 and 0.001? Overlapping metacells are not independent
    observations, so this is where pseudoreplication would show up.

The conclusion those two together support: oversample freely to improve
RANKING, never to compute significance. `aggregate_atac_sparse.py` records
`independent_metacell_equiv` for exactly that reason -- it is the honest
denominator.

Usage
-----
    python docs/benchmark_oversample.py             # full sweep
    python docs/benchmark_oversample.py --quick     # fewer points
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

HERE = Path(__file__).resolve().parent
CT = ["HSC", "MPP", "GMP", "MEP", "Mono", "Ery"]
N_PROG, N_NOISE, DIM = 40, 200, 30

# Shared latent geometry: the same centroids and program directions for both
# modalities, which is what a GLUE co-embedding approximates.
_RG = np.random.default_rng(0)
BASE = {c: _RG.normal(0, 1, DIM) * 2.0 for c in CT}
PROG = _RG.normal(0, 1, (N_PROG, DIM))
PROG /= np.linalg.norm(PROG, axis=1, keepdims=True)


def build(n_rna, n_atac, noise, seed, linked=True):
    """Two unpaired modalities on the shared geometry.

    linked=False severs the peak-gene correspondence, giving a pure null.
    """
    r = np.random.default_rng(seed)

    def side(n):
        lab = r.choice(CT, n)
        act = r.uniform(0, 1, (n, N_PROG))
        Z = np.vstack([BASE[l] for l in lab]) + act @ PROG \
            + r.normal(0, noise, (n, DIM))
        return lab, act, Z.astype(np.float32)

    lr, ar, Zr = side(n_rna)
    la, aa, Za = side(n_atac)

    def mat(act, shuffle_cols):
        M = r.normal(0, 1, (len(act), N_PROG + N_NOISE))
        drive = act[:, r.permutation(N_PROG)] if shuffle_cols else act
        M[:, :N_PROG] += 3.0 * drive
        return M.astype(np.float32)

    # Under the null the ATAC programs are permuted relative to RNA, so peak i
    # and gene i are driven by different programs while every marginal
    # distribution is unchanged.
    return (lr, ar, Zr, mat(ar, False)), (la, aa, Za, mat(aa, not linked))


def l2(Z):
    n = np.linalg.norm(Z, axis=1, keepdims=True)
    return Z / np.where(n == 0, 1, n)


def pair(rna, atac, k, oversample, seed=0):
    """GLUE-anchored pairing, matching build_metacells_for_group: farthest-point
    anchors, then per-modality nearest neighbours. Anchors may overlap, which is
    what makes oversampling possible at fixed k."""
    lr, _, Zr, Er = rna
    la, _, Za, Aa = atac
    Zr, Za = l2(Zr), l2(Za)
    rng = np.random.default_rng(seed)
    E, A = [], []
    for c in CT:
        ir = np.flatnonzero(lr == c)
        ia = np.flatnonzero(la == c)
        if len(ir) < k or len(ia) < k:
            continue
        n_mc = max(1, int(min(len(ir), len(ia)) / k * oversample))
        Zrg, Zag = Zr[ir], Za[ia]
        # farthest-point anchor selection in the shared space
        anchors = [rng.integers(len(Zrg))]
        d = np.linalg.norm(Zrg - Zrg[anchors[0]], axis=1)
        while len(anchors) < n_mc:
            nxt = int(np.argmax(d))
            anchors.append(nxt)
            d = np.minimum(d, np.linalg.norm(Zrg - Zrg[nxt], axis=1))
        for a in anchors:
            z = Zrg[a]
            nr = ir[np.argsort(np.linalg.norm(Zrg - z, axis=1))[:k]]
            na = ia[np.argsort(np.linalg.norm(Zag - z, axis=1))[:k]]
            E.append(Er[nr].mean(0))
            A.append(Aa[na].mean(0))
    return np.vstack(E), np.vstack(A)


def score(E, Acc, n_decoy=2000):
    true_rho = np.array([spearmanr(Acc[:, i], E[:, i]).statistic
                         for i in range(N_PROG)])
    off = [(i, j) for i in range(N_PROG) for j in range(N_PROG) if i != j][:n_decoy]
    decoy = np.array([spearmanr(Acc[:, i], E[:, j]).statistic for i, j in off])
    ok_t, ok_d = np.isfinite(true_rho), np.isfinite(decoy)
    auroc = (mannwhitneyu(true_rho[ok_t], decoy[ok_d], alternative="greater").statistic
             / (ok_t.sum() * ok_d.sum()))
    return (float(np.median(true_rho[ok_t])), float(np.median(decoy[ok_d])),
            float(auroc))


def sweep_ranking(noises, ks, oversamples, n_rna=20_000, n_atac=8_000):
    rows = []
    for noise in noises:
        data = build(n_rna, n_atac, noise=noise, seed=5)
        for k in ks:
            for ov in oversamples:
                E, A = pair(*data, k=k, oversample=ov)
                rt, rd, au = score(E, A)
                rows.append(dict(noise=noise, k=k, oversample=ov,
                                 n_metacells=E.shape[0],
                                 rho_true=round(rt, 4), rho_decoy=round(rd, 4),
                                 auroc=round(au, 4)))
                print(f"  noise={noise} k={k} ov={ov:2d} -> "
                      f"n={E.shape[0]:5d} rho={rt:+.3f} auroc={au:.3f}", flush=True)
    return pd.DataFrame(rows)


def sweep_null(oversamples, k=50, noise=1.0, n_rep=8):
    """Pure null: no true peak-gene links. A naive p-value treating n_metacells
    as the sample size should reject at ~5%; if the rate climbs with
    oversampling, those metacells are not independent observations."""
    rows = []
    for ov in oversamples:
        p05, p001, ns = [], [], []
        for rep in range(n_rep):
            data = build(6_000, 4_000, noise=noise, seed=100 + rep, linked=False)
            E, A = pair(*data, k=k, oversample=ov, seed=rep)
            ps = []
            for i in range(N_PROG):
                r = spearmanr(A[:, i], E[:, i])
                if np.isfinite(r.pvalue):
                    ps.append(r.pvalue)
            ps = np.array(ps)
            p05.append((ps < 0.05).mean())
            p001.append((ps < 0.001).mean())
            ns.append(E.shape[0])
        rows.append(dict(oversample=ov, k=k, noise=noise,
                         mean_n_metacells=int(np.mean(ns)),
                         frac_p_lt_0p05=round(float(np.mean(p05)), 4),
                         frac_p_lt_0p001=round(float(np.mean(p001)), 4)))
        print(f"  null ov={ov:2d} -> n={int(np.mean(ns)):5d} "
              f"p<0.05: {np.mean(p05):.3f}  p<0.001: {np.mean(p001):.3f}", flush=True)
    return pd.DataFrame(rows)


def make_figure(sw, nl, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8,
                         "axes.labelsize": 7, "legend.fontsize": 6,
                         "xtick.labelsize": 6, "ytick.labelsize": 6})
    GREY = "#6b7075"
    C = {1.0: "#9ecae1", 1.5: "#3b7dd8", 2.5: "#08306b"}
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.4, 3.1))

    h = sw[sw.k == 50]
    for nz in sorted(h.noise.unique()):
        d = h[h.noise == nz].sort_values("oversample")
        axA.plot(d.oversample, d.auroc, "-o", ms=3.6, lw=1.5,
                 color=C.get(nz, GREY), label=f"noise {nz}")
    axA.set_xscale("log", base=2)
    axA.set_xticks(sorted(h.oversample.unique()))
    axA.set_xticklabels([str(x) for x in sorted(h.oversample.unique())])
    axA.axhline(0.5, ls=":", lw=1.0, color=GREY)
    axA.set_xlabel("oversample factor (metacells per cell type, k=50 fixed)")
    axA.set_ylabel("AUROC, true vs. decoy peak\u2013gene links")
    axA.set_title("Ranking improves with oversampling", loc="left")
    axA.legend(frameon=False, title="integration noise", title_fontsize=6)

    axB.plot(nl.oversample, nl.frac_p_lt_0p05, "-o", ms=3.6, lw=1.5,
             color="#c44e52", label="p < 0.05")
    axB.plot(nl.oversample, nl.frac_p_lt_0p001, "-s", ms=3.6, lw=1.5,
             color="#e8a33d", label="p < 0.001")
    axB.axhline(0.05, ls=":", lw=1.0, color=GREY)
    axB.set_xscale("log", base=2)
    axB.set_xticks(list(nl.oversample))
    axB.set_xticklabels([str(x) for x in nl.oversample])
    axB.set_xlabel("oversample factor (pure null: no true links)")
    axB.set_ylabel("fraction of null links called significant")
    axB.set_title("Naive significance inflates \u2014 do not use for p-values",
                  loc="left")
    axB.legend(frameon=False)
    # Label the reference line at the RIGHT edge, below it: at the left edge it
    # sits directly on the p<0.05 series, which starts at 0.056.
    axB.text(nl.oversample.max() * 0.97, 0.045, "nominal 0.05",
             fontsize=6, color=GREY, va="top", ha="right")

    for ax, L in ((axA, "a"), (axB, "b")):
        ax.text(-0.14, 1.06, L, transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    noises = [1.0, 2.5] if args.quick else [1.0, 1.5, 2.5]
    ks = [50] if args.quick else [25, 50]
    ovs = [1, 4, 16] if args.quick else [1, 2, 4, 8, 16]

    print("ranking sweep (harder regimes -- the easy regime saturates at AUROC 1.0)")
    sw = sweep_ranking(noises, ks, ovs)
    sw.to_csv(HERE / "oversample_sweep.csv", index=False)
    print(f"wrote {HERE / 'oversample_sweep.csv'}")

    print("\nnull calibration sweep")
    nl = sweep_null(ovs, n_rep=4 if args.quick else 8)
    nl.to_csv(HERE / "oversample_null.csv", index=False)
    print(f"wrote {HERE / 'oversample_null.csv'}")

    make_figure(sw, nl, HERE / "oversampling_tradeoff.png")

    k50 = sw[sw.k == 50]
    print("\nAUROC at k=50, by noise x oversample:")
    print(k50.pivot_table(index="noise", columns="oversample",
                          values="auroc").to_string())
    print("\nmedian rho on true links (flat across oversample = no effect-size "
          "inflation):")
    print(k50.pivot_table(index="noise", columns="oversample",
                          values="rho_true").to_string())


if __name__ == "__main__":
    main()
