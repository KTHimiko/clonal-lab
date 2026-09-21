#!/usr/bin/env python3
"""
Stage C — score the parameter sweep against the observed cohort.

The cluster produced summary statistics for every point of the (s, mu) grid.
This script does the cheap half: compare each point with the observed data
and keep the ones that land close. The distribution of the kept parameters is
the approximate posterior.

Deliberately separate from the sweep. Re-scoring against another cohort or
another distance function costs seconds and re-runs nothing.

Usage:  .venv/bin/python analysis/04_abc.py [path/to/sweep.csv]
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
SWEEP = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sweep.csv"
CSV = (ROOT / "reference/data/watson2020/Maximum_likelihood_estimations"
            / "Maximum likelihood estimations - data files"
            / "all_studies_trimmed_all_genes.csv")
FIGS = ROOT / "analysis/figures"

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def clean(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

def section(t):
    print(f"\n{'-'*74}\n{t}\n{'-'*74}")


COHORT = "Coombs2017"
LIMIT = 0.0192
COMPARE_AGE = 70.0
QCOLS = ["vaf_q10", "vaf_q25", "vaf_q50", "vaf_q75", "vaf_q90"]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
ACCEPT = 0.05                      # keep the closest 5% of the grid


# ================================================================ observed
section("1. THE OBSERVED TARGET")
raw = pd.read_csv(CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
raw["age"] = pd.to_numeric(raw["age"], errors="coerce")
raw["VAF"] = pd.to_numeric(raw["VAF"], errors="coerce")
obs = raw.dropna(subset=["age", "VAF"])
obs = obs[(obs.study == COHORT) & (obs.VAF >= LIMIT)]

obs_q = np.quantile(obs.VAF.values, QUANTILES)
print(f"cohort {COHORT}: {len(obs)} variants above VAF {LIMIT}")
print(f"median age {obs.age.median():.0f}; comparing against simulated age {COMPARE_AGE:.0f}")
print("\nobserved quantiles:")
for q, v in zip(QUANTILES, obs_q):
    print(f"  q{int(q*100):02d}  {v:.4f}")


# ================================================================ the grid
section("2. THE SWEEP")
if not SWEEP.exists():
    print(f"sweep not found at {SWEEP}")
    sys.exit(1)

sw = pd.read_csv(SWEEP)
sw = sw[sw.age == COMPARE_AGE].copy()
print(f"grid points: {sw[['s','mu']].drop_duplicates().shape[0]}")
print(f"rows at age {COMPARE_AGE:.0f}: {len(sw)} (includes replicates)")
print(f"s  from {sw.s.min():.3f} to {sw.s.max():.3f}")
print(f"mu from {sw.mu.min():.2e} to {sw.mu.max():.2e}")

usable = sw.dropna(subset=QCOLS)
print(f"points with detectable clones: {len(usable)} of {len(sw)}")
if len(usable) < len(sw):
    print("  (points producing no clone above the limit cannot be scored —")
    print("   they are rejected, which is itself information)")


# =============================================================== the distance
section("3. DISTANCE AND ACCEPTANCE")
print("Quantiles are compared in log space, because VAF spans three orders of")
print("magnitude and an absolute difference would be dominated by the tail.")
print("Each quantile is scaled by its spread across the grid so that no single")
print("one dominates the distance.\n")

sim = np.log(usable[QCOLS].values)
target = np.log(obs_q)
scale = sim.std(axis=0)
scale[scale == 0] = 1.0
usable = usable.assign(distance=np.sqrt((((sim - target) / scale) ** 2).mean(axis=1)))

# average the replicates: the same parameters should be scored once
grid = (usable.groupby(["s", "mu"], as_index=False)
              .agg(distance=("distance", "mean"),
                   clones_per_person=("clones_per_person", "mean"),
                   carriers=("carriers", "mean"),
                   vaf_q50=("vaf_q50", "mean")))

cut = np.quantile(grid.distance, ACCEPT)
accepted = grid[grid.distance <= cut]
print(f"accepted the closest {ACCEPT:.0%}: {len(accepted)} of {len(grid)} grid points")
print(f"distance threshold: {cut:.4f}\n")

best = grid.loc[grid.distance.idxmin()]
print(f"closest point:  s = {best.s:.3f}   mu = {best.mu:.3e}   distance {best.distance:.4f}")
print(f"  its median VAF {best.vaf_q50:.4f} against observed {obs_q[2]:.4f}")

print(f"\nposterior for s:   {accepted.s.min():.3f} to {accepted.s.max():.3f}"
      f"   (median {accepted.s.median():.3f})")
print(f"posterior for mu:  {accepted.mu.min():.2e} to {accepted.mu.max():.2e}"
      f"   (median {accepted.mu.median():.2e})")

s_span = accepted.s.max() - accepted.s.min()
s_prior = grid.s.max() - grid.s.min()
mu_span = np.log10(accepted.mu.max()) - np.log10(accepted.mu.min())
mu_prior = np.log10(grid.mu.max()) - np.log10(grid.mu.min())
print(f"\nshrinkage against the grid:")
print(f"  s:   posterior spans {100*s_span/s_prior:.0f}% of the grid")
print(f"  mu:  posterior spans {100*mu_span/mu_prior:.0f}% of the grid (log scale)")


# =================================================================== figures
section("4. FIGURES")

fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.1),
                         gridspec_kw={"width_ratios": [1.35, 1, 1]})

# left: the distance surface over the grid
ax = axes[0]
piv = grid.pivot(index="mu", columns="s", values="distance")
im = ax.pcolormesh(piv.columns.values, piv.index.values, piv.values,
                   cmap=matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ[::-1]),
                   shading="nearest")
# Marked directly rather than with a legend: a legend box sits inside the
# plotting area here, and its sample marker reads as a second data point.
ax.scatter([best.s], [best.mu], s=90, facecolor="none", edgecolor=INK,
           linewidth=2, zorder=3)
ax.annotate(f"closest\ns={best.s:.2f}", xy=(best.s, best.mu),
            xytext=(10, 14), textcoords="offset points",
            color=INK, fontsize=9, fontweight="bold")
ax.set_yscale("log")
ax.set_xlabel("selection coefficient  s")
ax.set_ylabel("mutation rate  mu  (per cell per year)")
ax.set_title("Distance to the observed distribution", fontsize=11)
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label("distance", color=INK2)
cb.ax.tick_params(colors=MUTED)
ax.grid(False)
clean(ax)

# middle and right: the marginals
for ax, col, label, logx in ((axes[1], "s", "selection coefficient  s", False),
                             (axes[2], "mu", "mutation rate  mu", True)):
    vals = np.sort(grid[col].unique())
    kept = accepted[col].values
    counts = [np.sum(np.isclose(kept, v)) for v in vals]
    width = (np.diff(np.log10(vals)).mean() if logx else np.diff(vals).mean())
    if logx:
        ax.bar(vals, counts, width=[v*(10**(width/2)-10**(-width/2)) for v in vals],
               color=S1, edgecolor=SURF, linewidth=1)
        ax.set_xscale("log")
    else:
        ax.bar(vals, counts, width=width*0.85, color=S1, edgecolor=SURF, linewidth=1)
    ax.set_xlabel(label)
    ax.set_ylabel("accepted points")
    ax.set_title(f"Posterior for {col}", fontsize=11)
    clean(ax)

fig.suptitle("Stage C: approximate Bayesian computation over the sweep",
             x=0.008, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(FIGS/"08_abc_posterior.png"); plt.close(fig)
print("  08_abc_posterior.png")


# =================================================================== reading
section("5. READING THE RESULT")
print(f"""The fitness coefficient is constrained; the mutation rate is not, and the
reason is the limitation found back in stage 1.

The observed file lists only the variants that were DETECTED. It carries no
denominator — how many people were screened at each age — so the fraction of
people carrying a clone cannot be computed from it. That leaves only the
SHAPE of the VAF distribution to compare against.

Shape responds to fitness: a larger s grows clones faster and pushes the
distribution right. Shape barely responds to the mutation rate, which changes
how MANY clones exist rather than how big they are. With counts unavailable,
mu is close to unidentifiable from this dataset.

This is a real result, not a failure of the method. It says precisely what
extra data would be worth having: a cohort reporting how many people were
screened, not just which variants turned up.""")
