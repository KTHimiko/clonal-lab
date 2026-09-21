#!/usr/bin/env python3
"""
Stage K — does a fitness distribution, with N free, fit both observables?

Decision rule and identifiability checks were committed before the sweep ran:
analysis/PREREGISTRATION_STAGE_K.md (42e75b3).

    a point "fits both" if it is within 25% of the best distance achieved on
    EACH arm separately. If such a point exists, report it and stop. If not,
    the model cannot fit both even with N free — name it and stop.

Before any fit is reported the plan requires two checks, which stage J's plan
forgot and paid for:
    1. is the spectrum arm monotonic in the mean, at fixed shape and N?
    2. is the joint minimum unique?

Usage:  .venv/bin/python analysis/15_mixture_fit.py [path/to/sweep_mixture.csv]
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
GRID = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sweep_mixture.csv"
FABRE = ROOT / "reference/data/fabre2022/ALLvariants_exclSynonymous_Xadj.txt"
FIGS = ROOT / "analysis/figures"

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRIDC, BASE = "#e1e0d9", "#c3c2b7"
RAMP = ["#a9c9f2", "#6ba3e8", "#2a78d6", "#154e91"]
S2 = "#eb6834"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRIDC, "grid.linewidth": 0.8, "axes.grid": True,
    "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 9,
    "legend.labelcolor": INK2, "figure.dpi": 130,
})

def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def section(t):
    print(f"\n{'='*76}\n{t}\n{'='*76}")

QS = [0.10, 0.25, 0.50, 0.75, 0.90]
CS = [f"cs_vaf_q{int(q*100):02d}" for q in QS]
TR = [f"tr_rate_q{int(q*100):02d}" for q in QS]
FLOOR = 0.002
TOL = 0.25          # the pre-registered band


# ---------------------------------------------------------------- observed
d = pd.read_csv(FABRE, sep="\t")
d = d[d.VAF >= FLOOR].copy()
d["clone"] = (d.SardID.astype(str) + "|" + d.Gene + "|"
              + d.START.astype(str) + "|" + d.ALT.astype(str))
first = d.sort_values("Phase").groupby("SardID").Phase.min().rename("p0")
b = d.merge(first, on="SardID")
OBS_CS = np.quantile(b[b.Phase == b.p0].VAF.values, QS)

rates = []
for c, g in d.groupby("clone"):
    g = g.sort_values("Age")
    if len(g) < 3 or np.ptp(g.Age) < 8 or g.VAF.iloc[0] < FLOOR:
        continue
    rates.append(float(np.polyfit(g.Age, np.log(g.VAF), 1)[0]))
rates = np.array(rates)
OBS_TR = np.quantile(rates, QS)
OBS_WIDTH = float(OBS_TR[4] - OBS_TR[0])
OBS_NEG = float((rates < 0).mean())

section("1. THE TARGETS")
print(f"  spectrum  VAF quantiles   {np.round(OBS_CS, 4)}")
print(f"  trajectory rate quantiles {np.round(OBS_TR, 4)}")
print(f"  width {OBS_WIDTH:.4f}   negative {OBS_NEG:.1%}   n={rates.size}")


# -------------------------------------------------------------------- grid
sw = pd.read_csv(GRID).dropna(subset=CS + TR)
sw["d_cs"] = np.sqrt(((np.log(sw[CS].values) - np.log(OBS_CS)) ** 2).mean(axis=1))
sw["d_tr"] = np.sqrt(((sw[TR].values - OBS_TR) ** 2).mean(axis=1))
g = sw.groupby(["mean", "shape", "N"], as_index=False).agg(
    d_cs=("d_cs", "mean"), d_tr=("d_tr", "mean"),
    width=("tr_width", "mean"), neg=("tr_negative", "mean"),
    vaf50=("cs_vaf_q50", "mean"), rate50=("tr_rate_q50", "mean"))
print(f"\n  grid: {len(g)} points, N from {g.N.min():,} to {g.N.max():,}")


section("2. THE IDENTIFIABILITY CHECKS THE PLAN REQUIRES FIRST")

folds = []
for (sh, N), sub in g.groupby(["shape", "N"]):
    sub = sub.sort_values("mean")
    v = sub.vaf50.values
    folds.append(((np.diff(v) >= 0).all(), sh, N))
n_mono = sum(1 for f in folds if f[0])
print(f"  spectrum arm monotonic in the mean: {n_mono} of {len(folds)} (shape, N) slices")
if n_mono < len(folds):
    bad = [f"shape {f[1]}, N {f[2]:,}" for f in folds if not f[0]][:4]
    print(f"    folded in: {'; '.join(bad)}")

joint = np.sqrt((g.d_cs / g.d_cs.min()) ** 2 + (g.d_tr / g.d_tr.min()) ** 2)
g = g.assign(joint=joint)
best = g.loc[g.joint.idxmin()]
near = g[g.joint <= g.joint.min() * 1.05]
print(f"\n  joint minimum at mean {best['mean']:.2f}, shape {best['shape']}, N {int(best.N):,}")
print(f"  points within 5% of it: {len(near)}  -> "
      f"{'UNIQUE enough' if len(near) <= 3 else 'NOT unique, a plateau'}")


section("3. THE PRE-REGISTERED DECISION")

best_cs, best_tr = g.d_cs.min(), g.d_tr.min()
ok = g[(g.d_cs <= best_cs * (1 + TOL)) & (g.d_tr <= best_tr * (1 + TOL))]
print(f"  best distance achievable   spectrum {best_cs:.4f}   trajectory {best_tr:.4f}")
print(f"  band: within {TOL:.0%} of each  ->  {len(ok)} of {len(g)} points qualify\n")

if len(ok):
    print("  A JOINTLY-ACCEPTABLE REGION EXISTS.\n")
    print(f"  {'mean':>6}{'shape':>7}{'N':>10}{'d_cs':>8}{'d_tr':>8}"
          f"{'VAF q50':>10}{'width':>9}{'neg':>7}")
    for _, r in ok.sort_values("joint").head(10).iterrows():
        print(f"  {r['mean']:>6.2f}{r['shape']:>7.1f}{int(r.N):>10,}{r.d_cs:>8.4f}"
              f"{r.d_tr:>8.4f}{r.vaf50:>10.4f}{r.width:>9.4f}{r.neg:>7.1%}")
    print(f"""
  observed for comparison: VAF q50 {OBS_CS[2]:.4f}, width {OBS_WIDTH:.4f}, negative {OBS_NEG:.1%}
  N values in the accepted region: {sorted(ok.N.unique())}""")
else:
    print("  NO POINT FITS BOTH.\n")
    print("  The closest on each arm, and what it costs on the other:")
    print(f"  {'':<18}{'mean':>6}{'shape':>7}{'N':>10}{'d_cs':>8}{'d_tr':>8}")
    for lab, row in (("best spectrum", g.loc[g.d_cs.idxmin()]),
                     ("best trajectory", g.loc[g.d_tr.idxmin()]),
                     ("best joint", best)):
        print(f"  {lab:<18}{row['mean']:>6.2f}{row['shape']:>7.1f}{int(row.N):>10,}"
              f"{row.d_cs:>8.4f}{row.d_tr:>8.4f}")
    print(f"""
  The model cannot reproduce both the size distribution and the growth-rate
  distribution of one cohort, even with the stem-cell population free across an
  eightfold range. That is a structural inadequacy, and naming it is the
  pre-registered outcome — not tuning around it.""")

print("\n  Either way, per the pre-registration: STOP.")


section("4. WHAT N DID")
byN = g.groupby("N").agg(best_cs=("d_cs", "min"), best_tr=("d_tr", "min"),
                         best_joint=("joint", "min")).reset_index()
print(f"  {'N':>10}{'best d_cs':>12}{'best d_tr':>12}{'best joint':>12}")
for _, r in byN.iterrows():
    print(f"  {int(r.N):>10,}{r.best_cs:>12.4f}{r.best_tr:>12.4f}{r.best_joint:>12.3f}")
print("""
  N was taken from the literature in stage C and never fitted. This is the first
  time the project has let it move.""")


# ------------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.4))

ax = axes[0]
for i, N in enumerate(sorted(g.N.unique())):
    sub = g[g.N == N]
    ax.scatter(sub.d_cs, sub.d_tr, s=34, color=RAMP[i % len(RAMP)],
               label=f"N = {int(N):,}", zorder=3, edgecolor=SURF, linewidth=1)
ax.axvline(best_cs * (1 + TOL), color=S2, linewidth=2, zorder=2)
ax.axhline(best_tr * (1 + TOL), color=S2, linewidth=2, zorder=2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("distance, VAF spectrum")
ax.set_ylabel("distance, growth rates")
ax.set_title("Can any point be good at both?")
ax.legend(loc="upper right")
clean(ax)

ax = axes[1]
for i, N in enumerate(sorted(g.N.unique())):
    sub = g[(g.N == N) & (g["shape"] == 1.0)].sort_values("mean")
    ax.plot(sub["mean"], sub.vaf50, "-o", color=RAMP[i % len(RAMP)],
            linewidth=2, markersize=6, markeredgecolor=SURF, markeredgewidth=1.5,
            label=f"N = {int(N):,}")
ax.axhline(OBS_CS[2], color=S2, linewidth=2, zorder=1)
ax.text(g["mean"].max(), OBS_CS[2], " observed", color=S2, fontsize=9, va="bottom", ha="right")
ax.set_yscale("log")
ax.set_xlabel("mean of the fitness distribution")
ax.set_ylabel("median detected VAF")
ax.set_title("What N moves (shape = 1.0)")
ax.legend(loc="upper left")
clean(ax)

fig.tight_layout(); fig.savefig(FIGS / "17_mixture_fit.png"); plt.close(fig)
print("\n  17_mixture_fit.png")
