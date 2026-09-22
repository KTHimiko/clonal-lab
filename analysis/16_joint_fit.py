#!/usr/bin/env python3
"""
Stage O — the joint fit, decided by a rule fixed before the sweep ran.

analysis/PREREGISTRATION_STAGE_O.md, commit 0d4e080:

    a point "fits both" if it is within 25% of the best distance achieved on
    EACH arm. If one exists, report it and then check the three statistics the
    fit does NOT optimise for. If none does, four mechanisms and a joint search
    have failed and the stage N bracket stands. Either way, stop.

Usage:  .venv/bin/python analysis/16_joint_fit.py [path/to/sweep_joint.csv]
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
GRID = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sweep_joint.csv"
FIGS = ROOT / "analysis/figures"

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRIDC, BASE = "#e1e0d9", "#c3c2b7"
RAMP = ["#a9c9f2", "#6ba3e8", "#2a78d6", "#154e91"]
S2, S3 = "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "grid.color": GRIDC,
    "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

def section(t):
    print(f"\n{'='*76}\n{t}\n{'='*76}")

QS = [0.10, 0.25, 0.50, 0.75, 0.90]
CS = [f"cs_vaf_q{int(q*100):02d}" for q in QS]
TR = [f"tr_rate_q{int(q*100):02d}" for q in QS]
TOL = 0.25

OBS_CS = np.array([0.0029, 0.0045, 0.0072, 0.0147, 0.0401])
OBS_TR = np.array([-0.0451, -0.0002, 0.0531, 0.1006, 0.1548])
OBS = dict(declining=0.251, icc=0.129, pff=0.51, pfr=0.36)

sw = pd.read_csv(GRID)
# rows with no hazard ignore the onset, so their duplicates collapse
sw.loc[sw.rate == 0, "from_age"] = 0.0
sw = sw.dropna(subset=CS + TR)
g = sw.groupby(["mean", "shape", "N", "rate", "from_age"], as_index=False).agg(
    **{c: (c, "mean") for c in CS + TR},
    declining=("tr_negative", "mean"), icc=("tr_icc", "mean"),
    pff=("tr_pff", "mean"), pfr=("tr_pfr", "mean"), n=("tr_n", "mean"))
g["d_cs"] = np.sqrt(((np.log(g[CS].values) - np.log(OBS_CS)) ** 2).mean(axis=1))
g["d_tr"] = np.sqrt(((g[TR].values - OBS_TR) ** 2).mean(axis=1))

section("1. THE GRID")
print(f"  {len(g)} points after averaging replicates")
print(f"  hazard rates: {sorted(g.rate.unique())}")
print(f"  N: {[f'{int(v):,}' for v in sorted(g.N.unique())]}")


section("2. IDENTIFIABILITY, CHECKED BEFORE ANY FIT IS REPORTED")
folds = []
for key, sub in g.groupby(["shape", "N", "rate", "from_age"]):
    sub = sub.sort_values("mean")
    if len(sub) >= 3:
        folds.append((np.diff(sub["cs_vaf_q50"].values) >= 0).all())
print(f"  spectrum monotonic in the mean: {sum(folds)} of {len(folds)} slices")
if sum(folds) < len(folds):
    print("    NOT monotonic everywhere — a fit in a folded slice is not an estimate")


section("3. THE PRE-REGISTERED DECISION")
best_cs, best_tr = g.d_cs.min(), g.d_tr.min()
ok = g[(g.d_cs <= best_cs * (1 + TOL)) & (g.d_tr <= best_tr * (1 + TOL))]
print(f"  best achievable   spectrum {best_cs:.4f}   rates {best_tr:.4f}")
print(f"  within {TOL:.0%} of each: {len(ok)} of {len(g)} points\n")

if len(ok) == 0:
    print("  NO POINT FITS BOTH.")
    print("  Four mechanisms and a joint search have failed. The stage N bracket")
    print("  stands as the result, and the missing mechanism is none of them.")
    for lab, row in (("best spectrum", g.loc[g.d_cs.idxmin()]),
                     ("best rates", g.loc[g.d_tr.idxmin()])):
        print(f"    {lab:<16} mean {row['mean']:.2f} shape {row['shape']:.0f} "
              f"N {int(row.N):,} rate {row.rate:.2f} from {row.from_age:.0f}  "
              f"d_cs {row.d_cs:.3f}  d_tr {row.d_tr:.3f}")
else:
    ok = ok.assign(joint=np.hypot(ok.d_cs / best_cs, ok.d_tr / best_tr)).sort_values("joint")
    print("  A JOINTLY-ACCEPTABLE REGION EXISTS.\n")
    print(f"  {'mean':>6}{'shape':>7}{'N':>10}{'rate':>7}{'from':>6}"
          f"{'d_cs':>8}{'d_tr':>8}")
    for _, r in ok.head(8).iterrows():
        print(f"  {r['mean']:>6.2f}{r['shape']:>7.0f}{int(r.N):>10,}{r.rate:>7.2f}"
              f"{r.from_age:>6.0f}{r.d_cs:>8.4f}{r.d_tr:>8.4f}")
    frac_haz = float((ok.rate > 0).mean())
    print(f"\n  points in the region that carry a hazard: {frac_haz:.0%}")
    if frac_haz < 0.5:
        print("  MOST OF THE ACCEPTED REGION HAS NO HAZARD — the mechanism is not")
        print("  what put those points there, and section 4 is what decides it.")

    section("4. THE THREE STATISTICS THE FIT DOES NOT OPTIMISE FOR")
    print("  These were not used to choose any point. They are the test.\n")
    print(f"  {'':<22}{'declining':>11}{'ICC':>9}{'P(f|f)':>9}{'P(f|r)':>9}")
    print(f"  {'observed':<22}{OBS['declining']:>10.1%}{OBS['icc']:>9.3f}"
          f"{OBS['pff']:>9.2f}{OBS['pfr']:>9.2f}")
    base = g[(g.rate == 0)].sort_values("d_cs").iloc[0]
    print(f"  {'best with no hazard':<22}{base.declining:>10.1%}{base.icc:>9.3f}"
          f"{base.pff:>9.2f}{base.pfr:>9.2f}")
    for _, r in ok.head(3).iterrows():
        lab = f"rate {r.rate:.2f} from {r.from_age:.0f}"
        print(f"  {lab:<22}{r.declining:>10.1%}{r.icc:>9.3f}{r.pff:>9.2f}{r.pfr:>9.2f}")

print("\n  Either way, per the pre-registration: STOP.")


section("5. WHAT THE HAZARD BUYS, ACROSS THE WHOLE GRID")
by = g.groupby("rate").agg(best_cs=("d_cs", "min"), best_tr=("d_tr", "min"),
                           decl=("declining", "max"), icc=("icc", "median")).reset_index()
print(f"  {'rate':>7}{'best d_cs':>12}{'best d_tr':>12}{'max declining':>15}{'median ICC':>12}")
for _, r in by.iterrows():
    print(f"  {r.rate:>7.2f}{r.best_cs:>12.4f}{r.best_tr:>12.4f}"
          f"{r.decl:>14.1%}{r.icc:>12.3f}")

fig, ax = plt.subplots(figsize=(8.6, 4.6))
for i, rate in enumerate(sorted(g.rate.unique())):
    sub = g[g.rate == rate]
    ax.scatter(sub.d_cs, sub.d_tr, s=26, color=RAMP[i % len(RAMP)],
               label=f"hazard {rate}", zorder=3, edgecolor=SURF, linewidth=0.8)
ax.axvline(best_cs * (1 + TOL), color=S2, linewidth=2, zorder=2)
ax.axhline(best_tr * (1 + TOL), color=S2, linewidth=2, zorder=2)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("distance, VAF spectrum")
ax.set_ylabel("distance, growth rates")
ax.set_title("Population and hazard, swept together")
ax.legend(loc="upper right")
clean(ax)
fig.tight_layout(); fig.savefig(FIGS / "18_joint_fit.png"); plt.close(fig)
print("\n  18_joint_fit.png")
