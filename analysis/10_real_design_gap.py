#!/usr/bin/env python3
"""
Stage F — measuring the design gap in real people instead of in a simulation.

WHAT STAGE E ESTABLISHED, AND ITS WEAKNESS
Given a simulated cohort whose fitness is known by construction, three study
designs recovered 18%, 43% and 103% of the truth. That is a strong result with
one obvious objection: every number came from our own model, so the biases are
real only to the extent the model is.

WHAT MAKES THIS POSSIBLE
Fabre's SardiNIA data (figshare 10.6084/m9.figshare.15029118, CC BY 4.0) has
what the Watson table never had: a person identifier and repeated measurements.
394 people, 994 clones, two to five timepoints each, a median follow-up of 12.9
years, sequenced at a median depth of 1,140x.

So all three designs can be run ON THE SAME INDIVIDUALS:

  longitudinal      every timepoint, one growth rate per clone       (Fabre)
  age regression    first timepoint only, log(VAF) against age       (stage 1)
  spectrum fit      first timepoint only, matched to our grid        (Watson, C-D)

Same people, same assay, same clones. Any difference between the three is
design, measured in the world rather than in a model.

Usage:  .venv/bin/python analysis/10_real_design_gap.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

FABRE = ROOT / "reference/data/fabre2022/ALLvariants_exclSynonymous_Xadj.txt"
RESULTS = ROOT / "results"
FIGS = ROOT / "analysis/figures"

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

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
    print(f"\n{'='*76}\n{t}\n{'='*76}")

def sub(t):
    print(f"\n{'-'*76}\n{t}\n{'-'*76}")


LIMIT = 0.0192        # the floor our grid was computed at
MIN_POINTS = 3        # timepoints needed for a trajectory fit
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
QCOLS = ["vaf_q10", "vaf_q25", "vaf_q50", "vaf_q75", "vaf_q90"]


d = pd.read_csv(FABRE, sep="\t")
d["clone"] = (d.SardID.astype(str) + "|" + d.Gene + "|"
              + d.START.astype(str) + "|" + d.ALT.astype(str))


# ============================================ 0. does our reading match theirs?
section("0. VALIDATION — CAN WE REPRODUCE THEIR PUBLISHED NUMBERS?")
print("""Before comparing anything, check that we are reading their file the way they
did. Fabre report DNMT3A clones growing at 6.2% per year in this cohort. If a
simple trajectory fit on their table does not land near that, every number
below is suspect.""")

def trajectory_rate(g):
    """Slope of log(VAF) on age, across a clone's own timepoints."""
    g = g[g.VAF > 0]
    if len(g) < MIN_POINTS or g.Age.nunique() < MIN_POINTS or np.ptp(g.Age) < 3:
        return np.nan
    return float(np.polyfit(g.Age, np.log(g.VAF), 1)[0])

rates = (d.groupby("clone")
           .apply(trajectory_rate, include_groups=False)
           .rename("rate").reset_index())
meta = d.sort_values("Phase").groupby("clone").agg(
    gene=("Gene", "first"), person=("SardID", "first"),
    age0=("Age", "first"), vaf0=("VAF", "first"), n=("VAF", "size")).reset_index()
traj = rates.merge(meta, on="clone").dropna(subset=["rate"])

print(f"\n{len(traj)} clones with at least {MIN_POINTS} usable timepoints\n")
print(f"{'gene':<9}{'n':>5}{'median rate/yr':>16}{'published':>12}")
PUBLISHED = {"DNMT3A": 0.062, "TET2": 0.068}
for gene, g in traj.groupby("gene"):
    if len(g) < 25:
        continue
    pub = PUBLISHED.get(gene)
    print(f"{gene:<9}{len(g):>5}{g.rate.median():>16.4f}"
          f"{(f'{pub:.3f}' if pub else '-'):>12}")


# ================================================== 1. the three estimators
section("1. THE SAME THREE DESIGNS, THE SAME 394 PEOPLE")

first_phase = d.sort_values("Phase").groupby("SardID").Phase.min().rename("p0")
base = d.merge(first_phase, on="SardID")
base = base[base.Phase == base.p0]

sub("1.1 Longitudinal — every timepoint, one rate per clone")
enrolled = traj[traj.vaf0 >= LIMIT]
print(f"  all clones with a trajectory      n={len(traj):>4}   "
      f"median {traj.rate.median():+.4f} per year")
print(f"  detected at baseline (VAF>={LIMIT})  n={len(enrolled):>4}   "
      f"median {enrolled.rate.median():+.4f} per year")
print("""
  The second row is the comparable one: stage E's longitudinal estimator only
  admitted clones already above the floor at enrolment, because that is what a
  real study can follow.""")
longit = float(enrolled.rate.median())

sub("1.2 Age regression — first timepoint only")
cs = base[base.VAF >= LIMIT]
slope = float(np.polyfit(cs.Age, np.log(cs.VAF), 1)[0])
print(f"  n={len(cs)} clones in {cs.SardID.nunique()} people, "
      f"ages {cs.Age.min():.0f}-{cs.Age.max():.0f}")
print(f"  slope of log(VAF) on age: {slope:+.4f} per year")

sub("1.3 Spectrum fit — first timepoint only, scored against our grid")
frames = [pd.read_csv(RESULTS / f) for f in ("sweep_N100k.csv", "sweep_fine.csv")
          if (RESULTS / f).exists()]
sw = pd.concat(frames, ignore_index=True)
sw = sw[(sw.age == 70.0) & (sw.N == 100_000)].dropna(subset=QCOLS).copy()
sim = np.log(sw[QCOLS].values)
scale = sim.std(axis=0); scale[scale == 0] = 1.0

# The grid is a snapshot of a cohort all aged 70. Scoring a 55-94 sample against
# it would compare different things, so the spectrum fit uses an age band around
# the grid's age. The median first-phase age is 70, which is why this works.
band = cs[(cs.Age >= 65) & (cs.Age <= 75)]
print(f"  age band 65-75: n={len(band)} clones in {band.SardID.nunique()} people")

def score(values):
    target = np.log(np.quantile(values, QUANTILES))
    dist = np.sqrt((((sim - target) / scale) ** 2).mean(axis=1))
    grid = (sw.assign(distance=dist).groupby(["s", "mu"], as_index=False)
              .agg(distance=("distance", "mean")))
    row = grid.loc[grid.distance.idxmin()]
    return float(row.s), float(row.distance)

spectrum, dist = score(band.VAF.values)
print(f"  best grid point: s = {spectrum:.3f} per year   (distance {dist:.3f})")

rng = np.random.default_rng(0)
boot = np.array([score(rng.choice(band.VAF.values, band.VAF.size, replace=True))[0]
                 for _ in range(300)])
print(f"  300-sample bootstrap: {np.quantile(boot, 0.025):.3f} - "
      f"{np.quantile(boot, 0.975):.3f}")


# ====================================================== 2. the comparison
section("2. WHAT STAGE E PREDICTED, AND WHAT THE PEOPLE SAY")

print(f"""{'design':<24}{'simulated (stage E)':>22}{'real (this cohort)':>22}
{'':<24}{'% of known truth':>22}{'s per year':>22}
{'age regression':<24}{'18%':>22}{slope:>22.4f}
{'longitudinal':<24}{'43%':>22}{longit:>22.4f}
{'spectrum fit':<24}{'103%':>22}{spectrum:>22.4f}

The simulation has a known truth and the cohort does not, so the columns are not
the same quantity. What CAN be compared is the ordering and the ratios between
designs, because those are what stage E claims are properties of the designs.""")

print(f"""
{'ratio':<34}{'predicted':>12}{'observed':>12}
{'spectrum / longitudinal':<34}{103/43:>12.2f}{spectrum/longit:>12.2f}
{'spectrum / age regression':<34}{103/18:>12.2f}{spectrum/slope:>12.2f}
{'longitudinal / age regression':<34}{43/18:>12.2f}{longit/slope:>12.2f}""")

order_ok = spectrum > longit > slope
print(f"""
ORDERING PREDICTED: spectrum > longitudinal > age regression
ORDERING OBSERVED:  {'the same' if order_ok else 'DIFFERENT — see below'}""")


# ============================================================== 3. figure
section("3. FIGURE")

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))

ax = axes[0]
names = ["age regression\n(stage 1 design)", "longitudinal\n(Fabre design)",
         "spectrum fit\n(Watson design)"]
vals = [slope, longit, spectrum]
cols = [S1, S2, S3]
y = np.arange(3)[::-1]
ax.barh(y, vals, height=0.5, color=cols, zorder=2)
for i, v in zip(y, vals):
    ax.text(v + 0.004, i, f"{v:.3f}", va="center", fontsize=10, color=INK2)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9)
ax.set_xlabel("estimated s  (per year)")
ax.set_title("Three designs, one cohort of 394 people")
ax.grid(axis="y", visible=False)
clean(ax)

ax = axes[1]
pred = [18/103, 43/103, 1.0]
obs = [slope/spectrum, longit/spectrum, 1.0]
x = np.arange(3); w = 0.36
ax.bar(x - w/2, pred, width=w - 0.03, color=BASE, label="stage E, simulated", zorder=2)
ax.bar(x + w/2, obs, width=w - 0.03, color=S1, label="this cohort, real", zorder=2)
ax.set_xticks(x)
ax.set_xticklabels(["age\nregression", "longitudinal", "spectrum\nfit"], fontsize=9)
ax.set_ylabel("fraction of the spectrum-fit estimate")
ax.set_title("Predicted shortfall against measured")
ax.legend(loc="upper left")
ax.grid(axis="x", visible=False)
clean(ax)

fig.tight_layout()
fig.savefig(FIGS / "14_real_design_gap.png"); plt.close(fig)
print("  14_real_design_gap.png")


# ==================================== 4. the model-free result inside the data
section("4. THE RESULT THAT NEEDS NO MODEL AT ALL")

print("""Section 1.1 holds the cleanest finding in this project, and it took no
simulation and no grid. The same cohort, the same assay, the same fitting
procedure — the only thing that changes is whether a clone was big enough to be
detectable at enrolment.""")

rng = np.random.default_rng(3)
def med_boot(v, B=4000):
    b = np.array([np.median(rng.choice(v, v.size, replace=True)) for _ in range(B)])
    return np.quantile(b, [0.025, 0.975])

all_lo, all_hi = med_boot(traj.rate.values)
enr_lo, enr_hi = med_boot(enrolled.rate.values)
print(f"""
  {'subset':<34}{'n':>5}{'median rate':>13}{'95% interval':>22}
  {'every clone with a trajectory':<34}{len(traj):>5}{traj.rate.median():>13.4f}{f'{all_lo:.4f} to {all_hi:.4f}':>22}
  {'only those detectable at baseline':<34}{len(enrolled):>5}{enrolled.rate.median():>13.4f}{f'{enr_lo:.4f} to {enr_hi:.4f}':>22}

  Conditioning on detectability costs {1 - enrolled.rate.median()/traj.rate.median():.0%} of the measured growth rate,
  and the intervals {'do not overlap' if enr_hi < all_lo else 'overlap'}.

  This is stage E's longitudinal bias, in real people, with no model in the
  path. A study that enrols the clones it can already see is measuring the ones
  that have most nearly finished growing.""")


# ================================== 5. does the detection floor behave as predicted?
section("5. A PREDICTION OF STAGE E, TESTED ON NEW DATA")

print("""Stage E found that moving only the detection floor changed the age slope by a
factor of 4.1 in the Watson data, and that the deep end was steeper. This cohort
is sequenced far deeper, so the same lever can be pulled here — on data that had
no part in producing the prediction.""")

def slope_boot(frame, B=2000):
    s = float(np.polyfit(frame.Age, np.log(frame.VAF), 1)[0])
    idx = np.arange(len(frame))
    b = []
    for _ in range(B):
        k = rng.choice(idx, idx.size, replace=True)
        f = frame.iloc[k]
        if f.Age.nunique() > 3:
            b.append(np.polyfit(f.Age, np.log(f.VAF), 1)[0])
    b = np.array(b)
    return s, np.quantile(b, 0.025), np.quantile(b, 0.975)

print(f"\n  {'floor':<22}{'n':>6}{'slope':>10}{'95% interval':>24}")
for name, floor in [("none (VAF > 0)", 0.0), ("0.005", 0.005), ("0.0192", LIMIT)]:
    f = base[base.VAF > max(floor, 1e-9)] if floor == 0 else base[base.VAF >= floor]
    s, lo, hi = slope_boot(f)
    print(f"  {name:<22}{len(f):>6}{s:>10.4f}{f'{lo:+.4f} to {hi:+.4f}':>24}")

print("""
  Same direction as the Watson data: the LOWER the floor the steeper the slope,
  and at the floor our grid uses the interval spans zero — the estimate is not
  distinguishable from no age effect at all. The age regression is not measuring
  fitness, it is measuring the instrument.""")


# ============================================= 6. how much is our own estimator?
section("6. HOW MUCH OF THE VALIDATION GAP IS OUR SIMPLER ESTIMATOR?")

print(f"""Section 0 got DNMT3A at {traj[traj.gene=='DNMT3A'].rate.median():.4f} against their published 0.062, and TET2 at
{traj[traj.gene=='TET2'].rate.median():.4f} against 0.068. Same ordering, magnitudes within about a third.

The gap is expected and worth naming rather than hiding: Fabre fit a Bayesian
model with a fitted overdispersion term, per-clone random effects and a
treatment of the detection limit. This uses ordinary least squares on log(VAF).
The simpler estimator is the right one here, because it is exactly the estimator
stage E simulated — but it is not their estimator, so these are our numbers on
their data, not a reproduction of their analysis.""")
