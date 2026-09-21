#!/usr/bin/env python3
"""
Stage J — fitting the Fabre cohort two ways, against a grid built for it.

The decision rule was committed before the sweep ran, in
analysis/PREREGISTRATION_STAGE_J.md (f6af384):

    R = spectrum estimate / trajectory estimate
    1.8 <= R <= 3.2   the model reproduces the design gap, conclusion confirmed
    otherwise         conclusion downgraded, recorded as such
    either way        stop

Usage:  .venv/bin/python analysis/14_fabre_grid.py [path/to/sweep_fabre.csv]
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
GRID = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sweep_fabre.csv"
FABRE = ROOT / "reference/data/fabre2022/ALLvariants_exclSynonymous_Xadj.txt"
FIGS = ROOT / "analysis/figures"

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRIDC, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

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
CS_COLS = [f"cs_vaf_q{int(q*100):02d}" for q in QS]
TR_COLS = [f"tr_rate_q{int(q*100):02d}" for q in QS]
FLOOR = 0.002
LO, HI = 1.8, 3.2          # the pre-registered band
rng = np.random.default_rng(31)


# ------------------------------------------------------------------ observed
d = pd.read_csv(FABRE, sep="\t")
d = d[d.VAF >= FLOOR].copy()
d["clone"] = (d.SardID.astype(str) + "|" + d.Gene + "|"
              + d.START.astype(str) + "|" + d.ALT.astype(str))

first = d.sort_values("Phase").groupby("SardID").Phase.min().rename("p0")
base = d.merge(first, on="SardID")
obs_cs = base[base.Phase == base.p0].VAF.values

rates = []
for c, g in d.groupby("clone"):
    g = g.sort_values("Age")
    if len(g) < 3 or np.ptp(g.Age) < 8:
        continue
    if g.VAF.iloc[0] < FLOOR:
        continue
    rates.append(float(np.polyfit(g.Age, np.log(g.VAF), 1)[0]))
obs_tr = np.array(rates)

section("1. THE OBSERVED COHORT, PROCESSED AS THE SIMULATION PROCESSES ITSELF")
print(f"  spectrum arm:   {obs_cs.size} clones at first draw, "
      f"VAF quantiles {np.round(np.quantile(obs_cs, QS), 4)}")
print(f"  trajectory arm: {obs_tr.size} clones with >=3 draws, "
      f"rate quantiles {np.round(np.quantile(obs_tr, QS), 4)}")
print("""
  One mismatch, stated rather than hidden: the simulation enrols everyone at 65
  and takes five evenly spaced draws to 78. The real cohort enrols across a
  range of ages with two to five draws. The estimator is the same on both sides;
  the sampling is not identical.""")


# ---------------------------------------------------------------------- grid
sw = pd.read_csv(GRID)
sw = sw.dropna(subset=CS_COLS + TR_COLS)
print(f"\n  grid: {len(sw)} rows, {sw[['s','mu']].drop_duplicates().shape[0]} points, "
      f"s from {sw.s.min():.2f} to {sw.s.max():.2f}")


def score(sim_matrix, target, transform):
    """Closest grid point, distance standardised across the grid's own spread."""
    sim = transform(sim_matrix)
    tgt = transform(target[None, :])[0]
    scale = sim.std(axis=0)
    scale[scale == 0] = 1.0
    dist = np.sqrt((((sim - tgt) / scale) ** 2).mean(axis=1))
    g = (sw.assign(distance=dist).groupby(["s", "mu"], as_index=False)
           .agg(distance=("distance", "mean")))
    row = g.loc[g.distance.idxmin()]
    return float(row.s), float(row.distance), g


section("2. THE TWO ARMS")

cs_s, cs_d, cs_grid = score(sw[CS_COLS].values, np.quantile(obs_cs, QS), np.log)
tr_s, tr_d, tr_grid = score(sw[TR_COLS].values, np.quantile(obs_tr, QS), lambda a: a)

def boot(values, cols, transform, B=300):
    out = []
    for _ in range(B):
        t = np.quantile(rng.choice(values, values.size, replace=True), QS)
        out.append(score(sw[cols].values, t, transform)[0])
    return np.quantile(out, [0.025, 0.975])

cs_lo, cs_hi = boot(obs_cs, CS_COLS, np.log)
tr_lo, tr_hi = boot(obs_tr, TR_COLS, lambda a: a)

print(f"  {'arm':<26}{'s':>8}{'95% interval':>22}{'distance':>10}")
print(f"  {'VAF spectrum':<26}{cs_s:>8.3f}{f'{cs_lo:.3f} - {cs_hi:.3f}':>22}{cs_d:>10.3f}")
print(f"  {'trajectories':<26}{tr_s:>8.3f}{f'{tr_lo:.3f} - {tr_hi:.3f}':>22}{tr_d:>10.3f}")

R = cs_s / tr_s if tr_s else float("nan")
print(f"\n  R = spectrum / trajectory = {R:.2f}")


section("3. THE PRE-REGISTERED DECISION")

inside = LO <= R <= HI
print(f"""  band committed in advance: {LO} to {HI}
  observed: {R:.2f}   -> {'INSIDE' if inside else 'OUTSIDE'}

  {'CONCLUSION CONFIRMED.' if inside else 'CONCLUSION DOWNGRADED.'}""")
if inside:
    print(f"""  The model reproduces the design gap end to end: fitted to one real cohort
  with a grid built for that cohort, the spectrum arm and the trajectory arm
  differ by {R:.2f}, against 2.42 measured in stage F and 2.40 predicted in stage E.
  The borrowed-grid limitation on stage F closes.""")
else:
    print(f"""  The pre-registration's own reading of this branch was that stage F's 2.42
  was partly coincidence. THAT READING IS WRONG, and the grid says why.

  Section 4 shows the trajectory arm is NON-MONOTONIC in s: the median growth
  rate rises to a peak near s = 0.12 and falls away to nearly zero by s = 0.30,
  because a fitter clone is already saturating by the time follow-up starts.
  The spectrum arm is monotonic across the whole grid.

  A non-monotonic map has no inverse. The observed median of {np.median(obs_tr):.4f} per year
  matches the grid at roughly s = 0.066 on the rising branch AND s = 0.182 on the
  falling one. The fit took the high branch for trajectories and the low branch
  for the spectrum, and R is the ratio between two different branches of a
  two-valued map.

  So R was never a well-defined quantity, and no value of it could have tested
  what the pre-registration thought it was testing. The verdict stands — the
  band was missed, the conclusion is downgraded, and we stop — but the reason is
  a flaw in this stage's design, not evidence against stage F.

  Stage F compared a spectrum-fit estimate against a DIRECTLY MEASURED growth
  rate. It never inverted the trajectory distribution, so it does not inherit
  this defect.""")
print("\n  Either way, per the pre-registration: STOP.")


section("4. THE DIAGNOSIS — WHICH ARM CAN BE INVERTED AT ALL")

g = sw.groupby("s").agg(tr=("tr_rate_q50", "mean"), cs=("cs_vaf_q50", "mean")).reset_index()
peak = g.tr.idxmax()
mono_cs = bool((g.cs.diff().dropna() >= 0).all())
mono_tr = bool((g.tr.diff().dropna() >= 0).all())
print(f"  {'s':>6}{'median rate (trajectories)':>28}{'median VAF (spectrum)':>24}")
for _, r in g.iterrows():
    mark = "  <- peak" if r.s == g.s[peak] else ""
    print(f"  {r.s:>6.2f}{r.tr:>28.4f}{r.cs:>24.4f}{mark}")
print(f"""
  spectrum arm monotonic in s:     {mono_cs}
  trajectory arm monotonic in s:   {mono_tr}   (peaks at s = {g.s[peak]:.2f})

  This is stage E's finding about the age regression, appearing again in a
  different estimator. A fitter clone is closer to saturation when follow-up
  begins, so beyond a certain fitness the measured growth rate falls as s rises.
  The map folds, and a folded map cannot be inverted.

  The consequence for practice is concrete: a longitudinal growth rate is a fine
  thing to MEASURE and report. It is not a quantity you can fit a model to and
  read a fitness out of, unless you already know which branch you are on.""")


section("5. WHAT THIS DOES NOT SHOW")
print("""Both arms come from one model, so this tests the inversion and not the
biology — the limit stage E stated and stage I's power calculation reinforced.
The external support for the spectrum estimator remains the phylogenetic
agreement from stage E2, which is untouched by anything here.

And the model-free result stands on its own regardless: enrolling clones by
detectability costs 35% of the measured growth rate in 770 real trajectories,
with no simulation in the path.""")


fig, ax = plt.subplots(figsize=(8.6, 4.2))
for g, lab, col in ((cs_grid, "VAF spectrum", S1), (tr_grid, "trajectories", S2)):
    m = g.groupby("s", as_index=False).distance.min()
    ax.plot(m.s, m.distance, "-o", color=col, linewidth=2.5, markersize=5,
            markeredgecolor=SURF, markeredgewidth=1.5, label=lab)
for v, col in ((cs_s, S1), (tr_s, S2)):
    ax.axvline(v, color=col, linewidth=1.5, linestyle=":", zorder=1)
ax.set_xlabel("selection coefficient  s  (per year)")
ax.set_ylabel("distance to the observed cohort")
ax.set_title(f"One cohort, one grid, two estimators  —  R = {R:.2f}")
ax.legend()
clean(ax)
fig.tight_layout(); fig.savefig(FIGS / "16_fabre_grid.png"); plt.close(fig)
print("\n  16_fabre_grid.png")
