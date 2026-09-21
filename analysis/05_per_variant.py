#!/usr/bin/env python3
"""
Stage D — fitness per variant class, not one number for everything.

Stage C fitted a single selection coefficient to every clone at once. Stage 1
had already shown that is wrong: TET2 grows about twice as fast as DNMT3A, and
within DNMT3A the R882 hotspot grows faster than the rest.

No new simulation is needed. The sweep already spans a grid of s, and each
grid point is a cohort in which every clone carries that s. Scoring that same
grid against a variant-restricted slice of the observed data asks: what
fitness, on its own, reproduces the size distribution of THIS class of clone?

That is the payoff of having kept simulation and scoring apart.

Usage:  .venv/bin/python analysis/05_per_variant.py [path/to/sweep.csv]
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
SWEEP = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/sweep_N100k.csv"
CSV = (ROOT / "reference/data/watson2020/Maximum_likelihood_estimations"
            / "Maximum likelihood estimations - data files"
            / "all_studies_trimmed_all_genes.csv")
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
    print(f"\n{'-'*76}\n{t}\n{'-'*76}")


LIMIT = 0.0192          # common floor, the least sensitive cohort we keep
COMPARE_AGE = 70.0
QCOLS = ["vaf_q10", "vaf_q25", "vaf_q50", "vaf_q75", "vaf_q90"]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90]
ACCEPT = 0.05

# Published values, for checking against. Watson is per variant from a
# cross-sectional VAF spectrum; Fabre is per gene from 697 clones followed in
# the same people for a median of 13 years. They disagree, which is itself the
# useful context: any estimate landing between them is defensible.
PUBLISHED = {
    "DNMT3A R882":  {"Watson 2020": 0.148, "Fabre 2022": 0.050},   # R882H
    "DNMT3A other": {"Watson 2020": 0.120, "Fabre 2022": 0.050},   # non-hotspot median of top-20
    "TET2":         {"Fabre 2022": 0.068},
    "JAK2":         {"Watson 2020": 0.146, "Fabre 2022": None},    # V617F
}


# ================================================================== observed
section("1. THE OBSERVED DATA, SPLIT BY VARIANT CLASS")

raw = pd.read_csv(CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
raw["age"] = pd.to_numeric(raw["age"], errors="coerce")
raw["VAF"] = pd.to_numeric(raw["VAF"], errors="coerce")
obs = raw.dropna(subset=["age", "VAF"]).copy()

obs["klass"] = np.where(
    obs.gene == "DNMT3A",
    np.where(obs.variant.astype(str).str.startswith("R882"), "DNMT3A R882", "DNMT3A other"),
    obs.gene)

# A common floor is not enough on its own. A cohort whose own detection limit
# sits ABOVE the floor contributes only large clones, and pooling it skews
# every class it touches upward. ZinkWGS sees nothing below VAF 0.24, so
# including it put a third of the TET2 class at the top of the range and drove
# that estimate to an implausible 0.35 per year.
#
# The rule: keep only cohorts able to see everything above the common floor.
COHORT_LIMITS = {"Young2019": 0.0008, "Young2016": 0.0011, "Acuna2017": 0.0016,
                 "McKerrel2015": 0.0079, "Desai2018": 0.0182,
                 "Coombs2017": 0.0192, "ZinkWGS": 0.24}
eligible = [c for c, l in COHORT_LIMITS.items() if l <= LIMIT]
dropped = [c for c, l in COHORT_LIMITS.items() if l > LIMIT]

obs = obs[obs.study.isin(eligible) & (obs.VAF >= LIMIT)]
classes = [k for k in ["DNMT3A R882", "DNMT3A other", "TET2", "JAK2"]
           if (obs.klass == k).sum() >= 30]

print(f"common detection limit: VAF >= {LIMIT}")
print(f"cohorts kept:    {', '.join(sorted(eligible))}")
print(f"cohorts dropped: {', '.join(sorted(dropped))}  (limit above the floor)")
print("\nTwo separate corrections. The floor discards clones the least")
print("sensitive kept cohort could not have seen. Dropping the cohorts whose")
print("own limit sits above that floor removes studies that contribute only")
print("large clones and would skew every class they touch upward.\n")
print(f"{'class':<14} {'n':>5}  {'median VAF':>11}  {'median age':>11}  {'cohorts':>8}")
for k in classes:
    g = obs[obs.klass == k]
    print(f"{k:<14} {len(g):>5}  {g.VAF.median():>11.4f}  {g.age.median():>11.0f}  {g.study.nunique():>8}")

print(f"\nAll class median ages fall between 66 and 72, so all are scored")
print(f"against the simulated snapshot at age {COMPARE_AGE:.0f}.")


# =============================================================== the sweep
section("2. SCORING THE EXISTING SWEEP AGAINST EACH CLASS")

sw = pd.read_csv(SWEEP)
sw = sw[sw.age == COMPARE_AGE].dropna(subset=QCOLS).copy()
print(f"sweep: {sw[['s','mu']].drop_duplicates().shape[0]} grid points, "
      f"no new simulation run")

sim = np.log(sw[QCOLS].values)
scale = sim.std(axis=0)
scale[scale == 0] = 1.0

def best_s_for(values):
    """The grid point whose VAF quantiles land closest to these observations."""
    target = np.log(np.quantile(values, QUANTILES))
    d = np.sqrt((((sim - target) / scale) ** 2).mean(axis=1))
    grid = (sw.assign(distance=d)
              .groupby(["s", "mu"], as_index=False)
              .agg(distance=("distance", "mean")))
    row = grid.loc[grid.distance.idxmin()]
    return row.s, row.distance, grid


BOOTSTRAP = 400
rng = np.random.default_rng(0)

results = {}
for k in classes:
    vals = obs[obs.klass == k].VAF.values
    best, dist, grid = best_s_for(vals)

    # How much of the class difference survives the sample size? Each class
    # holds between 36 and 153 variants. Resampling them with replacement and
    # re-fitting shows how much the estimate moves for that reason alone —
    # which is the difference between "we got two numbers" and "these two
    # numbers differ".
    boot = np.array([best_s_for(rng.choice(vals, size=vals.size, replace=True))[0]
                     for _ in range(BOOTSTRAP)])

    results[k] = {"best_s": best, "distance": dist,
                  "s_lo": np.quantile(boot, 0.025),
                  "s_hi": np.quantile(boot, 0.975),
                  "boot": boot,
                  "n": int(vals.size)}


# ================================================================== results
section("3. INFERRED FITNESS, CLASS BY CLASS")
print(f"{'class':<14} {'n':>5}  {'s inferred':>11}  {'95% bootstrap':>17}  {'distance':>9}")
for k in classes:
    r = results[k]
    ci = f"{r['s_lo']:.3f} - {r['s_hi']:.3f}"
    print(f"{k:<14} {r['n']:>5}  {r['best_s']:>11.3f}  {ci:>17}  {r['distance']:>9.4f}")

print("\npairwise: does the difference survive resampling?")
for i, a in enumerate(classes):
    for b in classes[i+1:]:
        diff = results[a]["boot"] - results[b]["boot"]
        frac = float((diff > 0).mean())
        gap = results[a]["best_s"] - results[b]["best_s"]
        call = "separated" if frac > 0.95 or frac < 0.05 else "NOT separated"
        print(f"  {a:<13} vs {b:<13}  gap {gap:+.3f}   "
              f"P(first > second) = {frac:.2f}   {call}")

print(f"\n{'class':<14}  {'ours':>8}  {'Watson 2020':>12}  {'Fabre 2022':>11}   verdict")
for k in classes:
    ours = results[k]["best_s"]
    pub = PUBLISHED.get(k, {})
    w = pub.get("Watson 2020")
    f = pub.get("Fabre 2022")
    ws = f"{w:.3f}" if w else "-"
    fs = f"{f:.3f}" if f else "-"
    vals = [v for v in (w, f) if v]
    if vals and min(vals) <= ours <= max(vals):
        verdict = "between the two methods"
    elif vals and ours > max(vals):
        verdict = "above both"
    elif vals and ours < min(vals):
        verdict = "below both"
    else:
        verdict = "no published pair"
    print(f"{k:<14}  {ours:>8.3f}  {ws:>12}  {fs:>11}   {verdict}")


# ============================================================== the ordering
section("4. THE ORDERING, WHICH MATTERS MORE THAN THE VALUES")
order = sorted(classes, key=lambda k: -results[k]["best_s"])
print("inferred, fastest first:   " + "  >  ".join(order))
print("\nStage 1 measured growth rates directly from the raw data, by")
print("regressing log(VAF) on age. It found TET2 growing about twice as fast")
print("as DNMT3A, and R882 clones nearly twice the size of other DNMT3A")
print("variants at the same median age.")
print("\nThose were measurements of the data. The ordering above comes from")
print("fitting a mechanistic model. Agreement between the two is a check that")
print("the model is picking up real signal rather than fitting noise.")


# ================================================================== figure
section("5. FIGURE")

fig, ax = plt.subplots(figsize=(8.4, 4.4))
y = np.arange(len(classes))[::-1]

for i, k in zip(y, classes):
    r = results[k]
    ax.plot([r["s_lo"], r["s_hi"]], [i, i], color=BASE, linewidth=6,
            solid_capstyle="round", zorder=1,
            label="95% bootstrap" if i == y[0] else None)
    ax.scatter([r["best_s"]], [i], s=110, color=S1, zorder=3,
               label="this work" if i == y[0] else None)
    pub = PUBLISHED.get(k, {})
    if pub.get("Watson 2020"):
        ax.scatter([pub["Watson 2020"]], [i], s=90, marker="D", color=S2, zorder=3,
                   label="Watson 2020 (cross-sectional)" if i == y[0] else None)
    if pub.get("Fabre 2022"):
        ax.scatter([pub["Fabre 2022"]], [i], s=90, marker="^", color=S3, zorder=3,
                   label="Fabre 2022 (longitudinal)" if i == y[0] else None)

ax.set_yticks(y)
ax.set_yticklabels([f"{k}\n(n={results[k]['n']})" for k in classes])
ax.set_xlabel("selection coefficient  s  (per year)")
ax.set_title("Fitness per variant class, against two published methods")
ax.legend(loc="lower right")
ax.grid(axis="y", visible=False)
clean(ax)
fig.tight_layout(); fig.savefig(FIGS/"09_per_variant.png"); plt.close(fig)
print("  09_per_variant.png")


section("6. WHAT TO MAKE OF IT")
print("""The published methods disagree with each other by up to threefold on the
same gene, so 'correct' is not a single number here. An estimate that sits
between a cross-sectional and a longitudinal measurement of the same thing is
consistent with both, and one that sits outside both would be a signal that
something in the model is off.

The real test is not any single value but the ORDER: whether the model, given
only variant-restricted VAF distributions, recovers the same ranking that
direct regression on the raw data found in stage 1, and that the literature
reports. A model that reproduces an ordering it was not told about is doing
something more than curve fitting.""")
