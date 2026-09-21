#!/usr/bin/env python3
"""
Stage B — many clones, continuous mutation influx, and a VAF distribution.

Stage A had an exact answer to check against. This one does not, so it is
checked in two other ways: by reduction (with no influx and a single clone it
must return stage A's formula), and by confrontation (the distribution it
produces is placed next to the observed one).

Usage:  .venv/bin/python analysis/03_multiclone.py
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
from model.population import simulate_cohort, sizes_to_vaf

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


def log_ticks(ax, values=(0.02, 0.05, 0.1, 0.2, 0.5)):
    """Explicit ticks on a log axis.

    Matplotlib's automatic minor labels collide on a narrow log range — the
    default put five overlapping labels under each panel. Naming the ticks
    ourselves is the fix.
    """
    ax.set_xticks(values)
    ax.set_xticklabels([f"{v:g}" for v in values])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())

def section(t):
    print(f"\n{'-'*74}\n{t}\n{'-'*74}")


N_CELLS = 50_000          # stem cells; published estimates span 50k-200k
COHORT = "Coombs2017"     # the cohort we compare against
LIMIT = 0.0192            # its detection limit, measured in stage 1


# ============================================== check 1: reduction to stage A
section("CHECK 1 - REDUCTION TO THE SINGLE-CLONE CASE")
print("With no mutation influx and one seeded clone, the model must return the")
print("exact Moran result s/(1+s). If it does not, nothing downstream is worth")
print("reading.\n")
print(f"{'s':>6}  {'expected':>10}  {'simulated':>10}  {'+/-':>8}  {'sigmas':>7}")

failures = 0
for s in [0.05, 0.10, 0.20, 0.50]:
    r = simulate_cohort(N=N_CELLS, s=s, mu=0.0, years=150, people=20_000,
                        dt=0.5, max_clones=2, initial_clones=1, seed=11)
    surv = (r["sizes"][0][:, 0] > 0).mean()
    err = np.sqrt(surv * (1 - surv) / 20_000)
    exp = s / (1 + s)
    z = abs(surv - exp) / err
    failures += z >= 3
    print(f"{s:>6.2f}  {exp:>10.5f}  {surv:>10.5f}  {err:>8.5f}  {z:>6.2f}"
          f"  {'ok' if z < 3 else 'FAILED'}")

if failures:
    print("\nReduction check failed. Stopping.")
    sys.exit(1)


# ======================================================= the observed target
section("CHECK 2 - THE OBSERVED DISTRIBUTION")

raw = pd.read_csv(CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
raw["age"] = pd.to_numeric(raw["age"], errors="coerce")
raw["VAF"] = pd.to_numeric(raw["VAF"], errors="coerce")
obs = raw.dropna(subset=["age", "VAF"])
obs = obs[obs.study == COHORT]

print(f"cohort: {COHORT}")
print(f"  variants:       {len(obs)}")
print(f"  people's ages:  {obs.age.min():.0f} to {obs.age.max():.0f} "
      f"(median {obs.age.median():.0f})")
print(f"  detection limit: VAF >= {LIMIT}")
print(f"  median VAF:     {obs.VAF.median():.4f}")
print("\nWe simulate people at the same ages and apply the same limit. Without")
print("that truncation the model would be judged on clones the instrument")
print("could never have seen.")


# ================================================= the simulation, by mutation rate
section("CHECK 3 - HOW MANY CLONES THE MODEL PRODUCES")
print("The mutation rate sets how many clones appear; the fitness sets how big")
print("they get. Here we hold s fixed and vary the rate, to see which order of")
print("magnitude puts the model in the same range as the data.\n")

ages = np.array(sorted(obs.age.round().unique()))
print(f"{'mu (per cell/yr)':>18}  {'clones/person':>14}  {'median VAF':>12}")

results = {}
for mu in [5e-7, 2e-6, 5e-6]:
    r = simulate_cohort(N=N_CELLS, s=0.10, mu=mu, years=float(ages.max()),
                        people=2_000, dt=0.5, max_clones=64, seed=3,
                        record_ages=[float(obs.age.median())])
    vaf = sizes_to_vaf(r["sizes"][0], N_CELLS)
    seen = vaf >= LIMIT
    per_person = seen.sum(axis=1).mean()
    med = np.median(vaf[seen]) if seen.any() else float("nan")
    results[mu] = (per_person, med, vaf[seen])
    print(f"{mu:>18.1e}  {per_person:>14.3f}  {med:>12.4f}")

print(f"\nobserved, same cohort and limit:  median VAF {obs.VAF.median():.4f}")
print("\nThe clone count per person is not directly comparable: the data holds")
print("only people who had a detectable clone, not everyone screened. That")
print("missing denominator was the limitation found in stage 1, and it is why")
print("calibration has to run on the VAF distribution rather than on counts.")


# ================================================================= figures
section("FIGURES")

best_mu = 2e-6
sim_vaf = results[best_mu][2]

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))

# left: the distributions side by side, as cumulative curves
ax = axes[0]
for values, colour, label in ((np.sort(sim_vaf), S2, f"simulated (s=0.10, mu={best_mu:.0e})"),
                              (np.sort(obs.VAF.values), S1, f"observed ({COHORT}, n={len(obs)})")):
    if len(values):
        ax.step(values, np.arange(1, len(values)+1)/len(values), where="post",
                color=colour, linewidth=2, label=label)
ax.set_xscale("log")
ax.set_xlabel("variant allele frequency (VAF, log scale)")
ax.set_ylabel("cumulative proportion")
ax.set_title("Above the detection limit, shapes are comparable", fontsize=11)
ax.legend(loc="lower right")
clean(ax); log_ticks(ax)

# right: how the model's VAF distribution shifts with age
ax = axes[1]
snap_ages = [40.0, 60.0, 80.0]
r = simulate_cohort(N=N_CELLS, s=0.10, mu=best_mu, years=80.0, people=4_000,
                    dt=0.5, max_clones=64, seed=5, record_ages=snap_ages)
for age, snap, colour in zip(r["ages"], r["sizes"], (S1, S3, S2)):
    v = sizes_to_vaf(snap, N_CELLS)
    v = v[v >= LIMIT]
    if len(v):
        ax.step(np.sort(v), np.arange(1, len(v)+1)/len(v), where="post",
                color=colour, linewidth=2, label=f"age {age:.0f}  (n={len(v)})")
ax.set_xscale("log")
ax.set_xlabel("variant allele frequency (VAF, log scale)")
ax.set_ylabel("cumulative proportion")
ax.set_title("The distribution shifts right with age", fontsize=11)
ax.legend(loc="lower right")
clean(ax); log_ticks(ax)

fig.suptitle("Stage B: many clones, continuous mutation influx",
             x=0.012, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(FIGS/"07_vaf_distribution.png"); plt.close(fig)
print("  07_vaf_distribution.png")


# ================================================================== summary
section("WHERE THIS LEAVES US")
print("""The model now produces the right KIND of object: a distribution of
variant allele frequencies that shifts with age, truncated the way a real
instrument truncates it.

What it does not yet do is match. Fitness is still a single number shared by
every clone, and stage 1 showed that fitness is per variant — TET2 grows
about twice as fast as DNMT3A, and within DNMT3A the R882 hotspot grows
faster than the rest.

Stage C is the search: sweep the mutation rate and the fitness over a grid,
run thousands of cohorts, and keep the parameter pairs whose distribution
lands close to the observed one. That is the load the cluster exists for.""")
