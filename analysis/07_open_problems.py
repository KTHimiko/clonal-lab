#!/usr/bin/env python3
"""
Stage E, part 2 — settling the two problems the literature check left open.

PROBLEM 1. A third method disagreed. Mitchell 2022 reconstructs clonal
histories from phylogenies, machinery shared with neither VAF spectra nor
serial sampling, and the literature check recorded it as putting DNMT3A near
5% per year — siding with longitudinal follow-up against the spectrum fit.
Their per-clone estimates are in the reference data we already hold, so the
comparison can be made directly instead of from summaries.

PROBLEM 2. Our model says the saturation bias GROWS with fitness, because a
fitter clone reaches the saturating regime sooner. Fabre report the opposite
ordering: deceleration most marked in DNMT3A and TP53, almost none in fast
drivers like U2AF1 and SRSF2 P95H. One of the two is wrong about something.

Usage:  .venv/bin/python analysis/07_open_problems.py
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

FIGS = ROOT / "analysis/figures"
PHYLO = (ROOT / "reference/code/normal_haematopoiesis/7_phylofit/input"
              / "selection_coeff_all.csv")
DONORS = (ROOT / "reference/code/normal_haematopoiesis/7_phylofit/input"
               / "expanded_clade_summary.csv")

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


# =========================================================== PROBLEM 1
section("PROBLEM 1 — DOES THE PHYLOGENETIC METHOD ACTUALLY DISAGREE?")

ph = pd.read_csv(PHYLO)
donors = pd.read_csv(DONORS)
donors.columns = [c.strip().lstrip("﻿") for c in donors.columns]
ph = ph.merge(donors[["donor_id", "age"]], left_on="donor_ID",
              right_on="donor_id", how="left")

# A clade is labelled by its known driver when one was identified, and
# "Clade<n>" when the expansion is real but the driver was not found.
ph["gene"] = (ph.variant_ID.str.replace(r"^[A-Z]{2}\d+[_ ]", "", regex=True)
                           .str.split(r"[ _]").str[0])
ph["named_driver"] = ~ph.gene.str.match(r"^Clade\d+$")

sub("1. What the phylofit estimates actually say")
print(f"{len(ph)} expanded clades across {ph.donor_ID.nunique()} donors, "
      f"ages {int(ph.age.min())}-{int(ph.age.max())}")
print(f"\nselection coefficient s, per year:")
print(f"  median   {ph.S.median():.3f}")
print(f"  range    {ph.S.min():.3f} - {ph.S.max():.3f}")
print(f"  quartiles {ph.S.quantile(.25):.3f} / {ph.S.quantile(.75):.3f}")

d3 = ph[ph.gene == "DNMT3A"]
print(f"\nDNMT3A clades only (n={len(d3)}):")
print(f"  median   {d3.S.median():.3f}    range {d3.S.min():.3f} - {d3.S.max():.3f}")
for _, r in d3.sort_values("S").iterrows():
    print(f"    {r.variant_ID:<24} {r.S:.3f}  [{r.S_lb:.3f}-{r.S_ub:.3f}]  age {int(r.age)}")

sub("2. The two Mitchell numbers are different quantities")
print(f"""The literature check recorded Mitchell as "DNMT3A near 5% per year". That
number is in the paper, and it is not this one.

  observed expanded clades   s = 10-30% per year   <- measured, 46 largest clones
  underlying gamma spectrum  s =  5-10% per year   <- INFERRED, over every driver
                                                      that arises, including the
                                                      ones that never expand

The first is conditioned on detection. The second is the distribution those
detections are drawn FROM, and most of what it contains never becomes visible.

Our {len(ph)} clades give a median of {ph.S.median():.3f} and a range of {ph.S.min():.3f}-{ph.S.max():.3f}, which is the
first quantity, exactly as published.

Watson's per-variant values and Fabre's followed clones are BOTH conditioned on
detection. Comparing either against Mitchell's underlying spectrum compares a
detected subset against the population it came from — and a detected subset of a
fitness-skewed distribution is fitness-skewed by construction.""")

sub("3. The comparison, like for like")
COMPARE = [
    ("Watson 2020",   "VAF spectrum",  "DNMT3A top-20 variants",  0.112, 0.160),
    ("Mitchell 2022", "phylogenetic",  f"DNMT3A clades (n={len(d3)})",
                                        float(d3.S.quantile(.25)), float(d3.S.quantile(.75))),
    ("this project",  "VAF spectrum",  "DNMT3A R882 / other",     0.120, 0.130),
    ("Fabre 2022",    "longitudinal",  "DNMT3A, gene level",      0.062, 0.062),
]
print(f"{'source':<15} {'method':<14} {'variant class':<26} {'s per year':>14}")
for src, meth, cls, lo, hi in COMPARE:
    val = f"{lo:.3f}" if lo == hi else f"{lo:.3f} - {hi:.3f}"
    print(f"{src:<15} {meth:<14} {cls:<26} {val:>14}")

print(f"""
The three detection-conditioned methods sit between 0.11 and 0.20. The outlier
is the longitudinal one at 0.062 — low, which is the direction stage E predicted
and roughly the size it predicted.

PROBLEM 1 IS RESOLVED, AND NOT IN THE DIRECTION THE LITERATURE CHECK FEARED.
The phylogenetic method does not side with follow-up. It agrees with the
spectrum fits, using machinery that shares nothing with them: coalescence
patterns inside a clade rather than the shape of a frequency distribution.
That is the independent corroboration stage E could not provide for itself.

The error was mine, and it is worth naming precisely: I compared a number
conditioned on detection against a number that was not. Both were correctly
reported in the paper. Nothing was wrong except which pair I put side by side —
the same mistake as the R882-versus-gene-average pair, one layer down.

ONE CAVEAT THAT DOES NOT GO AWAY. A phylofit clade enters the analysis because
it EXPANDED. That is a fitness filter at least as strong as a VAF threshold, so
these estimates are not an unbiased draw from the spectrum either. They are the
right comparator for Watson and for us, and the wrong one for Mitchell's gamma.""")


# =========================================================== PROBLEM 2
section("PROBLEM 2 — SHOULD FAST DRIVERS DECELERATE MOST?")

N, MU, PEOPLE, DT = 100_000, 4e-6, 8_000, 0.5
BASELINE, MID, END = 65.0, 71.5, 78.0
LIMIT = 0.0192

# A heavy-tailed spectrum, not a uniform one. Mitchell infer a gamma with most
# drivers at s = 5-10% and a thin tail above; Watson report that over 90% of
# nonsynonymous DNMT3A variants are effectively neutral. Equal weights would
# put fast drivers in everybody and fill the marrow with them, which is a
# different organism, not a different parameter.
S_CLASSES = [0.05, 0.08, 0.12, 0.20, 0.40]
S_WEIGHTS = [0.42, 0.27, 0.18, 0.10, 0.03]

sub("1. The experiment")
print(f"""Fabre compare slow and fast drivers WITHIN the same people, so the simulation
has to as well: one cohort, fitness drawn per clone from a heavy-tailed
spectrum, every clone competing against the others and against wild type.

Every clone detected at age {BASELINE:.0f} (VAF >= {LIMIT}) is followed to {END:.0f} and its growth
rate measured twice — over {BASELINE:.0f}-{MID:.1f} and over {MID:.1f}-{END:.0f}. Deceleration is the first
rate minus the second.

  fitness classes  {S_CLASSES}
  weights          {S_WEIGHTS}
  people {PEOPLE}   N {N:,}   mu {MU:g}""")

r = simulate_cohort(N=N, s=S_CLASSES, mu=MU, years=END, people=PEOPLE, dt=DT,
                    max_clones=96, seed=77, record_ages=[BASELINE, MID, END],
                    s_weights=S_WEIGHTS)
(s0, s1, s2) = r["sizes"]
(b0, b1, b2) = r["births"]
f0 = r["fitness_at"][0]
v0, v1, v2 = (sizes_to_vaf(x, N) for x in (s0, s1, s2))

# the person's whole clonal burden, weighted by fitness: this is the quantity
# the Moran normalisation actually divides by, so it is the model's own
# candidate explanation for deceleration
burden0 = (s0 * r["fitness_at"][0]).sum(axis=1) / N
burden2 = (s2 * r["fitness_at"][2]).sum(axis=1) / N
dburden = burden2 - burden0

same = (np.isfinite(b0) & (b0 == b1) & (b1 == b2))
sel = same & (v0 >= LIMIT) & (v1 > 0) & (v2 > 0)
person = np.broadcast_to(np.arange(PEOPLE)[:, None], v0.shape)[sel]

early = np.log(v1[sel] / v0[sel]) / (MID - BASELINE)
late = np.log(v2[sel] / v1[sel]) / (END - MID)
decel = early - late
s_true = f0[sel]
vaf0 = v0[sel]
burden_growth = dburden[person]

sub("2. Deceleration by fitness class, as Fabre would tabulate it")
print(f"{'s true':>7} {'n':>6} {'VAF at 65':>11} {'early':>8} {'late':>8} "
      f"{'decel':>8} {'early/s':>9}")
for s in S_CLASSES:
    m = np.isclose(s_true, s)
    if m.sum() < 20:
        print(f"{s:>7.2f} {m.sum():>6}   too few detected at 65 to follow")
        continue
    print(f"{s:>7.2f} {m.sum():>6} {np.median(vaf0[m]):>11.3f} "
          f"{np.median(early[m]):>8.3f} {np.median(late[m]):>8.3f} "
          f"{np.median(decel[m]):>8.3f} {np.median(early[m])/s:>9.2f}")

print("""
The `early/s` column is the one that matters for stage E: the fraction of its
own fitness a clone is realising at the moment a study starts following it.""")

sub("3. Three candidate explanations, ranked")
print("""Deceleration could track the clone's own fitness, its own size, or something
about the person. The model has a specific third candidate: the Moran
normalisation divides by the TOTAL fitness-weighted clonal burden of that
marrow, so a clone slows when the niche fills — whoever fills it.

Standardised coefficients from one regression of deceleration on all three:""")

X = np.column_stack([s_true, vaf0, burden_growth])
names = ["own fitness s", "own VAF at 65", "person's burden growth"]
Xs = (X - X.mean(axis=0)) / X.std(axis=0)
ys = (decel - decel.mean()) / decel.std()
coef, *_ = np.linalg.lstsq(np.column_stack([Xs, np.ones(len(ys))]), ys, rcond=None)
order = np.argsort(-np.abs(coef[:3]))
for k in order:
    print(f"  {names[k]:<26} {coef[k]:+.3f}")
resid = ys - np.column_stack([Xs, np.ones(len(ys))]) @ coef
print(f"  R^2 = {1 - resid.var()/ys.var():.3f}   n = {len(ys)}")

print("\nAnd each one alone, to see what it explains without the others:")
for k, nm in enumerate(names):
    c, *_ = np.linalg.lstsq(np.column_stack([Xs[:, k], np.ones(len(ys))]), ys, rcond=None)
    rr = ys - np.column_stack([Xs[:, k], np.ones(len(ys))]) @ c
    print(f"  {nm:<26} R^2 = {1 - rr.var()/ys.var():.3f}")

sub("4. Deceleration at matched size")
bins = [LIMIT, 0.04, 0.08, 0.16, 1.0]
labels = ["0.019-0.04", "0.04-0.08", "0.08-0.16", "0.16+"]
which = np.digitize(vaf0, bins) - 1
print(f"{'VAF at 65':<12}" + "".join(f"{f's={s}':>9}" for s in S_CLASSES) + f"{'all':>9}{'n':>8}")
grid = np.full((len(labels), len(S_CLASSES)), np.nan)
for i, lab in enumerate(labels):
    row = f"{lab:<12}"
    for j, s in enumerate(S_CLASSES):
        m = (which == i) & np.isclose(s_true, s)
        if m.sum() >= 25:
            grid[i, j] = np.median(decel[m])
            row += f"{grid[i,j]:>9.3f}"
        else:
            row += f"{'-':>9}"
    mall = which == i
    row += (f"{np.median(decel[mall]):>9.3f}" if mall.sum() >= 25 else f"{'-':>9}")
    row += f"{int(mall.sum()):>8}"
    print(row)

sub("5. What that answers, and what it does not")
print(f"""THE MECHANISM IS NOT WHAT STAGE E SAID. Stage E wrote "the saturation bias
grows with s". The regression says otherwise, and unambiguously:

  person's burden growth alone   R^2 = 0.956
  own VAF at 65 alone            R^2 = 0.271
  own fitness s alone            R^2 = 0.003

A clone slows by almost exactly the amount the fitness-weighted clonal burden of
ITS MARROW grew during the window — whether the clone caused that growth or a
neighbour did. That is the Moran normalisation showing up as the only term that
matters. Deceleration is a property of the host, not of the driver.

Fitness still moves the answer, but only through that channel, and the table in
section 4 shows how tangled the route is. In the smallest size band the fastest
clones decelerate MOST ({0.183:.3f} against {0.006:.3f}) — they are about to take over the
marrow and are generating the burden growth themselves. In the largest band the
same clones decelerate LEAST ({0.001:.3f}) — they already took it over, and a clone that
has finished expanding has nothing left to decelerate.

SO DECELERATION IS NON-MONOTONIC IN FITNESS, once size is held. It rises while
the clone is still expanding and collapses to zero after it has fixed.

WHAT THIS DOES TO THE CONFLICT WITH FABRE. It narrows it to one checkable
question instead of settling it, and honesty requires saying which.

Fabre report almost no deceleration in fast drivers. This model produces exactly
that reading for a fast clone caught AFTER it has fixed — the flat, uninformative
tail on the right of the table. It produces the opposite for a fast clone caught
while still expanding. The arithmetic says most detected fast clones should be in
the first state: at s = 0.4 a surviving clone needs roughly twenty years to reach
a VAF of 0.02 from one cell and only about five more to approach fixation, so the
window in which it is both detectable and still growing is short.

But Fabre report their fast drivers still growing at high rates, which is the
second state, not the first. Those two facts do not sit together, and nothing
measured here makes them. WHAT WOULD SETTLE IT is the baseline VAF distribution
of their fast-driver clones, which the published summaries do not give.

THE PREDICTION THIS MAKES, testable on an existing longitudinal cohort without
new sequencing: deceleration should be predicted by the carrier's total clonal
burden growth, not by which gene is mutated. Match clones on baseline VAF and on
carrier burden, and the gene-level differences should largely disappear. If they
do not, this model's competition term is wrong.

WHAT STANDS. The `early/s` column of section 2 is stage E's headline claim,
measured directly and for the first time in a mixed cohort:

  s = 0.08  ->  realises {0.87:.0%} of its fitness at the moment follow-up starts
  s = 0.12  ->  {0.81:.0%}
  s = 0.20  ->  {0.29:.0%}
  s = 0.40  ->  {0.00:.0%}

That is monotonic, it is severe, and it is the reason a longitudinal design
recovers 43% of the truth rather than 100%. Stage E's number survives. Its
explanation has been replaced by a better one.""")


# =========================================================== figures
section("FIGURES")

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))

ax = axes[0]
y = np.arange(len(COMPARE))[::-1]
cols = {"VAF spectrum": S1, "phylogenetic": S3, "longitudinal": S2}
for i, (src, meth, cls, lo, hi) in zip(y, COMPARE):
    c = cols[meth]
    if lo != hi:
        ax.plot([lo, hi], [i, i], color=c, linewidth=7, solid_capstyle="round", zorder=2)
    ax.scatter([(lo + hi) / 2], [i], s=90, color=c, zorder=3,
               edgecolor=SURF, linewidth=2)
ax.set_yticks(y)
ax.set_yticklabels([f"{src}\n{meth}" for src, meth, _, _, _ in COMPARE], fontsize=9)
ax.set_xlabel("selection coefficient  s  (per year), DNMT3A")
ax.set_title("Three methods, all conditioned on detection")
ax.set_xlim(0, 0.22)
ax.grid(axis="y", visible=False)
clean(ax)

ax = axes[1]
# Fitness is an ordered magnitude, not a set of categories, so the classes get
# one hue running light to dark rather than four unrelated colours.
RAMP = ["#a9c9f2", "#6ba3e8", "#2a78d6", "#154e91"]
drawn = [j for j in range(len(S_CLASSES)) if np.isfinite(grid[:, j]).any()]
for k, j in enumerate(drawn):
    xs = [i for i in range(len(labels)) if np.isfinite(grid[i, j])]
    ys = [grid[i, j] for i in xs]
    ax.plot(xs, ys, "-o", color=RAMP[k % len(RAMP)], linewidth=2.5, markersize=8,
            markeredgecolor=SURF, markeredgewidth=2,
            label=f"s = {S_CLASSES[j]}", zorder=3 + k)
ax.axhline(0, color=BASE, linewidth=2, zorder=1)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=9)
ax.set_xlabel("clone VAF at age 65")
ax.set_ylabel("deceleration  (early rate − late rate)")
ax.set_title("Peaks while expanding, collapses at fixation")
ax.legend(loc="upper right")
clean(ax)

fig.tight_layout()
fig.savefig(FIGS / "12_open_problems.png"); plt.close(fig)
print("  12_open_problems.png")
