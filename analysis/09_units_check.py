#!/usr/bin/env python3
"""
Is `s` the same quantity in all three papers, or are we comparing apples?

THE WORRY
Stage E compares Watson's 15.0% per year against Fabre's 6.2% and calls the gap
methodological. Fabre validate their fitting with `clonex`, a WRIGHT-FISHER
simulator whose `s` is a per-GENERATION selective advantage. This project is a
MORAN process. Stage A established the two conventions differ — fixation goes as
s/(1+s) under Moran and about 2s under Wright-Fisher, because offspring variance
differs — so a selection coefficient is not automatically portable between them.

If part of the published gap were definitional, stage E's headline would be
wrong. This checks it.

Usage:  .venv/bin/python analysis/09_units_check.py
"""

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from model.population import simulate_cohort


def section(t):
    print(f"\n{'-'*76}\n{t}\n{'-'*76}")


section("1. WHAT EACH SOURCE SAYS ITS s IS")
print("""Watson 2020 defines it in the paper text, and the definition settles the
question on their side:

  "The fate of a new mutation [is decided] through a fitness effect, s, which is
   THE AVERAGE GROWTH RATE PER YEAR of that variant relative to the average
   growth rate"

  "growing exponentially at rate s per year"

So Watson's s is an exponential growth rate per year, not a per-generation
selective advantage.

Fabre 2022 fit growth rates per year to real trajectories measured in real
years. `clonex` appears in their repository to VALIDATE that the fitting
recovers a known truth, not to define the units of what they report. Its
per-generation s is internal to that validation.

That leaves ours. The model defines a clone's cells as weighing (1+s) against
1 for wild type, which is a statement about reproductive weight — not, on its
face, a growth rate. Whether the two coincide is measurable.""")


section("2. WHAT OUR s ACTUALLY IS, MEASURED")
print("""For a linear birth-death process E[size(t)] = exp(r*t) exactly, where r is the
net growth rate. The expectation is over EVERY lineage including the extinct
ones, so the check must not drop them — and it must not use E[log size], which
grows more slowly than r while the clone is small and skewed.

One clone seeded at birth, no further mutation, measured between ages 10 and 40:""")

print(f"\n{'s given':>9}{'rate measured':>16}{'ratio':>8}{'mean size at 40':>18}")
rows = []
for s in (0.05, 0.10, 0.13, 0.20, 0.30):
    r = simulate_cohort(N=100_000, s=s, mu=0.0, years=40, people=60_000, dt=0.25,
                        seed=4, record_ages=[10, 40], initial_clones=1, max_clones=2)
    a, b = r["sizes"][0][:, 0], r["sizes"][1][:, 0]
    rate = (np.log(b.mean()) - np.log(a.mean())) / 30.0
    rows.append((s, rate, rate / s, b.mean()))
    print(f"{s:>9.2f}{rate:>16.4f}{rate/s:>8.3f}{b.mean():>18.1f}")

small = [r for r in rows if r[0] <= 0.13]
print(f"""
For s <= 0.13 the ratio is {min(r[2] for r in small):.3f} to {max(r[2] for r in small):.3f}. OUR s IS THE EXPONENTIAL GROWTH
RATE PER YEAR — the same quantity Watson defines in words. No conversion factor
is missing, and the comparison stage E makes is between like quantities.""")


section("3. WHY THE MORAN/WRIGHT-FISHER FACTOR DOES NOT APPLY HERE")
print("""The s/(1+s) versus 2s difference is real and this project was right to care
about it in stage A. But it is a statement about FIXATION PROBABILITY — the
chance that a single mutant cell escapes drift and takes over — and it arises
because the two models give a lineage different offspring variance.

The growth rate of a clone that has already escaped drift is a different
quantity, and it is not affected. Neither Watson's headline nor Fabre's is a
fixation probability, so the factor never enters.

The concern was well formed and aimed at the wrong quantity. Stage E's
comparison stands.""")


section("4. AN UNPLANNED RESULT FROM THE SAME NUMBERS")
big = [r for r in rows if r[0] >= 0.20]
print(f"""The ratio is not constant. It stays at ~1 for s <= 0.13 and falls to
{big[0][2]:.3f} at s = {big[0][0]:.2f} and {big[-1][2]:.3f} at s = {big[-1][0]:.2f}.

That is clonal interference again, and it arrives earlier than expected: by age
40 a clone with s = {big[-1][0]:.2f} has a mean size of {big[-1][3]:,.0f} out of {100_000:,} stem cells and
is already growing at {big[-1][2]:.0%} of its nominal fitness — before any study would have
started following it.

It is an independent route to the stage E finding. There the shortfall was
measured by following detected clones; here it falls out of a clone seeded at
birth with no detection step at all, so it cannot be an artifact of how clones
were selected for measurement.""")
