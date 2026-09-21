#!/usr/bin/env python3
"""
One grid point, tailored to the Fabre cohort, emitting TWO summary vectors.

WHY A SECOND SWEEP EXISTS
Stage C's grid was built for the Watson table: ages 50/60/70/80, a floor of
VAF 0.0192, and a mutation-rate range guessed before anything was known. Stage F
scored the Fabre cohort against it and had to say so in its limitations —
"half of the 2.42 is ours rather than the data's", because the numerator came
from a ruler made for a different cohort.

This grid uses Fabre's own sampling: their baseline age spread, their 13-year
follow-up, their floor, their read depth.

WHY TWO VECTORS
The same simulated people can be measured the way a cross-sectional study
measures them (one draw, the shape of the VAF distribution) and the way a
longitudinal study does (repeated draws, one growth rate per clone). Recording
both from one simulation means the same model can be fitted to the same cohort
two ways, and the two answers compared. Stages E and F predict they will differ
by about 2.4; if our own estimator family reproduces that gap on real data, the
result stops depending on a grid borrowed from elsewhere.

Usage:
    sweep_fabre.py --s 0.10 --mu 2e-6 --people 5000 --seed 1 --out point.csv
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from model.population import simulate_cohort, sizes_to_vaf

QUANTILES = (0.10, 0.25, 0.50, 0.75, 0.90)

# --- the cohort this grid is shaped to ---------------------------------------
# Fabre's first-phase ages run 55 to 94 with a median of 70; the spectrum arm
# assigns each simulated person one of these, weighted to match.
CROSS_AGES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
CROSS_WEIGHTS = (0.06, 0.13, 0.20, 0.22, 0.18, 0.12, 0.06, 0.03)

# The trajectory arm: enrolled at 65, five draws over 13 years, as in stage H.
TRAJ_AGES = (65.0, 68.25, 71.5, 74.75, 78.0)

FLOOR = 0.002      # Fabre's usable floor after filtering
DEPTH = 1140       # their median sequencing depth


def observe(vaf, depth, rng):
    """What the sequencer reports: Binomial(depth, VAF)/depth."""
    out = np.zeros_like(vaf)
    m = vaf > 0
    out[m] = rng.binomial(depth, np.clip(vaf[m], 0.0, 1.0)) / depth
    return out


def quantile_row(values, prefix):
    row = {}
    for q in QUANTILES:
        row[f"{prefix}_q{int(q*100):02d}"] = (float(np.quantile(values, q))
                                              if values.size else float("nan"))
    return row


def cross_sectional(result, index, N, rng):
    """One draw per person, at an age drawn from the cohort's distribution."""
    people = result["sizes"][0].shape[0]
    assigned = rng.choice(len(CROSS_AGES), size=people, p=CROSS_WEIGHTS)
    kept = []
    per_person = np.zeros(people, dtype=np.int64)
    for k, age in enumerate(CROSS_AGES):
        who = np.flatnonzero(assigned == k)
        if who.size == 0:
            continue
        v = observe(sizes_to_vaf(result["sizes"][index[age]][who], N), DEPTH, rng)
        d = v >= FLOOR
        per_person[who] = d.sum(axis=1)
        kept.append(v[d])
    values = np.concatenate(kept) if kept else np.array([])
    row = {"cs_n_detected": int(values.size),
           "cs_carriers": float((per_person > 0).mean()),
           "cs_clones_per_person": float(per_person.mean())}
    row.update(quantile_row(values, "cs_vaf"))
    row["cs_vaf_mean_log"] = float(np.log(values).mean()) if values.size else float("nan")
    return row


def trajectories(result, index, N, rng):
    """
    Clones detectable at enrolment, followed for 13 years.

    The birth stamp is checked at every draw: slots are recycled when a clone
    dies, and following 'the clone in slot 7' without checking splices two
    different clones into one invented trajectory.
    """
    snaps = [observe(sizes_to_vaf(result["sizes"][index[a]], N), DEPTH, rng)
             for a in TRAJ_AGES]
    births = [result["births"][index[a]] for a in TRAJ_AGES]

    same = np.isfinite(births[0])
    for b in births[1:]:
        same &= np.isfinite(b) & (b == births[0])
    enrolled = same & (snaps[0] >= FLOOR)
    for s in snaps[1:]:
        enrolled &= s > 0

    idx = np.flatnonzero(enrolled.ravel())
    if idx.size < 20:
        row = {"tr_n": int(idx.size), "tr_median": float("nan"),
               "tr_early": float("nan"), "tr_late": float("nan"),
               "tr_decel": float("nan")}
        row.update(quantile_row(np.array([]), "tr_rate"))
        return row

    ages = np.asarray(TRAJ_AGES)
    logv = np.stack([np.log(s.ravel()[idx]) for s in snaps])       # (T, clones)
    x = ages - ages.mean()
    rates = (x[:, None] * logv).sum(axis=0) / (x * x).sum()

    half = len(ages) // 2 + 1
    xe = ages[:half] - ages[:half].mean()
    xl = ages[half - 1:] - ages[half - 1:].mean()
    early = (xe[:, None] * logv[:half]).sum(axis=0) / (xe * xe).sum()
    late = (xl[:, None] * logv[half - 1:]).sum(axis=0) / (xl * xl).sum()

    row = {"tr_n": int(idx.size), "tr_median": float(np.median(rates)),
           "tr_early": float(np.median(early)), "tr_late": float(np.median(late)),
           "tr_decel": float(np.median(early - late))}
    row.update(quantile_row(rates, "tr_rate"))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s", type=float, required=True)
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--N", type=int, default=100_000)
    ap.add_argument("--people", type=int, default=5_000)
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=str, required=True)
    a = ap.parse_args()

    ages = sorted(set(CROSS_AGES) | set(TRAJ_AGES))
    result = simulate_cohort(N=a.N, s=a.s, mu=a.mu, years=max(ages),
                             people=a.people, dt=a.dt, max_clones=96,
                             seed=a.seed, record_ages=ages)
    index = {float(age): i for i, age in enumerate(result["ages"])}
    rng = np.random.default_rng(a.seed + 7919)

    row = {"s": a.s, "mu": a.mu, "N": a.N, "seed": a.seed,
           "people": a.people, "floor": FLOOR, "depth": DEPTH}
    row.update(cross_sectional(result, index, a.N, rng))
    row.update(trajectories(result, index, a.N, rng))
    row["overflow"] = result["overflow"]

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    main()
