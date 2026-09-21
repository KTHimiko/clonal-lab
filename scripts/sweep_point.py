#!/usr/bin/env python3
"""
Simulate one point of the parameter grid and emit summary statistics.

This is what a single cluster job runs. It deliberately does NOT compute the
distance to the observed data: the expensive part (simulation) belongs on the
cluster, the cheap part (comparison) belongs wherever the analysis happens.
Keeping them apart means the grid can be re-scored against different cohorts
or different distance functions without re-running anything.

Usage:
    sweep_point.py --s 0.10 --mu 2e-6 --people 5000 --seed 1 --out point.csv
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from model.population import simulate_cohort, sizes_to_vaf

QUANTILES = (0.10, 0.25, 0.50, 0.75, 0.90)


def summarise(vaf_matrix, limit):
    """
    Reduce one cohort to a handful of numbers.

    ABC compares summary statistics rather than raw distributions, because
    the raw output of a stochastic model never matches anything exactly. The
    statistics have to be few enough to compare and rich enough to
    distinguish parameter values.

    We keep the shape of the detected VAF distribution (quantiles) and how
    many clones there are (count). Stage B showed these respond to different
    parameters: the mutation rate moves the count, the fitness moves the
    shape — which is what makes them separately identifiable.
    """
    detected = vaf_matrix >= limit
    per_person = detected.sum(axis=1)
    values = vaf_matrix[detected]

    row = {
        "n_people": vaf_matrix.shape[0],
        "n_detected": int(detected.sum()),
        "clones_per_person": float(per_person.mean()),
        "carriers": float((per_person > 0).mean()),
    }
    for q in QUANTILES:
        row[f"vaf_q{int(q*100):02d}"] = float(np.quantile(values, q)) if values.size else float("nan")
    row["vaf_mean_log"] = float(np.log(values).mean()) if values.size else float("nan")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s", type=float, required=True, help="selection coefficient")
    ap.add_argument("--mu", type=float, required=True, help="driver mutation rate per cell per year")
    ap.add_argument("--N", type=int, default=50_000, help="stem cell population")
    ap.add_argument("--people", type=int, default=5_000)
    ap.add_argument("--ages", type=str, default="50,60,70,80")
    ap.add_argument("--limit", type=float, default=0.0192, help="VAF detection limit")
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=str, required=True)
    a = ap.parse_args()

    ages = [float(x) for x in a.ages.split(",")]

    result = simulate_cohort(N=a.N, s=a.s, mu=a.mu, years=max(ages),
                             people=a.people, dt=a.dt, max_clones=64,
                             seed=a.seed, record_ages=ages)

    rows = []
    for age, snapshot in zip(result["ages"], result["sizes"]):
        row = {"s": a.s, "mu": a.mu, "N": a.N, "age": float(age),
               "limit": a.limit, "seed": a.seed}
        row.update(summarise(sizes_to_vaf(snapshot, a.N), a.limit))
        row["overflow"] = result["overflow"]
        rows.append(row)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
