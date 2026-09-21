#!/usr/bin/env python3
"""
One batch of the mixture grid: a fitness DISTRIBUTION instead of a single s.

WHY THIS EXISTS
Stage J found the trajectory arm could not be fitted, and the diagnosis was not
the estimator: the observed per-clone growth rates span 0.200 from the 10th to
the 90th percentile, and a single-s model produces at most 0.132 anywhere on its
grid. Real clones differ from each other in fitness; our model's clones did not.

A spot check then showed a gamma mixture reaches 0.2015 — essentially exact —
but that every mixture wide enough to do so makes clones far too large in the
VAF spectrum. The two observables pull in opposite directions.

The obvious suspect is N, which was taken from the literature in stage C and
never fitted. More stem cells means the same clone is a smaller fraction of the
blood, which lowers VAF without touching growth rates — exactly the direction
needed. So N joins the grid.

BATCHED ON PURPOSE
Stage J ran 252 jobs of 7.5 seconds and spent 36 minutes on 3.5 CPU-hours: the
scheduler cost more than the work. Each job here sweeps a whole row of means,
which cuts the submission overhead by the length of that row.

Usage:
    sweep_mixture.py --N 100000 --shape 1.0 --seed 1 \
        --mean-list 0.04,0.06,0.08,0.10,0.12,0.14,0.16 --out batch.csv
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from model.population import simulate_cohort, sizes_to_vaf

QS = (0.10, 0.25, 0.50, 0.75, 0.90)
CROSS_AGES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
CROSS_WEIGHTS = (0.06, 0.13, 0.20, 0.22, 0.18, 0.12, 0.06, 0.03)
TRAJ_AGES = (65.0, 68.25, 71.5, 74.75, 78.0)
FLOOR, DEPTH = 0.002, 1140
CLASSES = 8


def mixture(mean, shape, k=CLASSES):
    """
    k equally weighted fitness classes at the midpoint quantiles of a gamma.

    The quantiles come from a large deterministic sample rather than from
    scipy's inverse CDF: the cluster nodes are a minimal cloud image with numpy
    and no scipy, and removing a dependency is cheaper and more durable than
    installing one on three machines. A million draws under a fixed seed puts
    these quantiles well inside the noise of the simulation they feed.
    """
    q = (np.arange(k) + 0.5) / k
    draws = np.random.default_rng(20260921).gamma(shape, mean / shape, 1_000_000)
    return list(np.round(np.quantile(draws, q), 5)), [1.0 / k] * k


def observe(v, rng):
    out = np.zeros_like(v)
    m = v > 0
    out[m] = rng.binomial(DEPTH, np.clip(v[m], 0.0, 1.0)) / DEPTH
    return out


def summarise(mean, shape, N, mu, people, dt, seed):
    s, w = mixture(mean, shape)
    ages = sorted(set(CROSS_AGES) | set(TRAJ_AGES))
    r = simulate_cohort(N=N, s=s, s_weights=w, mu=mu, years=max(ages),
                        people=people, dt=dt, max_clones=96, seed=seed,
                        record_ages=ages)
    idx = {float(a): i for i, a in enumerate(r["ages"])}
    rng = np.random.default_rng(seed + 7919)

    row = {"mean": mean, "shape": shape, "N": N, "mu": mu,
           "people": people, "seed": seed, "floor": FLOOR, "depth": DEPTH}

    # --- spectrum arm: one draw per person, age from the cohort distribution
    assigned = rng.choice(len(CROSS_AGES), size=people, p=CROSS_WEIGHTS)
    kept, per_person = [], np.zeros(people, dtype=np.int64)
    for k, age in enumerate(CROSS_AGES):
        who = np.flatnonzero(assigned == k)
        if who.size == 0:
            continue
        v = observe(sizes_to_vaf(r["sizes"][idx[age]][who], N), rng)
        d = v >= FLOOR
        per_person[who] = d.sum(axis=1)
        kept.append(v[d])
    vals = np.concatenate(kept) if kept else np.array([])
    row["cs_n"] = int(vals.size)
    row["cs_carriers"] = float((per_person > 0).mean())
    for q in QS:
        row[f"cs_vaf_q{int(q*100):02d}"] = float(np.quantile(vals, q)) if vals.size else float("nan")

    # --- trajectory arm: enrolled at 65, five draws to 78
    snaps = [observe(sizes_to_vaf(r["sizes"][idx[a]], N), rng) for a in TRAJ_AGES]
    births = [r["births"][idx[a]] for a in TRAJ_AGES]
    same = np.isfinite(births[0])
    for b in births[1:]:
        same &= np.isfinite(b) & (b == births[0])
    sel = same & (snaps[0] >= FLOOR)
    for sn in snaps[1:]:
        sel &= sn > 0
    i = np.flatnonzero(sel.ravel())
    row["tr_n"] = int(i.size)
    if i.size >= 30:
        a = np.asarray(TRAJ_AGES)
        x = a - a.mean()
        logv = np.stack([np.log(sn.ravel()[i]) for sn in snaps])
        rates = (x[:, None] * logv).sum(axis=0) / (x * x).sum()
        for q in QS:
            row[f"tr_rate_q{int(q*100):02d}"] = float(np.quantile(rates, q))
        row["tr_negative"] = float((rates < 0).mean())
        row["tr_width"] = float(np.quantile(rates, 0.9) - np.quantile(rates, 0.1))
    else:
        for q in QS:
            row[f"tr_rate_q{int(q*100):02d}"] = float("nan")
        row["tr_negative"] = float("nan")
        row["tr_width"] = float("nan")
    row["overflow"] = r["overflow"]
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mean-list", required=True)
    ap.add_argument("--shape", type=float, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--mu", type=float, default=2e-6)
    ap.add_argument("--people", type=int, default=6000)
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [summarise(float(m), a.shape, a.N, a.mu, a.people, a.dt, a.seed)
            for m in a.mean_list.split(",")]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
