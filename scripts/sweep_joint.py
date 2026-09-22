#!/usr/bin/env python3
"""
Stage O — one batch of the joint grid: clone population AND hazard together.

Every earlier run of the hazard fixed the clone population at values fitted with
no mechanism at all. That understates the mechanism, because a hazard that
removes clones changes which population fits best. This frees them together.

Emits five numbers per point. Two are fitted against (the VAF spectrum and the
growth-rate distribution); three are consequences the fit does not optimise for
and are therefore the actual test — the declining fraction, the within-person
correlation, and whether decline persists within a clone.

Plan and decision rule: analysis/PREREGISTRATION_STAGE_O.md

Usage:
    sweep_joint.py --N 200000 --shape 2.0 --rate 0.03 --from-age 55 --seed 1 \
        --mean-list 0.05,0.07,0.09,0.11 --out batch.csv
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
TRAJ = (65.0, 68.25, 71.5, 74.75, 78.0)
ALL = sorted(set(CROSS_AGES) | set(TRAJ))
FLOOR, DEPTH, DELTA, CLASSES = 0.002, 1140, 0.20, 8


def mixture(mean, shape, k=CLASSES):
    """Gamma quantiles from a deterministic sample: the nodes have no scipy."""
    q = (np.arange(k) + 0.5) / k
    draws = np.random.default_rng(20260921).gamma(shape, mean / shape, 1_000_000)
    return list(np.round(np.quantile(draws, q), 5)), [1.0 / k] * k


def fit_masked(logv, mask, ages):
    """Slope per clone over the draws that cleared the floor, `ages` explicit."""
    ages = np.asarray(ages, dtype=float)
    n = mask.sum(axis=0)
    safe = np.where(n > 0, n, 1)
    w = mask.astype(float)
    xbar = (w * ages[:, None]).sum(axis=0) / safe
    xc = (ages[:, None] - xbar) * w
    y = np.where(mask, logv, 0.0)
    num = (xc * y).sum(axis=0) - xc.sum(axis=0) * (y.sum(axis=0) / safe)
    den = (xc * (ages[:, None] - xbar) * w).sum(axis=0)
    ok = (n >= 2) & (den > 0)
    return np.where(ok, num / np.where(ok, den, 1.0), np.nan)


def summarise(mean, shape, N, rate, from_age, mu, people, dt, seed):
    s, w = mixture(mean, shape)
    r = simulate_cohort(N=N, s=s, s_weights=w, mu=mu, years=max(ALL),
                        people=people, dt=dt, max_clones=96, seed=seed,
                        record_ages=ALL, switch_rate=rate,
                        switch_delta=DELTA, switch_from=from_age)
    idx = {float(a): i for i, a in enumerate(r["ages"])}
    rng = np.random.default_rng(seed + 7919)

    def observe(age):
        v = sizes_to_vaf(r["sizes"][idx[age]], N)
        o = np.zeros_like(v)
        m = v > 0
        o[m] = rng.binomial(DEPTH, np.clip(v[m], 0.0, 1.0)) / DEPTH
        return o

    snaps_all = {a: observe(a) for a in ALL}
    row = {"mean": mean, "shape": shape, "N": N, "rate": rate,
           "from_age": from_age, "delta": DELTA, "mu": mu,
           "people": people, "seed": seed}

    assigned = rng.choice(len(CROSS_AGES), size=people, p=CROSS_WEIGHTS)
    kept = []
    for j, age in enumerate(CROSS_AGES):
        who = np.flatnonzero(assigned == j)
        if who.size:
            v = snaps_all[age][who]
            kept.append(v[v >= FLOOR])
    vals = np.concatenate(kept) if kept else np.array([])
    row["cs_n"] = int(vals.size)
    for q in QS:
        row[f"cs_vaf_q{int(q*100):02d}"] = float(np.quantile(vals, q)) if vals.size else float("nan")

    snaps = np.stack([snaps_all[a] for a in TRAJ])
    births = [r["births"][idx[a]] for a in TRAJ]
    same = np.isfinite(births[0])
    for b in births[1:]:
        same &= np.isfinite(b) & (b == births[0])
    above = snaps >= FLOOR
    enrolled = same & above[0] & (above.sum(axis=0) >= 3)

    nan5 = {f"tr_rate_q{int(q*100):02d}": float("nan") for q in QS}
    if enrolled.sum() < 30:
        row.update(nan5)
        row.update(tr_n=0, tr_negative=float("nan"), tr_icc=float("nan"),
                   tr_pff=float("nan"), tr_pfr=float("nan"), overflow=r["overflow"])
        return row

    flat = enrolled.ravel()
    logv = np.log(np.where(above, snaps, 1.0)).reshape(len(TRAJ), -1)[:, flat]
    mask = above.reshape(len(TRAJ), -1)[:, flat]
    rates = fit_masked(logv, mask, TRAJ)
    ok = np.isfinite(rates)
    rates = rates[ok]

    row["tr_n"] = int(rates.size)
    for q in QS:
        row[f"tr_rate_q{int(q*100):02d}"] = float(np.quantile(rates, q))
    row["tr_negative"] = float((rates < 0).mean())

    half = len(TRAJ) // 2 + 1
    a = np.asarray(TRAJ)
    early = fit_masked(logv[:half][:, ok], mask[:half][:, ok], a[:half])
    late = fit_masked(logv[half - 1:][:, ok], mask[half - 1:][:, ok], a[half - 1:])
    good = np.isfinite(early) & np.isfinite(late)
    fell, rose = good & (early < 0), good & (early >= 0)
    row["tr_pff"] = float((late[fell] < 0).mean()) if fell.any() else float("nan")
    row["tr_pfr"] = float((late[rose] < 0).mean()) if rose.any() else float("nan")

    pid = np.broadcast_to(np.arange(people)[:, None], enrolled.shape).ravel()[flat][ok]
    uniq, counts = np.unique(pid, return_counts=True)
    keep = np.isin(pid, uniq[counts >= 2])
    row["tr_icc"] = float("nan")
    if keep.sum() > 50:
        pg, rg = pid[keep], rates[keep]
        groups = [rg[pg == q] for q in np.unique(pg)]
        k = np.mean([len(g) for g in groups])
        msb = np.var([g.mean() for g in groups], ddof=1) * k
        msw = np.mean([g.var(ddof=1) for g in groups])
        row["tr_icc"] = float((msb - msw) / (msb + (k - 1) * msw))
    row["overflow"] = r["overflow"]
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mean-list", required=True)
    ap.add_argument("--shape", type=float, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--from-age", type=float, required=True)
    ap.add_argument("--mu", type=float, default=2e-6)
    ap.add_argument("--people", type=int, default=7000)
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [summarise(float(m), a.shape, a.N, a.rate, a.from_age, a.mu,
                      a.people, a.dt, a.seed)
            for m in a.mean_list.split(",")]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
