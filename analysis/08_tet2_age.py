#!/usr/bin/env python3
"""
The last open problem — TET2 and age-dependent fitness.

THE OBSERVATION
Fabre report that "age was a significant factor specifically for TET2-mutant
clones, which grew faster in older individuals." Our model has no age-dependent
fitness at all, and its competition term makes every clone slow down with age,
so this runs opposite to everything the simulator does.

THE MECHANISM IS NOT THE OBVIOUS ONE
Mouse work attributes the effect mainly to "the aging-associated reduction in
fitness of aged competitor non-mutant HSCs" — the TET2 clone does not speed up,
the wild type slows down. In a Moran process fitness is relative, so those are
different models with different consequences:

  s_ramp       the mutation itself gets better with host age.
               Lifts only the classes it is given for.
  wt_decline   the wild type gets worse with host age.
               Lifts EVERY clone at once, whatever its driver.

That difference is testable, and the test does not need new data: if the wild
type were degrading, DNMT3A clones would accelerate alongside TET2 ones.

Usage:  .venv/bin/python analysis/08_tet2_age.py
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
OBS_CSV = (ROOT / "reference/data/watson2020/Maximum_likelihood_estimations"
                / "Maximum likelihood estimations - data files"
                / "all_studies_trimmed_all_genes.csv")
ELIGIBLE = ["Young2019", "Young2016", "Acuna2017", "McKerrel2015",
            "Desai2018", "Coombs2017"]

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


N, MU, DT, LIMIT = 100_000, 2e-6, 0.5, 0.0192
AGES = [50., 55., 60., 65., 70., 75., 80., 85.]
FROM = 50.0

S_SLOW, S_FAST = 0.10, 0.13          # DNMT3A-like and TET2-like base fitness


def run(people, seed, s_ramp=None, wt_decline=0.0, s=(S_SLOW, S_FAST)):
    return simulate_cohort(N=N, s=list(s), s_weights=[0.5, 0.5], mu=MU,
                           years=max(AGES), people=people, dt=DT, max_clones=64,
                           seed=seed, record_ages=AGES, s_ramp=s_ramp,
                           wt_decline=wt_decline, age_effects_from=FROM)


def slope_by_class(r, s_value, ages=AGES, limit=LIMIT, rng=None, n_sample=None):
    """
    Slope of log(VAF) on age for one fitness class, one draw per person.

    This is the estimator stage E showed cannot recover s. It is used here for
    what it CAN do: compare two worlds measured the same way. A biased ruler
    still tells you which of two objects is longer.
    """
    people = r["sizes"][0].shape[0]
    rng = rng or np.random.default_rng(0)
    groups = np.array_split(rng.permutation(people), len(ages))
    xs, ys = [], []
    for k, (age, idx) in enumerate(zip(ages, groups)):
        v = sizes_to_vaf(r["sizes"][k][idx], N)
        f = r["fitness_at"][k][idx]
        m = (v >= limit) & np.isclose(f, s_value)
        xs.append(np.full(m.sum(), age)); ys.append(v[m])
    x = np.concatenate(xs); y = np.concatenate(ys)
    if n_sample is not None and x.size > n_sample:
        pick = rng.choice(x.size, size=n_sample, replace=False)
        x, y = x[pick], y[pick]
    if x.size < 12 or np.ptp(x) == 0:
        return np.nan, x.size
    return float(np.polyfit(x, np.log(y), 1)[0]), int(x.size)


# ================================================== 1. the two mechanisms
section("1. TWO MECHANISMS, AND WHY THEY ARE NOT INTERCHANGEABLE")

PEOPLE = 6_000
worlds = {
    "no age effect":                 dict(),
    "TET2 fitness ramps 3%/yr":      dict(s_ramp=[0.0, 0.03]),
    "wild type decays 1%/yr":        dict(wt_decline=0.01),
}
print(f"{'world':<28} {'DNMT3A-like':>26} {'TET2-like':>26}")
print(f"{'':<28} {'n':>7} {'med VAF':>9} {'slope':>9} {'n':>7} {'med VAF':>9} {'slope':>9}")
sig = {}
for name, kw in worlds.items():
    r = run(PEOPLE, 31, **kw)
    row = f"{name:<28}"
    rec = {}
    for s in (S_SLOW, S_FAST):
        k = AGES.index(70.0)
        v = sizes_to_vaf(r["sizes"][k], N); f = r["fitness_at"][k]
        m = (v >= LIMIT) & np.isclose(f, s)
        sl, n = slope_by_class(r, s, rng=np.random.default_rng(5))
        rec[s] = (int(m.sum()), float(np.median(v[m])), sl)
        row += f"{m.sum():>7} {np.median(v[m]):>9.3f} {sl:>9.4f}"
    sig[name] = rec
    print(row)

base = sig["no age effect"]
ramp = sig["TET2 fitness ramps 3%/yr"]
decl = sig["wild type decays 1%/yr"]
print(f"""
THE SIGNATURE, stated as a ratio to the no-effect world at age 70:

                          DNMT3A-like median VAF   TET2-like median VAF
  TET2 fitness ramps            x{ramp[S_SLOW][1]/base[S_SLOW][1]:.2f}                   x{ramp[S_FAST][1]/base[S_FAST][1]:.2f}
  wild type decays              x{decl[S_SLOW][1]/base[S_SLOW][1]:.2f}                   x{decl[S_FAST][1]/base[S_FAST][1]:.2f}

A driver-specific ramp leaves the other gene alone. A decaying wild type drags
it along — here it inflates the untouched DNMT3A-like class by a factor of
{decl[S_SLOW][1]/base[S_SLOW][1]:.1f}. That is the discriminator, and it needs no new sequencing: it asks
whether the age effect is specific to TET2 or shared with every driver.""")


# ============================================ 2. what the real data shows
section("2. WHAT THE OBSERVED DATA SHOWS, AND HOW FAR IT GOES")

raw = pd.read_csv(OBS_CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
for c in ("age", "VAF"):
    raw[c] = pd.to_numeric(raw[c], errors="coerce")
obs = raw.dropna(subset=["age", "VAF"])
obs = obs[obs.study.isin(ELIGIBLE) & (obs.VAF >= LIMIT)]

sub("2.1 The two genes, measured the same way")
def obs_slope(frame):
    return float(np.polyfit(frame.age, np.log(frame.VAF), 1)[0])

genes = {}
for g in ("TET2", "DNMT3A"):
    x = obs[obs.gene == g]
    genes[g] = x
    print(f"  {g:<8} n={len(x):>4}  cohorts={x.study.nunique()}  "
          f"median age {x.age.median():.0f}  slope {obs_slope(x):+.4f}")

gap = obs_slope(genes['TET2']) - obs_slope(genes['DNMT3A'])
print(f"\n  difference in slope, TET2 - DNMT3A: {gap:+.4f}")

sub("2.2 Does that gap survive the sample size?")
rng = np.random.default_rng(11)
B = 4000
boot = np.empty(B)
for b in range(B):
    a = genes["TET2"].sample(len(genes["TET2"]), replace=True, random_state=int(rng.integers(1e9)))
    c = genes["DNMT3A"].sample(len(genes["DNMT3A"]), replace=True, random_state=int(rng.integers(1e9)))
    boot[b] = obs_slope(a) - obs_slope(c)
lo, hi = np.quantile(boot, [0.025, 0.975])
p_gt = float((boot > 0).mean())
print(f"  bootstrap 95% interval: {lo:+.4f} to {hi:+.4f}")
print(f"  P(TET2 steeper than DNMT3A) = {p_gt:.2f}")
print(f"""
  {'SEPARATED' if p_gt > 0.95 else 'NOT SEPARATED'} at the 95% level.

  TET2's age slope is about twice DNMT3A's in this sample, which is the
  direction Fabre report. With {len(genes['TET2'])} TET2 variants from {genes['TET2'].study.nunique()} cohorts the interval
  spans zero, so the data is consistent with a TET2-specific effect and equally
  consistent with none.""")


# ====================================== 3. could this dataset ever tell?
section("3. COULD A DATASET THIS SIZE DETECT THE RAMP AT ALL?")

print(f"""Rather than argue about it, measure it. Build two worlds that are deliberately
hard to tell apart: one with constant fitness, one where TET2 fitness ramps,
with the ramping world's base fitness lowered so that BOTH produce the same
median detected VAF at age 70. An estimator that only sees age 70 cannot
separate them by construction; the question is whether the age profile can.

Then draw {len(genes['TET2'])} variants — the real TET2 sample size — from each world, many
times, and ask how often the slope from the ramping world exceeds the slope
from the constant one.""")

sub("3.1 Calibrating the two worlds to agree at age 70")
RAMP = 0.03
def median_at_70(s_fast, ramp):
    r = run(3_000, 44, s_ramp=[0.0, ramp], s=(S_SLOW, s_fast))
    k = AGES.index(70.0)
    v = sizes_to_vaf(r["sizes"][k], N); f = r["fitness_at"][k]
    m = (v >= LIMIT) & np.isclose(f, s_fast)
    return float(np.median(v[m])), int(m.sum())

target, n_t = median_at_70(S_FAST, 0.0)
print(f"  constant world: s = {S_FAST}, median VAF at 70 = {target:.4f}  (n={n_t})")
best = None
for cand in np.round(np.arange(0.08, 0.135, 0.005), 3):
    med, n = median_at_70(cand, RAMP)
    d = abs(med - target)
    if best is None or d < best[2]:
        best = (cand, med, d, n)
print(f"  ramping world:  s = {best[0]}, ramp {RAMP:.0%}/yr, "
      f"median VAF at 70 = {best[1]:.4f}  (n={best[3]})")
print(f"  the two agree at age 70 to within {abs(best[1]-target)/target:.1%}")

sub("3.2 Power at the real sample size")
S_RAMPED = float(best[0])
REPS = 60
N_TET2 = len(genes["TET2"])

# Simulate each replicate ONCE and reuse it for every sample size. The first
# version re-ran the whole pair of cohorts for each n, which quadrupled the
# work to produce identical worlds.
print(f"  simulating {REPS} replicate pairs", end="", flush=True)
pairs = []
for i in range(REPS):
    pairs.append((run(2_500, 200 + i),
                  run(2_500, 900 + i, s_ramp=[0.0, RAMP], s=(S_SLOW, S_RAMPED))))
    if (i + 1) % 10 == 0:
        print(f" {100*(i+1)//REPS}%", end="", flush=True)
print()

def power_at(n):
    a, b = [], []
    for i, (rc, rr) in enumerate(pairs):
        x, _ = slope_by_class(rc, S_FAST, rng=np.random.default_rng(1000 + i), n_sample=n)
        y, _ = slope_by_class(rr, S_RAMPED, rng=np.random.default_rng(1000 + i), n_sample=n)
        if np.isfinite(x): a.append(x)
        if np.isfinite(y): b.append(y)
    a, b = np.array(a), np.array(b)
    return float(np.mean(b[:, None] > a[None, :])), a, b

auc, sl_const, sl_ramp = power_at(N_TET2)
print(f"  constant world, slope at n={N_TET2}:  {sl_const.mean():+.4f} +/- {sl_const.std():.4f}")
print(f"  ramping  world, slope at n={N_TET2}:  {sl_ramp.mean():+.4f} +/- {sl_ramp.std():.4f}")
print(f"  P(ramping sample looks steeper than constant sample) = {auc:.2f}")
print("""
  0.50 would mean the two worlds are indistinguishable at this sample size;
  1.00 would mean always told apart. Both were built to agree exactly at age
  70, so this measures only what the AGE PROFILE adds on top.""")

sub("3.3 How much data would be needed")
print(f"{'n TET2 variants':>16} {'P(correct call)':>17}")
POWER = {}
for n in (N_TET2, 100, 300, 1000, 3000):
    POWER[n] = power_at(n)[0]
    print(f"{n:>16} {POWER[n]:>17.2f}", flush=True)


# ================================ 4. the comparison that removes composition
section("4. TET2 AGAINST DNMT3A, WITHIN EACH COHORT")

print(f"""Section 2 compared the two genes across a pooled set and found the gap does not
survive resampling. Section 3 says why: at n={N_TET2} even a real 3%/yr ramp is called
correctly only {POWER[N_TET2]:.0%} of the time. The floor at VAF {LIMIT} is what costs the sample —
it cuts TET2 from 75 variants in 4 cohorts to {N_TET2} in 2.

Dropping the floor recovers the sample and reintroduces the confounder stage E
identified: cohorts differ in sensitivity AND in age distribution, so a pooled
slope partly measures which studies survived the filter.

Unless the comparison is made WITHIN each cohort. A cohort has one protocol,
one detection limit and one age distribution, so TET2 and DNMT3A inside it face
identical measurement. The difference of slopes is then free of composition by
construction, and the cohorts can be combined afterwards.""")

sub("4.1 Per cohort")
pool = raw.dropna(subset=["age", "VAF"])
pool = pool[pool.study.isin(ELIGIBLE)]

rows = []
print(f"{'cohort':<14} {'n TET2':>7} {'n DNMT3A':>9} {'slope TET2':>11} "
      f"{'slope DNMT3A':>13} {'difference':>11}")
for study, g in pool.groupby("study"):
    tt, dd = g[g.gene == "TET2"], g[g.gene == "DNMT3A"]
    if len(tt) < 8 or len(dd) < 8 or tt.age.nunique() < 4:
        continue
    a, b = obs_slope(tt), obs_slope(dd)
    rows.append((study, len(tt), len(dd), a, b, a - b))
    print(f"{study:<14} {len(tt):>7} {len(dd):>9} {a:>11.4f} {b:>13.4f} {a-b:>11.4f}")

per = pd.DataFrame(rows, columns=["study", "n_t", "n_d", "s_t", "s_d", "diff"])
print(f"\n  cohorts usable: {len(per)}   TET2 variants: {per.n_t.sum()}   "
      f"DNMT3A variants: {per.n_d.sum()}")
print(f"  every cohort positive: {bool((per['diff'] > 0).all())}")

sub("4.2 Combined, by resampling within cohorts")
rng = np.random.default_rng(7)
B = 4000
boot = np.empty(B)
for b in range(B):
    diffs, weights = [], []
    for _, row in per.iterrows():
        g = pool[pool.study == row.study]
        tt = g[g.gene == "TET2"].sample(int(row.n_t), replace=True,
                                        random_state=int(rng.integers(1e9)))
        dd = g[g.gene == "DNMT3A"].sample(int(row.n_d), replace=True,
                                          random_state=int(rng.integers(1e9)))
        if tt.age.nunique() < 3 or dd.age.nunique() < 3:
            continue
        diffs.append(obs_slope(tt) - obs_slope(dd))
        weights.append(row.n_t)
    boot[b] = np.average(diffs, weights=weights) if diffs else np.nan

boot = boot[np.isfinite(boot)]
lo, hi = np.quantile(boot, [0.025, 0.975])
p_gt = float((boot > 0).mean())
point = float(np.average(per["diff"], weights=per.n_t))
print(f"  weighted mean difference: {point:+.4f} per year")
print(f"  bootstrap 95% interval:   {lo:+.4f} to {hi:+.4f}")
print(f"  P(TET2 steeper than DNMT3A) = {p_gt:.3f}")
print(f"\n  {'SEPARATED' if p_gt > 0.95 else 'NOT SEPARATED'} at the 95% level.")

sub("4.3 What this settles, and what it does not")
print(f"""{'The TET2 age effect is present in this data once composition is removed.' if p_gt > 0.95
   else 'Even within cohorts the effect does not reach the 95% level here.'}

WHAT IT CANNOT SAY. Section 1 showed the two mechanisms differ in whether the
OTHER gene comes along, and this test measures exactly that difference — so a
positive result rules out a purely wild-type explanation, in which both genes
would rise together and the difference would be zero. It does NOT distinguish
"TET2 gets intrinsically better" from "TET2 is protected while the wild type
degrades". In a Moran process those are the same model: fitness is relative,
and "everyone else got worse except me" is arithmetically identical to "I got
better". The mouse work and Fabre are not in conflict — they are describing the
same relative change from opposite sides.

WHAT WOULD DISTINGUISH THEM is something this model cannot see, because it has
no absolute clock: the ABSOLUTE output of wild-type haematopoiesis with age.
That is a measurement, not an inference — and it is why the mouse
transplantation experiments were done.""")


# ================================================================== figure
section("FIGURE")

fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))

ax = axes[0]
names = list(worlds)
x = np.arange(len(names))
w = 0.36
for k, (s, col, lab) in enumerate([(S_SLOW, S1, "DNMT3A-like"),
                                   (S_FAST, S2, "TET2-like")]):
    vals = [sig[n][s][1] / sig["no age effect"][s][1] for n in names]
    ax.bar(x + (k - 0.5) * w, vals, width=w - 0.03, color=col, label=lab, zorder=2)
ax.axhline(1.0, color=INK2, linewidth=2, zorder=3)
ax.set_xticks(x)
ax.set_xticklabels(["no age\neffect", "TET2 fitness\nramps", "wild type\ndecays"], fontsize=9)
ax.set_ylabel("median VAF at 70, ratio to no effect")
ax.set_title("Only one mechanism spares the other gene")
ax.legend(loc="upper left")
ax.grid(axis="x", visible=False)
clean(ax)

ax = axes[1]
y = np.arange(len(per))[::-1]
ax.axvline(0, color=INK2, linewidth=2, zorder=3)
ax.scatter(per["diff"], y, s=per.n_t * 4 + 40, color=S3, zorder=4,
           edgecolor=SURF, linewidth=2)
ax.set_yticks(y)
ax.set_yticklabels([f"{r.study}\n(n={int(r.n_t)})" for _, r in per.iterrows()], fontsize=8)
ax.set_xlabel("slope TET2 − slope DNMT3A  (per year)")
ax.set_title("Within each cohort, same protocol")
ax.grid(axis="y", visible=False)
clean(ax)

ax = axes[2]
ns = sorted(POWER)
ax.plot(ns, [POWER[n] for n in ns], "-o", color=S1, linewidth=2.5, markersize=8,
        markeredgecolor=SURF, markeredgewidth=2, zorder=3)
ax.axhline(0.95, color=S2, linewidth=2, zorder=2)
ax.text(ns[-1], 0.955, "95%", color=S2, fontsize=9, ha="right", va="bottom")
ax.axhline(0.5, color=BASE, linewidth=2, zorder=1)
ax.set_xscale("log")
ax.set_xticks(ns); ax.set_xticklabels([str(n) for n in ns], fontsize=9)
ax.set_ylim(0.4, 1.05)
ax.set_xlabel("TET2 variants available")
ax.set_ylabel("P(a real 3%/yr ramp is called correctly)")
ax.set_title("What the sample size allows")
clean(ax)

fig.tight_layout()
fig.savefig(FIGS / "13_tet2_age.png"); plt.close(fig)
print("  13_tet2_age.png")
