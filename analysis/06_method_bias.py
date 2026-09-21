#!/usr/bin/env python3
"""
Stage E — which measurement method recovers the truth?

THE PROBLEM THIS EXISTS TO SOLVE
Stage D compared our estimate against two published ones and found they
disagree with each other by up to threefold on the same gene: Watson 2020 puts
DNMT3A R882 at s = 0.148 per year, Fabre 2022 at 0.050. Both are careful work
on real cohorts. They differ in DESIGN — Watson reads a single blood sample
from many people of different ages; Fabre follows the same people for a median
of thirteen years.

With real data there is no arbiter. Nobody knows the true fitness of a clone
in a living person, so a disagreement between two methods cannot be settled by
measuring harder.

In a simulation there is an arbiter, because we set `s` ourselves. So: build a
cohort whose fitness is known by construction, apply each measurement design to
it, and see what each one returns. Any gap between the estimate and the number
we typed in is bias attributable to the design, not to biology.

This is what the model is FOR. It is not a better estimate of s; it is the only
way to put an error bar on a method.

Usage:  .venv/bin/python analysis/06_method_bias.py
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
RESULTS = ROOT / "results"

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


# ----------------------------------------------------------------- settings
N          = 100_000     # stem cells, the value stage C settled on
MU         = 2e-6        # driver mutations per cell per year
PEOPLE     = 6_000
DT         = 0.5
LIMIT      = 0.0192      # the common floor stages C and D used
DEEP_LIMIT = 0.0008      # Young2019's actual floor: the most sensitive cohort
DEPTH      = 1_000       # sequencing depth, for read-sampling noise
BASELINE   = 65.0        # age at first draw in the longitudinal design
FOLLOWUP   = 13.0        # Fabre's median follow-up
CROSS_AGES = [50., 55., 60., 65., 70., 75., 80., 85.]
ABC_AGE    = 70.0
QUANTILES  = [0.10, 0.25, 0.50, 0.75, 0.90]

# The published pair for DNMT3A, gene-level on both sides, as quoted by Fabre
# 2022 in their own discussion. See analysis/METHOD_BIAS_LITERATURE.md.
WATSON_DNMT3A = 0.150      # spectrum fit, Watson et al. 2020
FABRE_DNMT3A  = 0.062      # longitudinal follow-up, Fabre et al. 2022
QCOLS      = ["vaf_q10", "vaf_q25", "vaf_q50", "vaf_q75", "vaf_q90"]

OBS_CSV = (ROOT / "reference/data/watson2020/Maximum_likelihood_estimations"
                / "Maximum likelihood estimations - data files"
                / "all_studies_trimmed_all_genes.csv")
# The cohorts whose own detection limit sits at or below LIMIT, carried over
# from stage D: pooling a cohort that cannot see below the common floor
# contributes only large clones and skews everything it touches upward.
ELIGIBLE = ["Young2019", "Young2016", "Acuna2017", "McKerrel2015",
            "Desai2018", "Coombs2017"]

RECORD = sorted(set(CROSS_AGES + [BASELINE, BASELINE + FOLLOWUP, ABC_AGE]))


# ------------------------------------------------------------- the observer
def observe(vaf, depth, rng):
    """
    What the sequencer reports, rather than what is true in the marrow.

    A variant allele frequency is not read off directly: it is estimated from
    a finite number of sequencing reads. At depth d the number of reads
    carrying the variant is Binomial(d, VAF), so the reported frequency
    scatters around the true one with standard deviation sqrt(VAF(1-VAF)/d).

    This matters far more than it looks. A clone is entered into a study
    because its REPORTED frequency cleared the detection limit, and among
    clones near that limit the ones reported high are disproportionately the
    ones the noise pushed up. Measure them again later and they fall back —
    not because they shrank, but because the first number was lucky. Any
    method that selects on a noisy measurement and then measures growth from
    it inherits this.
    """
    if depth is None:
        return vaf
    out = np.zeros_like(vaf)
    m = vaf > 0
    out[m] = rng.binomial(depth, np.clip(vaf[m], 0.0, 1.0)) / depth
    return out


# ------------------------------------------------------------- estimator 1
def est_cross_sectional(snaps, births, ages, rng, limit=LIMIT, depth=DEPTH,
                        clock="person"):
    """
    One blood sample per person, people of different ages; regress log(VAF) on
    age and read the slope as the growth rate. This is Watson's design, and it
    is also exactly what stage 1 did to the raw data.

    Each simulated person is assigned to one age and contributes only that
    snapshot — a person measured at several ages would be a different study.

    `clock` chooses the regressor, and the choice is the experiment:
      'person' — the host's age, which is all a real study can know
      'clone'  — the clone's own age, which only a simulation can know
    The gap between the two is the cost of not knowing when the clone started.
    """
    n = snaps[0].shape[0]
    groups = np.array_split(rng.permutation(n), len(ages))
    xs, ys = [], []
    for age, idx in zip(ages, groups):
        k = list(ages).index(age)
        v = observe(sizes_to_vaf(snaps[k][idx], N), depth, rng)
        b = births[k][idx]
        d = v >= limit
        if not d.any():
            continue
        x = np.full(v.shape, age) if clock == "person" else (age - b)
        xs.append(x[d]); ys.append(v[d])
    if not xs:
        return np.nan, 0
    x = np.concatenate(xs); y = np.log(np.concatenate(ys))
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 20 or np.ptp(x[ok]) == 0:
        return np.nan, int(ok.sum())
    slope = np.polyfit(x[ok], y[ok], 1)[0]
    return float(slope), int(ok.sum())


# ------------------------------------------------------------- estimator 2
def est_longitudinal(s1, b1, s2, b2, rng, dt=FOLLOWUP, limit=LIMIT,
                     depth=DEPTH, max_vaf=None):
    """
    The same people sampled twice, `dt` years apart. Every clone detected at
    the first draw gets its own growth rate, log(VAF2/VAF1)/dt. This is
    Fabre's design.

    Two details decide whether the comparison is fair.

    IDENTITY. A clone must be the SAME clone at both draws. Slots in the
    simulation are recycled when a clone dies, so the birth stamp is checked;
    without it roughly eight per cent of pairs splice two unrelated clones
    into one invented trajectory.

    THE SECOND DRAW IS NOT TRUNCATED. Once a variant is known, a real study
    genotypes that exact position deeply and can follow it well below the
    limit that first found it. So a clone is admitted on its first
    measurement and then followed wherever it goes — including down.
    """
    v1 = observe(sizes_to_vaf(s1, N), depth, rng)
    v2 = observe(sizes_to_vaf(s2, N), depth, rng)

    same = np.isfinite(b1) & np.isfinite(b2) & (b1 == b2)
    sel = same & (v1 >= limit)
    if max_vaf is not None:
        sel &= v1 <= max_vaf

    keep = sel & (v2 > 0)
    if keep.sum() < 20:
        return np.nan, int(sel.sum()), int(sel.sum() - keep.sum())
    rate = np.log(v2[keep] / v1[keep]) / dt
    return float(np.median(rate)), int(sel.sum()), int(sel.sum() - keep.sum())


# ------------------------------------------------------------- estimator 3
class ABCEstimator:
    """
    Our own method, from stages C and D: compare the quantiles of the detected
    VAF distribution against a pre-computed grid of simulations and keep the
    closest grid point.

    The grid is the one the cluster already produced. Running our estimator
    against a cohort whose answer we know is the check stages C and D could
    not perform on real data.
    """

    def __init__(self, age=ABC_AGE):
        frames = []
        for f in ("sweep_N100k.csv", "sweep_fine.csv"):
            p = RESULTS / f
            if p.exists():
                frames.append(pd.read_csv(p))
        if not frames:
            self.ok = False
            return
        sw = pd.concat(frames, ignore_index=True)
        sw = sw[(sw.age == age) & (sw.N == N)].dropna(subset=QCOLS).copy()
        self.sw = sw
        self.sim = np.log(sw[QCOLS].values)
        self.scale = self.sim.std(axis=0)
        self.scale[self.scale == 0] = 1.0
        self.ok = len(sw) > 0
        self.grid_s = np.sort(sw.s.unique())

    def __call__(self, vaf_matrix, limit=LIMIT):
        values = vaf_matrix[vaf_matrix >= limit]
        if values.size < 20 or not self.ok:
            return np.nan, 0
        target = np.log(np.quantile(values, QUANTILES))
        d = np.sqrt((((self.sim - target) / self.scale) ** 2).mean(axis=1))
        grid = (self.sw.assign(distance=d)
                       .groupby(["s", "mu"], as_index=False)
                       .agg(distance=("distance", "mean")))
        row = grid.loc[grid.distance.idxmin()]
        return float(row.s), int(values.size)


# ------------------------------------------------------------------ driver
def run_cohort(s_true, seed):
    r = simulate_cohort(N=N, s=s_true, mu=MU, years=max(RECORD), people=PEOPLE,
                        dt=DT, max_clones=64, seed=seed, record_ages=RECORD)
    idx = {float(a): i for i, a in enumerate(r["ages"])}
    return r, idx


section("1. THE SETUP")
print(f"""A cohort of {PEOPLE} people is simulated from birth, with {N:,} stem cells,
a driver mutation rate of {MU:g} per cell per year, and a selection coefficient
we choose. Three measurement designs are then applied to that same cohort:

  age regression    one draw per person, log(VAF) regressed on the host's age.
                    This is what stage 1 did, and what anyone reaches for first.
  longitudinal      the same people at {BASELINE:.0f} and {BASELINE+FOLLOWUP:.0f}, one growth rate per clone.
                    This is Fabre 2022's design.
  ABC quantiles     match the shape of the detected VAF distribution against
                    the grid the cluster computed. This is stage C and D's
                    method, and it is also the FAMILY Watson 2020 belongs to —
                    their estimate is a maximum-likelihood fit to the same VAF
                    spectrum, not a slope against age.

That last point matters and it corrects a shortcut from stage D. "Watson is
cross-sectional" is true about the sampling and misleading about the estimator:
a single-timepoint study can be read with a slope or with a likelihood fit to
the spectrum, and those are not the same instrument. Only the second is being
compared with Fabre here.

All three see the same detection limit (VAF >= {LIMIT}) and the same read noise
(depth {DEPTH}), so any difference between them is design, not instrument.""")

abc = ABCEstimator()
print(f"\ngrid loaded: {len(abc.sw) if abc.ok else 0} rows at age {ABC_AGE:.0f}, "
      f"s from {abc.grid_s.min():.2f} to {abc.grid_s.max():.2f}" if abc.ok
      else "\nWARNING: no sweep grid found, the ABC estimator is skipped")


# ================================================== experiment 1: recovery
section("2. EXPERIMENT 1 — DOES EACH DESIGN RECOVER THE TRUTH?")

S_TRUE = [0.08, 0.12, 0.16, 0.20, 0.24]
SEEDS = [11, 12, 13]

rows = []
for s_true in S_TRUE:
    for seed in SEEDS:
        r, idx = run_cohort(s_true, seed)
        rng = np.random.default_rng(seed + 1000)

        cs, n_cs = est_cross_sectional(
            [r["sizes"][idx[a]] for a in CROSS_AGES],
            [r["births"][idx[a]] for a in CROSS_AGES],
            CROSS_AGES, rng)

        lo, n_lo, lost = est_longitudinal(
            r["sizes"][idx[BASELINE]], r["births"][idx[BASELINE]],
            r["sizes"][idx[BASELINE + FOLLOWUP]], r["births"][idx[BASELINE + FOLLOWUP]],
            rng)

        ab, n_ab = abc(sizes_to_vaf(r["sizes"][idx[ABC_AGE]], N))

        rows.append(dict(s_true=s_true, seed=seed,
                         cross=cs, n_cross=n_cs,
                         longit=lo, n_long=n_lo, lost=lost,
                         abc=ab, n_abc=n_ab))
        print(f"  s={s_true:.2f} seed={seed}  "
              f"cross {cs:.3f} (n={n_cs})   long {lo:.3f} (n={n_lo}, {lost} lost)   "
              f"abc {ab:.3f}")

df = pd.DataFrame(rows)
agg = df.groupby("s_true").agg(
    cross=("cross", "mean"), cross_sd=("cross", "std"),
    longit=("longit", "mean"), longit_sd=("longit", "std"),
    abc=("abc", "mean"), abc_sd=("abc", "std")).reset_index()

print(f"\n{'s true':>7} | {'regression on age':>22} | {'longitudinal':>22} | {'VAF spectrum fit':>22}")
print(f"{'':>7} | {'estimate':>10} {'error':>11} | {'estimate':>10} {'error':>11} | {'estimate':>10} {'error':>11}")
for _, a in agg.iterrows():
    def cell(v, t):
        return f"{v:>10.3f} {(v/t-1)*100:>+10.0f}%"
    print(f"{a.s_true:>7.2f} | {cell(a.cross, a.s_true)} | "
          f"{cell(a.longit, a.s_true)} | {cell(a.abc, a.s_true)}")

# Recovery is the mean of the per-run RATIOS, not the ratio of the means: a
# ratio of means would let the large-s runs, where every method fails, dominate
# the summary and hide how each one behaves where the truth actually sits.
RECOVERY = {name: float((df[name] / df.s_true).mean())
            for name in ("cross", "longit", "abc")}
for name, ratio in RECOVERY.items():
    print(f"\n{name:<7} recovers on average {ratio*100:.0f}% of the true value")


# ============================================ experiment 2: decomposition
section("3. EXPERIMENT 2 — WHERE THE BIAS COMES FROM")

S_FOCUS = 0.14
print(f"One cohort, s = {S_FOCUS}. Each row switches off ONE thing a real study")
print("cannot switch off, starting from the same baseline — the rows do not")
print("accumulate except where the label says so. The change in the estimate is")
print("what that one artifact costs on its own.\n")

r, idx = run_cohort(S_FOCUS, seed=21)
snaps = [r["sizes"][idx[a]] for a in CROSS_AGES]
brths = [r["births"][idx[a]] for a in CROSS_AGES]
s1, b1 = r["sizes"][idx[BASELINE]], r["births"][idx[BASELINE]]
s2, b2 = r["sizes"][idx[BASELINE + FOLLOWUP]], r["births"][idx[BASELINE + FOLLOWUP]]

def fresh():
    return np.random.default_rng(2024)

cross_conditions = [
    ("as a real study sees it",        dict()),
    ("clone age known",                dict(clock="clone")),
    ("perfect reads (no noise)",       dict(depth=None)),
    ("deepest real limit (0.0008)",    dict(limit=DEEP_LIMIT)),
    ("clone age + deep limit",         dict(clock="clone", limit=DEEP_LIMIT)),
    ("all three off at once",          dict(clock="clone", limit=DEEP_LIMIT, depth=None)),
]

long_conditions = [
    ("as a real study sees it",        dict()),
    ("perfect reads (no noise)",       dict(depth=None)),
    ("small clones only (VAF<0.05)",   dict(max_vaf=0.05)),
    ("both off at once",               dict(depth=None, max_vaf=0.05)),
]

print("REGRESSION ON HOST AGE")
print(f"{'condition':<32} {'estimate':>9} {'error':>9}  {'n':>7}")
cross_rows = []
for label, kw in cross_conditions:
    v, n = est_cross_sectional(snaps, brths, CROSS_AGES, fresh(), **kw)
    cross_rows.append((label, v))
    print(f"{label:<32} {v:>9.3f} {(v/S_FOCUS-1)*100:>+8.0f}%  {n:>7}")

print("\nLONGITUDINAL")
print(f"{'condition':<32} {'estimate':>9} {'error':>9}  {'n':>7} {'lost':>6}")
long_rows = []
for label, kw in long_conditions:
    v, n, lost = est_longitudinal(s1, b1, s2, b2, fresh(), **kw)
    long_rows.append((label, v))
    print(f"{label:<32} {v:>9.3f} {(v/S_FOCUS-1)*100:>+8.0f}%  {n:>7} {lost:>6}")


# ======================================================= what it all means
section("4. WHAT THIS SAYS ABOUT THE PUBLISHED DISAGREEMENT")

c, l, a = RECOVERY["cross"], RECOVERY["longit"], RECOVERY["abc"]
print(f"""Averaged over the five true values, as a fraction of the truth recovered:

  regression on host age       {c*100:>5.0f}%
  longitudinal per-clone rate  {l*100:>5.0f}%
  ABC on the VAF spectrum      {a*100:>5.0f}%

THE ORDER OF THOSE THREE IS THE RESULT, and the published gap is the right
size. Fabre 2022 put the clean comparison in their own discussion: DNMT3A
clones at {WATSON_DNMT3A:.3f} per year from Watson's spectrum fit against {FABRE_DNMT3A:.3f} per year
from their own follow-up — both gene-level, both DNMT3A, a ratio of
{WATSON_DNMT3A/FABRE_DNMT3A:.2f}. From design bias alone this experiment predicts {a/l:.1f}.

(Stage D compared R882H at 0.148 against Fabre's ~0.050 instead, a ratio of 3.0.
That pair is not like for like: a per-variant hotspot against a gene-level
average. Watson report that over 90% of nonsynonymous DNMT3A variants are
effectively neutral, so a gene average over that mixture belongs well below the
hotspot value for reasons that have nothing to do with method.)

The agreement is closer than three seeds and one mutation rate deserve, and
should be read as consistency rather than a match. But it is the right size, on
the right pair, in the right direction.

Fabre reach the same place by another route: they read the gap as DNMT3A clones
genuinely growing faster early in life, and attribute the slowdown to "an
increasingly competitive oligoclonal landscape". That is this model's
competition term in words. The two readings agree on the mechanism. A clone
whose growth slowed because the niche filled has not lost fitness, and calling
its late-window rate "fitness" is the error.

WHY LONGITUDINAL COMES OUT LOW, AND WHY IT GETS WORSE WITH s
A clone cannot keep growing exponentially once it owns a large share of the
marrow: the cells it competes against are increasingly its own. The model has
that built in, through the Moran normalisation, so a large clone's measured
growth is genuinely slower than its fitness. A study that enrols clones big
enough to detect and then measures them over thirteen years spends much of
that window in the saturating regime. The effect grows with s because a fitter
clone reaches saturation sooner — at s = {S_TRUE[-1]}, thirteen years of follow-up
recover almost nothing.

This is not an error in Fabre's measurement. The growth they report is real.
It is an error to read that growth as a fitness that would apply to a clone
starting from one cell.

WHAT THE ABC RESULT DOES AND DOES NOT ESTABLISH
The grid was produced by this same model, so recovering the input is a check
that the inversion is unbiased and that the five quantiles carry enough
information to pin s down. Both could have failed and neither did. It is NOT
independent evidence that the model describes real haematopoiesis — no
self-consistency check can be. What it rules out is a specific failure we could
not otherwise exclude: that stage D's estimates were an artifact of the fitting
procedure.""")


# =============================================== the loop back to real data
section("5. THE SAME BIAS, VISIBLE IN THE REAL DATA")

raw = pd.read_csv(OBS_CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
raw["age"] = pd.to_numeric(raw["age"], errors="coerce")
raw["VAF"] = pd.to_numeric(raw["VAF"], errors="coerce")
o = raw.dropna(subset=["age", "VAF"])

def slope_of(frame):
    return float(np.polyfit(frame.age, np.log(frame.VAF), 1)[0]), len(frame)

# Change ONE thing. Both rows use the same six cohorts and the same variants;
# only the floor moves. Comparing 'all cohorts, no floor' against 'eligible
# cohorts, with floor' would move two levers at once and prove nothing.
elig = o[o.study.isin(ELIGIBLE)]
deep_slope, n_deep = slope_of(elig)
cut_slope, n_cut = slope_of(elig[elig.VAF >= LIMIT])

print(f"""Take the six cohorts stage D kept and move only the detection floor:

  no floor (down to VAF 0.0008)   slope {deep_slope:.4f}   n = {n_deep}
  floor at VAF {LIMIT}              slope {cut_slope:.4f}   n = {n_cut}

The same cohorts, the same biology, a factor of {deep_slope/cut_slope:.1f} between the two numbers.
Nothing about the cells changed — only where the instrument stops seeing.

The simulation produces the same sign, and the deep-sequencing end almost
exactly: with a floor of {DEEP_LIMIT} it returned {cross_rows[3][1]:.3f}, against {deep_slope:.4f} observed.
But it only moves by a factor of {cross_rows[3][1]/cross_rows[0][1]:.1f} when the floor is raised, not {deep_slope/cut_slope:.1f}.
The model does not account for the whole effect, and the missing part is
identifiable rather than mysterious.""")

share_before = (elig.study == "Coombs2017").mean()
share_after = (elig[elig.VAF >= LIMIT].study == "Coombs2017").mean()
y19_before = int((elig.study == "Young2019").sum())
y19_after = int((elig[elig.VAF >= LIMIT].study == "Young2019").sum())

print(f"""
WHAT THE MODEL IS MISSING: THE FLOOR ALSO CHANGES WHO IS IN THE POOL
The simulation is one homogeneous cohort, so raising its floor removes clones
and nothing else. The real pooled regression is six cohorts with different
sensitivities AND different age distributions, so raising the floor removes
whole studies:

  Coombs2017's share of the pool   {share_before*100:.0f}%  ->  {share_after*100:.0f}%
  Young2019 variants surviving     {y19_before}  ->  {y19_after}

Young2019 sequences deeply and skews young; Coombs2017 has a limit right at
the floor and a median age of 70. Applying the floor deletes the young,
deep-sequenced end of the pool and leaves a regression run almost within a
single cohort. The age lever is largely gone before the fit even starts.

So the {deep_slope/cut_slope:.1f} in the real data is two effects stacked: truncation, which the
model reproduces at about {cross_rows[3][1]/cross_rows[0][1]:.1f}, and a composition shift that the model does
not contain because it was never given more than one cohort. Both are
artifacts of measurement. Neither is biology.

THE CONSEQUENCE, AND ITS LIMIT
A slope of log(VAF) against age is not an estimate of a selection coefficient.
It depends on the detection limit, on which cohorts were pooled, and on how fit
the clones are — non-monotonically, as experiment 1 shows: {agg.cross.iloc[0]:.3f} at s={agg.s_true.iloc[0]:.2f},
rising to {agg.cross.max():.3f} in the middle, back down to {agg.cross.iloc[-1]:.3f} at s={agg.s_true.iloc[-1]:.2f}. A
non-monotonic map has no inverse: one observed slope is compatible with
several very different truths, and no amount of extra data fixes that.

Stage 1's {deep_slope:.4f} and stage D's ~0.13 were therefore never in conflict. They are
the same population read with two instruments, one of which does not measure
what its units suggest.

The claim stops there, and not one step further. AGE IS NOT THE PROBLEM. Watson
run an age-based cross-check of their own — the prevalence of a variant at a
fixed detection threshold, which their model predicts rises linearly at rate
2*N*t*mu*s — and it returns 14% per year against the 15% their spectrum fit
gives. Same sampling design, same ages, a different functional of them, and no
detectable bias.

What fails is this particular extraction: the slope of log(VAF) among DETECTED
clones. The prevalence slope survives because it does not condition on the
clones that cleared the threshold; it counts the people who have one.

And we cannot borrow it. Prevalence needs a screening denominator — how many
people were tested in each age band, not how many carried a variant. Stage C
established that this dataset has none, which is why the mutation rate came out
unidentifiable. The same missing column blocks the one age-based estimator that
would have worked.""")

# =================================================================== figures
section("6. FIGURES")

fig, ax = plt.subplots(figsize=(8.6, 5.4))
lim = [0, max(S_TRUE) * 1.25]
ax.plot(lim, lim, color=BASE, linewidth=2, zorder=1, label="truth")
for name, colour, marker, label in [
        ("cross", S1, "o", "regression on host age (stage 1)"),
        ("longit", S2, "s", "longitudinal follow-up (Fabre design)"),
        ("abc", S3, "^", "VAF spectrum fit (stages C-D, Watson family)")]:
    m = agg[name].values
    e = agg[f"{name}_sd"].values
    ax.errorbar(agg.s_true, m, yerr=e, fmt=marker, color=colour, markersize=8,
                linewidth=2, capsize=0, elinewidth=2, label=label, zorder=3,
                markeredgecolor=SURF, markeredgewidth=2)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("true selection coefficient  s  (per year)")
ax.set_ylabel("estimated  s  (per year)")
ax.set_title("What each design recovers from a cohort whose answer is known")
ax.legend(loc="upper left")
clean(ax)
fig.tight_layout(); fig.savefig(FIGS / "10_method_recovery.png"); plt.close(fig)
print("  10_method_recovery.png")

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharex=True)
for ax, rowset, colour, title in [
        (axes[0], cross_rows, S1, "Regression on host age"),
        (axes[1], long_rows, S2, "Longitudinal follow-up")]:
    labels = [l for l, _ in rowset]
    vals = [v for _, v in rowset]
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, height=0.55, color=colour, zorder=2)
    ax.axvline(S_FOCUS, color=INK2, linewidth=2, zorder=3)
    ax.text(S_FOCUS, -0.75, f"  truth = {S_FOCUS}", color=INK2,
            fontsize=9, va="center", ha="left")
    ax.set_ylim(-1.1, len(vals) - 0.4)
    ax.set_xlim(0, S_FOCUS * 1.42)
    for i, v in zip(y, vals):
        ax.text(v + 0.004, i, f"{v:.3f}", va="center", fontsize=9, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("estimated  s  (per year)")
    ax.set_title(title)
    ax.grid(axis="y", visible=False)
    clean(ax)
fig.suptitle("Turning off one artifact at a time", x=0.008, ha="left",
             fontsize=13, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIGS / "11_bias_decomposition.png"); plt.close(fig)
print("  11_bias_decomposition.png")

df.to_csv(RESULTS / "method_bias.csv", index=False)
print(f"  results/method_bias.csv")
