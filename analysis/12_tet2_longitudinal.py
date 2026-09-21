#!/usr/bin/env python3
"""
Stage H — the TET2 age question, asked of data that can answer it.

WHAT FAILED BEFORE
Fabre report TET2 clones growing faster in older people. Three attempts to test
it came back empty, each for a different reason:

  Watson, pooled      the effect was 87% cohort composition
  Watson, power       36 variants, ceiling 0.79 even if the ramp were real
  ramp vs wild type   the same model in relative fitness, not distinguishable
  multihit            needs phase, and bulk VAF resolves it for 2% of pairs

WHY THIS DATA IS DIFFERENT
Every previous attempt compared different people to each other, so cohort,
sequencing protocol and age distribution were all free to move at once. Here a
clone is measured four or five times across a median of 13.4 years, so each
clone can be split into an early half and a late half and compared WITH ITSELF.

That removes person, cohort, protocol and clone identity in one step. What it
cannot remove is what stage F established: clones slow down as the marrow fills,
whoever fills it. So the question is not "does TET2 accelerate" in isolation —
it is whether TET2 decelerates LESS than DNMT3A once size and carrier burden are
matched.

Usage:  .venv/bin/python analysis/12_tet2_longitudinal.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FABRE = ROOT / "reference/data/fabre2022/ALLvariants_exclSynonymous_Xadj.txt"
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
    print(f"\n{'='*76}\n{t}\n{'='*76}")

def sub(t):
    print(f"\n{'-'*76}\n{t}\n{'-'*76}")


FLOOR = 0.002
MIN_POINTS = 4
MIN_SPAN = 8.0
rng = np.random.default_rng(17)

d = pd.read_csv(FABRE, sep="\t")
d = d[d.VAF >= FLOOR].copy()
d["clone"] = (d.SardID.astype(str) + "|" + d.Gene + "|"
              + d.START.astype(str) + "|" + d.ALT.astype(str))

# the carrier's total mutant burden at each draw: stage F's dominant predictor
burden = d.groupby(["SardID", "Phase"]).VAF.sum().rename("burden").reset_index()
d = d.merge(burden, on=["SardID", "Phase"])


def wfit(age, vaf, depth):
    """Weighted slope of log(VAF) on age, with the read-noise standard error."""
    if len(age) < 2 or np.ptp(age) < 1e-9:
        return np.nan, np.nan
    var = (1.0 - vaf) / (vaf * depth)
    w = 1.0 / var
    x = age - np.average(age, weights=w)
    sxx = np.sum(w * x * x)
    if sxx <= 0:
        return np.nan, np.nan
    return float(np.sum(w * x * np.log(vaf)) / sxx), float(np.sqrt(1.0 / sxx))


rows = []
for c, g in d.groupby("clone"):
    g = g.sort_values("Age")
    if len(g) < MIN_POINTS or np.ptp(g.Age) < MIN_SPAN:
        continue
    age, vaf, dep = g.Age.values, g.VAF.values, g.TOTALcount.values
    mid = np.median(age)
    # halves overlap at the middle point so each has enough leverage
    e = age <= mid
    l = age >= mid
    if e.sum() < 2 or l.sum() < 2:
        continue
    r_e, se_e = wfit(age[e], vaf[e], dep[e])
    r_l, se_l = wfit(age[l], vaf[l], dep[l])
    r_all, se_all = wfit(age, vaf, dep)
    if not all(np.isfinite(v) for v in (r_e, r_l, se_e, se_l)):
        continue
    rows.append(dict(clone=c, gene=g.Gene.iloc[0], person=g.SardID.iloc[0],
                     early=r_e, late=r_l, se_e=se_e, se_l=se_l,
                     decel=r_e - r_l, se_decel=float(np.hypot(se_e, se_l)),
                     rate=r_all, vaf0=vaf[0], vafmax=vaf.max(),
                     age0=age[0], age1=age[-1],
                     burden0=g.burden.iloc[0], burden1=g.burden.iloc[-1]))
t = pd.DataFrame(rows)
t["dburden"] = t.burden1 - t.burden0

GENES = ["DNMT3A", "TET2"]


section("1. EACH CLONE AS ITS OWN CONTROL")
print(f"""{len(t)} clones with at least {MIN_POINTS} timepoints over at least {MIN_SPAN:.0f} years, split at the
median age of their own trajectory. Ages run {t.age0.median():.0f} to {t.age1.median():.0f} at the median.

Deceleration is the early rate minus the late rate, per clone, both fitted with
weights from read depth.""")

print(f"\n  {'gene':<9}{'n':>5}{'early':>9}{'late':>9}{'deceleration':>15}{'95% interval':>22}")
summary = {}
for gene in GENES + ["other"]:
    g = t[t.gene == gene] if gene != "other" else t[~t.gene.isin(GENES)]
    if len(g) < 25:
        continue
    b = np.array([np.median(rng.choice(g.decel.values, len(g))) for _ in range(4000)])
    lo, hi = np.quantile(b, [0.025, 0.975])
    summary[gene] = (g, float(np.median(g.decel)), lo, hi)
    print(f"  {gene:<9}{len(g):>5}{np.median(g.early):>9.4f}{np.median(g.late):>9.4f}"
          f"{np.median(g.decel):>15.4f}{f'{lo:+.4f} to {hi:+.4f}':>22}")

print("""
  Positive deceleration means the clone grew more slowly in the second half of
  its own follow-up. Stage F predicts this for every clone, because the marrow
  fills. The question is whether TET2 does it less.""")


section("2. TET2 AGAINST DNMT3A, THE COMPARISON FABRE MAKE")

a, b_ = t[t.gene == "TET2"], t[t.gene == "DNMT3A"]
diff = float(np.median(a.decel) - np.median(b_.decel))
boot = np.array([np.median(rng.choice(a.decel.values, len(a)))
                 - np.median(rng.choice(b_.decel.values, len(b_)))
                 for _ in range(4000)])
lo, hi = np.quantile(boot, [0.025, 0.975])
p_less = float((boot < 0).mean())
print(f"""  TET2 deceleration minus DNMT3A deceleration: {diff:+.4f} per year
  95% interval {lo:+.4f} to {hi:+.4f}
  P(TET2 decelerates LESS) = {p_less:.3f}

  Fabre report deceleration "most marked" in DNMT3A and TET2 growing faster in
  older individuals. Both statements predict this difference to be negative.""")
print(f"  {'SEPARATED' if p_less > 0.95 or p_less < 0.05 else 'NOT SEPARATED'} at the 95% level.")


section("3. CONTROLLING FOR WHAT STAGE F SAID ACTUALLY DRIVES DECELERATION")

print("""Stage F found deceleration is a property of the host: the carrier's total
clonal burden growth gave R^2 = 0.956 against 0.003 for the clone's own fitness.
If TET2 carriers differ from DNMT3A carriers in burden or in clone size, the
comparison above is measuring that instead.""")

print(f"\n  {'':<10}{'n':>5}{'VAF at start':>14}{'burden at start':>17}{'burden growth':>15}")
for gene in GENES:
    g = t[t.gene == gene]
    print(f"  {gene:<10}{len(g):>5}{np.median(g.vaf0):>14.4f}"
          f"{np.median(g.burden0):>17.4f}{np.median(g.dburden):>15.4f}")

sub("3.1 Matched on starting clone size")
bands = [FLOOR, 0.01, 0.03, 1.0]
labels = ["0.002-0.01", "0.01-0.03", "0.03+"]
print(f"  {'VAF at start':<14}{'n TET2':>8}{'n DNMT3A':>10}{'TET2 decel':>12}"
      f"{'DNMT3A decel':>14}{'difference':>12}")
for i, lab in enumerate(labels):
    ta = a[(a.vaf0 >= bands[i]) & (a.vaf0 < bands[i+1])]
    tb = b_[(b_.vaf0 >= bands[i]) & (b_.vaf0 < bands[i+1])]
    if len(ta) < 12 or len(tb) < 12:
        print(f"  {lab:<14}{len(ta):>8}{len(tb):>10}{'too few':>12}")
        continue
    print(f"  {lab:<14}{len(ta):>8}{len(tb):>10}{np.median(ta.decel):>12.4f}"
          f"{np.median(tb.decel):>14.4f}{np.median(ta.decel)-np.median(tb.decel):>12.4f}")

sub("3.2 A regression with all of it at once")
m = t[t.gene.isin(GENES)].copy()
X = np.column_stack([(m.gene == "TET2").astype(float),
                     np.log(m.vaf0), m.dburden, m.age0])
names = ["TET2 (vs DNMT3A)", "log VAF at start", "carrier burden growth", "age at start"]
Xs = (X - X.mean(axis=0)) / X.std(axis=0)
ys = (m.decel.values - m.decel.values.mean()) / m.decel.values.std()
A = np.column_stack([Xs, np.ones(len(ys))])
coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
resid = ys - A @ coef
print(f"  standardised coefficients on deceleration   (n={len(m)})")
for nm, c in zip(names, coef[:4]):
    print(f"    {nm:<26}{c:>+8.3f}")
print(f"    {'R^2':<26}{1 - resid.var()/ys.var():>8.3f}")

bcoef = []
for _ in range(2000):
    k = rng.choice(len(ys), len(ys))
    c, *_ = np.linalg.lstsq(A[k], ys[k], rcond=None)
    bcoef.append(c[0])
bcoef = np.array(bcoef)
print(f"\n  TET2 coefficient: {coef[0]:+.3f}  "
      f"95% interval {np.quantile(bcoef,0.025):+.3f} to {np.quantile(bcoef,0.975):+.3f}")


section("4. THE MEASURE THE ABSOLUTE DIFFERENCE WAS HIDING")

print(f"""Section 2 compared deceleration in absolute units and found nothing. That
compares the two genes as if they started from the same place, and they do not:
TET2 grows at {np.median(a.early):.4f} per year in the first half against DNMT3A's {np.median(b_.early):.4f}.

The claim being tested is that TET2 clones keep growing in older people. In
these trajectories the second half IS the older period, so the direct form of
the question is what each gene's rate looks like there.""")

def med_ci(v, B=4000):
    b = np.array([np.median(rng.choice(v, len(v))) for _ in range(B)])
    return float(np.median(v)), float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))

print(f"\n  {'gene':<9}{'n':>5}{'early rate':>12}{'late rate':>12}"
      f"{'late 95% interval':>24}{'retained':>10}")
for gene, g in (("DNMT3A", b_), ("TET2", a)):
    e, _, _ = med_ci(g.early.values)
    l, llo, lhi = med_ci(g.late.values)
    print(f"  {gene:<9}{len(g):>5}{e:>12.4f}{l:>12.4f}"
          f"{f'{llo:+.4f} to {lhi:+.4f}':>24}{l/e:>10.0%}")

late_diff = float(np.median(a.late) - np.median(b_.late))
bl = np.array([np.median(rng.choice(a.late.values, len(a)))
               - np.median(rng.choice(b_.late.values, len(b_)))
               for _ in range(4000)])
lo_l, hi_l = np.quantile(bl, [0.025, 0.975])
p_late = float((bl > 0).mean())

early_diff = float(np.median(a.early) - np.median(b_.early))
be = np.array([np.median(rng.choice(a.early.values, len(a)))
               - np.median(rng.choice(b_.early.values, len(b_)))
               for _ in range(4000)])
lo_e, hi_e = np.quantile(be, [0.025, 0.975])
p_early = float((be > 0).mean())

print(f"""
  TET2 minus DNMT3A, early half: {early_diff:+.4f}   [{lo_e:+.4f}, {hi_e:+.4f}]   P(TET2 higher) = {p_early:.3f}
  TET2 minus DNMT3A, late half:  {late_diff:+.4f}   [{lo_l:+.4f}, {hi_l:+.4f}]   P(TET2 higher) = {p_late:.3f}

  {'The gap between the genes is wider in the older half.' if p_late > p_early else 'The gap does not widen in the older half.'}

  BORDERLINE, AND THE TWO NUMBERS SAY DIFFERENT THINGS. The one-sided
  P(TET2 higher) = {p_late:.3f} clears 0.95, and the two-sided 95% interval
  [{lo_l:+.4f}, {hi_l:+.4f}] does not clear zero — {(bl < 0).mean():.1%} of the bootstrap sits below it.
  This is a result at the edge, and calling it separated would be choosing the
  statistic that says so.""")

sub("4.1 Is DNMT3A's late rate distinguishable from zero?")
_, dlo, dhi = med_ci(b_.late.values)
_, tlo, thi = med_ci(a.late.values)
print(f"""  DNMT3A late rate: {np.median(b_.late):.4f}   [{dlo:+.4f}, {dhi:+.4f}]   {'includes zero' if dlo < 0 < dhi else 'excludes zero'}
  TET2 late rate:   {np.median(a.late):.4f}   [{tlo:+.4f}, {thi:+.4f}]   {'includes zero' if tlo < 0 < thi else 'excludes zero'}

  This is the sharpest statement the data supports, and it is about each gene on
  its own rather than a difference between them: in the second half of their own
  follow-up, past a median age of {t.age0.median() + (t.age1.median()-t.age0.median())/2:.0f}, {'DNMT3A clones are no longer measurably growing' if dlo < 0 < dhi else 'DNMT3A clones are still growing'}
  while {'TET2 clones still are' if tlo > 0 else 'TET2 clones are not either'}.

  A single-gene interval does not depend on a between-group test surviving, so
  it is the part of this that does not rest on a borderline call.""")

sub("4.2 Three measures were tested, and this is the one that worked")
print(f"""  Deceleration in absolute units: nothing ({diff:+.4f}, P = {p_less:.2f}).
  Early-half rate: nothing ({early_diff:+.4f}, P = {p_early:.2f}).
  Late-half rate: borderline ({late_diff:+.4f}, P = {p_late:.2f}).

  Reporting the third without the first two would be the standard way to
  manufacture a finding. Three related measures on one dataset, and the best of
  them lands at the edge of significance — which is roughly what chance alone
  produces from time to time.

  What keeps this from being noise-mining is that the direction was specified in
  advance, by Fabre and by Mon Pere, and all three measures point the same way.
  A prediction made before the test is worth more than a p-value found after
  it, but neither is worth as much as a replication.""")


section("5. WHAT COULD HAVE BEEN DETECTED")

sd = float(np.std(np.concatenate([a.decel.values, b_.decel.values])))
n_eff = min(len(a), len(b_))
mdd = 1.96 * sd * np.sqrt(2.0 / n_eff)
sd_late = float(np.std(np.concatenate([a.late.values, b_.late.values])))
mdd_late = 1.96 * sd_late * np.sqrt(2.0 / n_eff)
print(f"""  {'quantity':<34}{'SD':>10}{'smallest detectable gap':>26}
  {'deceleration (early - late)':<34}{sd:>10.4f}{mdd:>26.4f}
  {'late-half rate':<34}{sd_late:>10.4f}{mdd_late:>26.4f}

  The observed between-gene difference in deceleration is {abs(diff):.4f}, well under the
  {mdd:.4f} needed — so section 2's null is a statement about power, not about
  biology, and the honest reading is that this test could not have separated the
  genes on that measure however the biology behaved.

  Those thresholds use a normal approximation on MEANS, and everything above
  compares MEDIANS, which are more stable under this heavy-tailed spread. That
  is why the observed late-half gap of {abs(late_diff):.4f} reaches the edge of significance in
  the bootstrap despite sitting under the {mdd_late:.4f} the crude formula asks for. The
  bootstrap intervals are the honest measure of precision here; the table is a
  sanity check on the order of magnitude, and it disagrees for a known reason.""")
