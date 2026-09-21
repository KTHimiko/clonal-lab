#!/usr/bin/env python3
"""
Stage G — is the TET2 age effect a second mutation rather than a rising fitness?

THE HYPOTHESIS, FROM MON PERE ET AL. 2026
DNMT3A is enriched for single-hit clones; TET2, ASXL1, JAK2, SF3B1 and SRSF2
for multihit evolution, and "the fittest clones emerge predominantly later in
life in accordance with a multistep evolutionary process". A TET2 clone that
picks up a second driver accelerates, and under a single-mutation model that is
indistinguishable from "TET2 fitness rises with age" — gene-specific in exactly
the way the observation is.

THE DISTINCTION THAT DECIDES THE TEST
Two drivers in one person are not one thing. They are two, with opposite
consequences:

  CIS    both mutations in the same cells. Fitness adds, the clone ACCELERATES.
  TRANS  two independent clones. They compete, and stage F measured what that
         does: everyone slows.

A test that counts "people with two drivers" adds those signals and measures
their sum, which is nothing in particular. They have to be separated first.

HOW BULK VAF CAN SEPARATE THEM
If two variants sit in the same cells, their VAFs move together: the ratio
between them stays constant as the clone grows, whatever that ratio is. Two
independent clones grow at their own rates, so the log of their VAF ratio drifts
with time at a rate equal to the difference in growth rates.

So the statistic is the SLOPE of log(VAF_a / VAF_b) against age. Near zero means
same clone. The separation is probabilistic, not definitive — two independent
clones of near-equal fitness also give a flat ratio — which is why the null below
is built from pairs that are independent by construction.

Usage:  .venv/bin/python analysis/11_multihit.py
"""

import sys
from pathlib import Path
from itertools import combinations
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


MIN_SHARED = 3        # timepoints a pair must share
MIN_SPAN = 5.0        # years the pair must span
FLOOR = 0.002         # below this the VAF ratio is mostly read noise

d = pd.read_csv(FABRE, sep="\t")
d["clone"] = (d.SardID.astype(str) + "|" + d.Gene + "|"
              + d.START.astype(str) + "|" + d.ALT.astype(str))
d = d[d.VAF >= FLOOR]

# per-clone trajectories, carrying sequencing depth so the fit can be weighted
TRAJ = {c: g[["Age", "VAF", "TOTALcount"]].sort_values("Age").values
        for c, g in d.groupby("clone")}
INFO = d.sort_values("Phase").groupby("clone").agg(
    person=("SardID", "first"), gene=("Gene", "first"),
    pos=("START", "first"), age0=("Age", "first"),
    vafmax=("VAF", "max")).to_dict("index")
traj, info = TRAJ, INFO


def weighted_rate(c):
    """
    A clone's growth rate and the uncertainty the sequencer imposes on it.

    A VAF is estimated from a finite number of reads, so log(VAF) at one
    timepoint carries variance (1-VAF)/(VAF*depth). That variance is enormous
    for small clones and small for large ones, which is why the fit is weighted
    and why the standard error below decides what this analysis can conclude.
    """
    a = TRAJ[c]
    if len(a) < MIN_SHARED or np.ptp(a[:, 0]) < MIN_SPAN:
        return np.nan, np.nan
    age, vaf, depth = a[:, 0], a[:, 1], a[:, 2]
    var = (1.0 - vaf) / (vaf * depth)
    w = 1.0 / var
    x = age - np.average(age, weights=w)
    sxx = np.sum(w * x * x)
    if sxx <= 0:
        return np.nan, np.nan
    slope = np.sum(w * x * np.log(vaf)) / sxx
    return float(slope), float(np.sqrt(1.0 / sxx))


RATE, RATE_SE = {}, {}
for c in TRAJ:
    RATE[c], RATE_SE[c] = weighted_rate(c)


def pair_stat(ca, cb):
    """Rate difference and its standard error, propagated from read noise."""
    ra, rb = RATE[ca], RATE[cb]
    sa, sb = RATE_SE[ca], RATE_SE[cb]
    if not all(np.isfinite(v) for v in (ra, rb, sa, sb)):
        return None
    diff = ra - rb
    se = float(np.hypot(sa, sb))
    return diff, se


section("1. HOW PRECISE CAN THIS TEST BE, BEFORE ASKING WHAT IT FINDS")

print("""The first version of this analysis called a pair "same clone" when the two
growth rates differed by less than 0.01 per year, and found nothing. That was
guaranteed, and the reason is measurable rather than biological.

A VAF read at depth d carries variance (1-VAF)/(VAF*d) on the log scale. With a
median depth of ~1,140x and four timepoints over ~13 years, the standard error
on the DIFFERENCE of two clones' growth rates runs:

  VAF 0.002-0.01   SE ~ 0.074 per year     <- threshold was 7x smaller than noise
  VAF 0.01-0.05    SE ~ 0.031
  VAF 0.05-0.20    SE ~ 0.014
  VAF 0.20-1.00    SE ~ 0.006

Two thirds of the measurements sit in the top row. A threshold of 0.01 applied
there is not a test, it is a coin. So the question has to be asked the other way
round: which pairs are measured precisely enough to answer it at all?""")

within = []
for person, g in d.groupby("SardID"):
    clones = sorted(g.clone.unique())
    for ca, cb in combinations(clones, 2):
        r = pair_stat(ca, cb)
        if r is None:
            continue
        diff, se = r
        within.append(dict(person=person, ca=ca, cb=cb, diff=diff, se=se,
                           z=diff / se,
                           gene_a=INFO[ca]["gene"], gene_b=INFO[cb]["gene"],
                           same_gene=INFO[ca]["gene"] == INFO[cb]["gene"],
                           vaf_min=min(INFO[ca]["vafmax"], INFO[cb]["vafmax"])))
within = pd.DataFrame(within)

rng = np.random.default_rng(11)
clones = [c for c in TRAJ if np.isfinite(RATE[c]) and np.isfinite(RATE_SE[c])]
null = []
tries = 0
while len(null) < 6000 and tries < 400_000:
    tries += 1
    i, j = rng.choice(len(clones), 2, replace=False)
    ca, cb = clones[i], clones[j]
    if INFO[ca]["person"] == INFO[cb]["person"]:
        continue
    r = pair_stat(ca, cb)
    if r:
        null.append((r[0], r[1]))
null = pd.DataFrame(null, columns=["diff", "se"])
null["z"] = null["diff"] / null.se

sub("1.1 What fraction of pairs can resolve a real difference at all")
print(f"  {'SE of the rate difference':<30}{'within person':>16}{'null':>10}")
for cut in (0.10, 0.05, 0.02, 0.01):
    print(f"  {'better than ' + format(cut, '.2f'):<30}"
          f"{(within.se < cut).sum():>16}{(null.se < cut).sum():>10}")
print(f"  {'total pairs':<30}{len(within):>16}{len(null):>10}")

PRECISE = 0.02
wp = within[within.se < PRECISE]
np_ = null[null.se < PRECISE]
print(f"""
  Only {len(wp)} of {len(within)} within-person pairs ({len(wp)/len(within):.0%}) are measured well enough
  that "the two rates agree" means anything. Everything below uses those.""")

sub("1.2 Among the precise pairs, is there an excess of agreement?")
print(f"  {'':<22}{'n':>6}{'|z| < 2':>10}{'rate':>9}")
print(f"  {'within person':<22}{len(wp):>6}{(wp.z.abs() < 2).sum():>10}{(wp.z.abs() < 2).mean():>9.3f}")
print(f"  {'between people':<22}{len(np_):>6}{(np_.z.abs() < 2).sum():>10}{(np_.z.abs() < 2).mean():>9.3f}")
exc = (wp.z.abs() < 2).mean() - (np_.z.abs() < 2).mean()
# named distinctly: `b` is reused by a later bootstrap, and reading the wrong
# one here reported P = 0.01 for a quantity that is 0.90
excess_boot = np.array([rng.choice((wp.z.abs() < 2).values, len(wp)).mean()
                        - rng.choice((np_.z.abs() < 2).values, len(np_)).mean()
                        for _ in range(4000)])
P_EXCESS = float((excess_boot > 0).mean())
lo, hi = np.quantile(excess_boot, [0.025, 0.975])
print(f"""
  excess within people: {exc:+.3f}   95% interval {lo:+.3f} to {hi:+.3f}
  P(excess > 0) = {P_EXCESS:.3f}

  Two mutations in the same cells must give z near zero, because their VAFs are
  the same measurement of the same cells. An excess of |z| < 2 within people,
  over pairs that are independent by construction, is the signature of cis
  pairs existing and being detectable.""")

within["cis"] = (within.z.abs() < 2) & (within.se < PRECISE)

section("2. IS TET2 ENRICHED FOR THEM, AS MON PERE REPORT?")

null_rate = float(((null.se < PRECISE) & (null.z.abs() < 2)).mean())

rows = []
for gene in ("TET2", "DNMT3A"):
    m = (within.gene_a == gene) | (within.gene_b == gene)
    g = within[m]
    if len(g) < 20:
        continue
    rows.append((gene, len(g), int(g.cis.sum()), g.cis.mean(),
                 int(g[g.same_gene].shape[0]), int(g[g.same_gene].cis.sum())))

print(f"  {'gene':<9}{'pairs':>7}{'flat':>6}{'rate':>8}{'expected':>10}"
      f"{'same-gene pairs':>17}{'of those flat':>15}")
for gene, n, c, rate, sg, sgc in rows:
    print(f"  {gene:<9}{n:>7}{c:>6}{rate:>8.3f}{null_rate:>10.3f}{sg:>17}{sgc:>15}")

if len(rows) == 2:
    (g1, n1, c1, r1, _, _), (g2, n2, c2, r2, _, _) = rows
    boot = []
    a = within[(within.gene_a == g1) | (within.gene_b == g1)].cis.values
    b = within[(within.gene_a == g2) | (within.gene_b == g2)].cis.values
    for _ in range(4000):
        boot.append(rng.choice(a, a.size).mean() - rng.choice(b, b.size).mean())
    boot = np.array(boot)
    lo, hi = np.quantile(boot, [0.025, 0.975])
    print(f"""
  {g1} minus {g2}: {r1-r2:+.3f}   95% interval {lo:+.3f} to {hi:+.3f}
  P({g1} higher) = {(boot > 0).mean():.3f}   """
          f"{'SEPARATED' if (boot>0).mean() > 0.95 or (boot>0).mean() < 0.05 else 'NOT SEPARATED'}")


section("3. DO FLAT-RATIO PAIRS GROW FASTER, AS THE HYPOTHESIS REQUIRES?")

rates = RATE
in_pair = set(within.ca) | set(within.cb)
cis_clones = set(within[within.cis].ca) | set(within[within.cis].cb)
trans_clones = (set(within[~within.cis].ca) | set(within[~within.cis].cb)) - cis_clones
solo = {c for c in traj if c not in in_pair
        and sum(1 for x in traj if info[x]["person"] == info[c]["person"]) == 1}

def describe(name, keys):
    v = np.array([rates[c] for c in keys if np.isfinite(rates.get(c, np.nan))])
    if v.size < 15:
        print(f"  {name:<34}{v.size:>5}   too few")
        return None
    b = np.array([np.median(rng.choice(v, v.size)) for _ in range(4000)])
    lo, hi = np.quantile(b, [0.025, 0.975])
    print(f"  {name:<34}{v.size:>5}{np.median(v):>11.4f}   {lo:.4f} to {hi:.4f}")
    return v

print(f"  {'group':<34}{'n':>5}{'median rate':>11}   95% interval")
v_solo = describe("only driver in that person", solo)
v_cis = describe("in a flat-ratio pair (cis-like)", cis_clones)
v_trans = describe("in a drifting pair (independent)", trans_clones)

print(f"""
  THE HYPOTHESIS PREDICTS cis-like > solo, because two drivers in one cell add
  fitness. STAGE F PREDICTS independent < solo, because competing clones slow
  each other. The two predictions point in opposite directions, which is what
  makes the split worth making.""")


section("4. DOES IT RISE WITH AGE — THE MECHANISM FOR APPARENT AGE-DEPENDENCE?")

first = d.sort_values("Phase").groupby("SardID").Age.first().rename("age0")
ppl = within.groupby("person").cis.any().rename("has_cis").reset_index()
ppl = ppl.merge(first, left_on="person", right_index=True)
bands = pd.cut(ppl.age0, [50, 65, 70, 75, 110],
               labels=["<65", "65-70", "70-75", "75+"])
tab = ppl.groupby(bands, observed=True).agg(
    people=("person", "size"), with_cis=("has_cis", "sum"))
tab["rate"] = tab.with_cis / tab.people
print(f"  {'age at entry':<14}{'people':>8}{'with a flat pair':>18}{'rate':>8}")
for idx, r in tab.iterrows():
    print(f"  {str(idx):<14}{int(r.people):>8}{int(r.with_cis):>18}{r.rate:>8.2f}")


section("5. WHAT THIS ESTABLISHES")

precise_frac = len(wp) / len(within)
print(f"""THE TEST CANNOT BE RUN ON THIS DATA, and the reason is specific enough to act
on rather than a shrug.

Distinguishing two mutations in one cell from two clones in one marrow requires
their growth rates to agree to better than the difference competition would
produce — about 0.02 per year. Read noise at a median depth of 1,140x gives that
precision only when a clone sits above VAF ~0.05 at most timepoints. **{len(wp)} of {len(within)}
within-person pairs ({precise_frac:.0%}) clear that bar.**

Everything downstream inherits it:

  TET2 vs DNMT3A enrichment   n too small, interval spans zero
  do cis pairs grow faster    12 clones, not reportable
  does it rise with age       1 to 2 people per age band

The suggestive part is worth recording without leaning on it: among the 13
pairs that ARE precise, {(wp.z.abs() < 2).mean():.0%} agree within error against {(np_.z.abs() < 2).mean():.0%} of independent
pairs, P(excess > 0) = {P_EXCESS:.2f}. That is the shape the hypothesis predicts, at a
sample size that cannot support it.

WHAT WOULD ANSWER IT, in increasing order of how well it works:

  1. Deeper sequencing. The standard error falls as 1/sqrt(depth), so moving
     from 1,140x to ~11,000x would bring the VAF 0.01-0.05 band from 0.031 to
     0.010 per year and make most pairs usable.
  2. More timepoints over a longer window. Both enter through the spread of
     ages, which is already near its practical limit here at 13 years.
  3. PHASING, which is the actual answer. Whether two mutations sit in the same
     cell is not a question about correlated frequencies — it is a question
     about whether they are on the same DNA molecule. Single-cell colonies or
     long reads settle it directly, and inferring it from bulk VAF correlation
     only works in the highest-frequency corner of the data.

That third point is not a limitation of this analysis, it is the reason Mon Pere
used single-cell-derived colonies to make the claim in the first place. We tried
to reach a phasing conclusion with an instrument that does not measure phase.""")
