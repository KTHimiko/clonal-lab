#!/usr/bin/env python3
"""
Stage I — the Lothian replication, attempted and abandoned for a stated reason.

The plan, the power and the known weaknesses were committed before any rate
from this cohort was computed (analysis/PREREGISTRATION_STAGE_I.md, 643bd60).

WHAT HAPPENED
The pre-registration named the GEO accession but not which of its two
supplementary files to use. The first run took the 1% VAF file, produced an
estimate, and the pre-specified sanity check rejected it: median clone growth in
LBC1921 came out at -0.045 per year, clones shrinking, with NF1, KMT2A and
BCORL1 as the most frequent genes. Those are not CHIP drivers at those
frequencies, and real clones do not shrink en masse.

The 1% file is the low-threshold caller output. The 2% file is the CHIP-grade
set: 23 genes led by DNMT3A and TET2, median depth 2,188x, matching the 2,238x
the paper reports. That is the one the authors analysed.

On the correct file the test is not underpowered. It is impossible.

Usage:  .venv/bin/python analysis/13_lothian_replication.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LBC2 = (ROOT / "reference/data/robertson2022"
             / "GSE178936_LBC_ARCHER.2PCT_VAF.Feb22.non-synonymous.tsv")
LBC1 = (ROOT / "reference/data/robertson2022"
             / "GSE178936_LBC_ARCHER.1PCT_VAF.Feb22.non-synonymous.tsv")

def section(t):
    print(f"\n{'='*76}\n{t}\n{'='*76}")

FLOOR, MIN_WAVES = 0.002, 3


def load(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d[d.AF >= FLOOR].copy()
    d["cohort"] = d.participant_id.str.extract(r"(LBC\d+)")
    d["clone"] = (d.participant_id + "|" + d.PreferredSymbol + "|"
                  + d.chromosome.astype(str) + ":" + d.position.astype(str)
                  + d.mutation.astype(str))
    return d


section("1. WHICH FILE IS THE COHORT, AND HOW THE FIRST RUN GOT IT WRONG")

for tag, path in (("1% VAF (caller output)", LBC1), ("2% VAF (CHIP-grade)", LBC2)):
    d = load(path)
    top = ", ".join(d.PreferredSymbol.value_counts().head(4).index)
    print(f"  {tag:<26} {len(d):>6} rows  {d.PreferredSymbol.nunique():>3} genes  "
          f"depth {d.DP.median():>6,.0f}x   top: {top}")

print("""
  The paper reports a median depth of 2,153x and a mean of 2,238x. Only the
  second file matches, and only the second has the gene composition of clonal
  haematopoiesis. The first run used the first file and the sanity check caught
  it: LBC1921 clones came out shrinking at 4.5% per year.""")


section("2. ON THE CORRECT FILE, WHAT IS THERE")

d = load(LBC2)
w = d.groupby("clone").agg(waves=("wave", "nunique"),
                           gene=("PreferredSymbol", "first"),
                           cohort=("cohort", "first"))
print(f"  distinct clones: {len(w)}   median depth {d.DP.median():,.0f}x")
print(f"\n  clones by number of waves:")
for k, v in w.waves.value_counts().sort_index().items():
    print(f"    {k} waves: {v}")

u = w[w.waves >= MIN_WAVES]
print(f"\n  with at least {MIN_WAVES} waves, by gene:")
for gene, n in u.groupby("gene").size().sort_values(ascending=False).head(6).items():
    print(f"    {gene:<10}{n:>4}")

n_t = int((u.gene == "TET2").sum())
n_d = int((u.gene == "DNMT3A").sum())


section("3. THE TEST IS NOT UNDERPOWERED, IT IS NOT RUNNABLE")

print(f"""  TET2 clones available: {n_t}
  DNMT3A clones available: {n_d}

  The pre-registration set power at 0.25 assuming 49 TET2 and 129 DNMT3A clones,
  which is what the 1% file appeared to offer. The CHIP-grade file has {n_t}.

  A median of {n_t} values has no useful interval, and a between-group comparison
  with {n_t} on one side is not an estimate of anything. NO NUMBER IS REPORTED HERE,
  because reporting one invites it being quoted.

  The estimate the first run produced on the wrong file — -0.0116, interval
  -0.0443 to +0.0772 — is DISCARDED. It is recorded in the git history and in
  this comment so the discard is visible, not so the number is available.""")


section("4. WHAT THE ATTEMPT ESTABLISHED ANYWAY")

print("""  1. STAGE H HAD 0.56 POWER. Computed for the first time while planning this,
     by resampling stage H's own distributions at its own sample sizes. A
     borderline positive from a study with 56% power is the profile of a result
     that does not replicate. That is a retrospective correction to stage H and
     is recorded in analysis/TET2_LONGITUDINAL.md.

  2. THE PRE-REGISTRATION WORKED, by failing in a visible way. It fixed the
     estimator, the thresholds and the bootstrap, and it did not fix which file
     to read. The sanity check it also specified is what caught that. A plan
     that cannot fail visibly is not a plan.

  3. WHAT WOULD BE NEEDED. Power scales roughly with the square root of sample
     size, so moving stage H's 0.56 to a conventional 0.80 needs something near
     200 to 230 TET2 clones with four or more timepoints — about twice Fabre's
     114, and twenty times what Lothian offers.

  The honest position on the TET2 question is unchanged from stage H, and is now
  held more loosely: TET2 clones appear to still be growing after 75 when DNMT3A
  clones are not, in one cohort, at 56% power, with three measures tested and one
  reaching the edge. It has not been replicated, and this cohort cannot do it.""")
