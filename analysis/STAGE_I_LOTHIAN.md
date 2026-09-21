# Stage I — a replication that could not be run

Reproduce with `.venv/bin/python analysis/13_lothian_replication.py`.
Plan fixed in advance: [`PREREGISTRATION_STAGE_I.md`](PREREGISTRATION_STAGE_I.md),
committed as `643bd60` before any rate from this cohort was computed.

Stage H found TET2 clones still growing after 75 while DNMT3A clones were not,
at a borderline significance and with three measures tested. It needed a
replication. Robertson 2022's Lothian Birth Cohorts looked like the one.

**It is not.** This document records why, and two things the attempt established
anyway.

---

## 1. The pre-registration failed in a visible way, which is the point

It fixed the estimator, the thresholds, the time axis and the bootstrap. It named
the GEO accession, **GSE178936**, and did not name which of its two supplementary
files to read.

The first run took the 1% VAF file and produced an estimate: −0.0116, interval
−0.0443 to +0.0772. The pre-specified sanity check rejected it:

- LBC1921 clones had a **median growth rate of −0.045 per year** — shrinking
- the most frequent genes were **NF1, BCORL1, KMT2A**

Real clones do not shrink en masse, and those are not clonal-haematopoiesis
drivers at those frequencies.

| file | rows | genes | median depth | leading genes |
|---|---|---|---|---|
| 1% VAF | 8,566 | 69 | 1,110× | NF1, BCORL1, KMT2A, DNMT3A |
| **2% VAF** | **255** | **23** | **2,153×** | **DNMT3A, TET2, JAK2, CDKN2A** |

The paper reports a median depth of 2,153×. Only the second file matches, and
only the second has the gene composition of CHIP. The 1% file is caller output;
the 2% file is the cohort.

**That estimate is discarded.** It is recorded here so the discard is visible,
not so the number is available.

---

## 2. On the right file the test is not underpowered, it is not runnable

96 distinct clones. With at least three waves:

| gene | clones |
|---|---|
| DNMT3A | 25 |
| **TET2** | **8** |
| JAK2 | 5 |

The pre-registration computed power as 0.25 assuming 49 TET2 and 129 DNMT3A
clones, which is what the wrong file appeared to offer. There are eight.

A median of eight values has no useful interval, and a between-group comparison
with eight on one side is not an estimate of anything. **No number is reported**,
because reporting one invites it being quoted.

---

## 3. What the attempt established anyway

**Stage H had 0.56 power.** Computed while planning this, by resampling stage H's
own distributions at its own sample sizes. A borderline positive from a
56%-powered study is the profile of a result that does not replicate. Recorded as
a retrospective correction in
[`TET2_LONGITUDINAL.md`](TET2_LONGITUDINAL.md).

**A plan that cannot fail visibly is not a plan.** The pre-registration did not
prevent the mistake — it made the mistake catchable, by committing a sanity check
alongside the test. Without it, −0.0116 would have been reported as a failed
replication of stage H, and the failure would have been an artifact of reading
the wrong file.

**What would be needed.** Power scales roughly with the square root of sample
size, so moving stage H's 0.56 to a conventional 0.80 needs on the order of 200
to 230 TET2 clones with four or more timepoints — about twice Fabre's 114, and
twenty-five times what Lothian offers.

---

## Where the TET2 question stands

Unchanged from stage H, and now held more loosely.

TET2 clones appear to still be growing after 75 when DNMT3A clones are not, in
**one cohort, at 56% power, with three measures tested and one reaching the
edge.** It has not been replicated, and this cohort cannot do it.
