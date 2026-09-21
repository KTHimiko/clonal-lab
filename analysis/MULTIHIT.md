# Stage G — multihit, and an instrument that does not measure phase

Reproduce with `.venv/bin/python analysis/11_multihit.py`.

Mon Père et al. 2026 report that DNMT3A is enriched for single-hit clones while
TET2 is enriched for multihit evolution, and that the fittest clones emerge late
in life. That supplies a third explanation for the TET2 age effect
([`TET2_AGE.md`](TET2_AGE.md)): a clone that picks up a second driver
accelerates, and under a single-mutation model that is indistinguishable from
rising fitness — gene-specific in exactly the way the observation is.

Fabre's data has person identifiers, so the test looked cheap.

---

## The distinction the test lives or dies on

Two drivers in one person are not one thing. They are two, with **opposite**
consequences:

| | | |
|---|---|---|
| **cis** | both mutations in the same cells | fitness adds, the clone **accelerates** |
| **trans** | two independent clones | they compete, and stage F measured what that does: **everyone slows** |

Counting "people with two drivers" adds those signals and measures their sum,
which is nothing in particular.

**The proposed separation.** If two variants sit in the same cells their VAFs
move together, so the log of their ratio stays flat as the clone grows. Two
independent clones drift apart at a rate equal to the difference in their growth
rates. So: fit each clone's rate, take the difference, and call a pair *cis* when
it is near zero.

---

## Why the first version was guaranteed to find nothing

It called a pair cis when the rates differed by less than 0.01 per year, and the
excess over a between-person null was +0.006. That looked like a clean negative.

It was not a negative, it was a **coin**. A VAF read at depth *d* carries
variance `(1-VAF)/(VAF·d)` on the log scale, and with a median depth of 1,140×
over four timepoints the standard error on the *difference* of two rates is:

| VAF band | measurements | SE of the rate difference |
|---|---|---|
| 0.002 – 0.01 | 1,610 | **0.074** /yr |
| 0.01 – 0.05 | 1,000 | 0.031 |
| 0.05 – 0.20 | 263 | 0.014 |
| 0.20 – 1.00 | **82** | 0.006 |

**The threshold was seven times smaller than the noise in the band holding two
thirds of the data.** Nothing could have been detected there, and the null
result carried no information about biology.

---

## The properly weighted version

Each clone's rate is refit by weighted least squares with per-point variance from
read counts, so every pair carries its own standard error and can be asked
whether it is precise enough to answer anything.

| SE of the rate difference | within-person pairs | null pairs |
|---|---|---|
| better than 0.10 | 580 | 5,661 |
| better than 0.05 | 259 | 2,473 |
| **better than 0.02** | **13** | 150 |
| better than 0.01 | 1 | 9 |
| total | 613 | 6,000 |

**13 of 613 pairs — 2% — are measured well enough that "these two rates agree"
means anything.**

Among those 13: **46% agree within error against 28% of pairs that are
independent by construction.** Excess +0.182, 95% interval −0.089 to +0.472,
P(excess > 0) = 0.90.

That is the shape the hypothesis predicts, at a sample size that cannot support
it. Recorded, not leaned on.

Everything downstream inherits the same limit:

| question | what the data gives |
|---|---|
| is TET2 enriched for cis pairs? | +0.007, interval spans zero, P = 0.81 |
| do cis clones grow faster? | 12 clones, not reportable |
| does cis prevalence rise with age? | 1–2 people per age band |

---

## What would answer it

1. **Deeper sequencing.** The standard error falls as `1/√depth`, so 1,140× →
   ~11,000× moves the 0.01–0.05 band from 0.031 to 0.010 per year and makes most
   pairs usable.
2. **More timepoints, longer window.** Both act through the spread of ages,
   already near its practical limit at 13 years.
3. **Phasing — the actual answer.** Whether two mutations sit in the same cell is
   not a question about correlated frequencies. It is a question about whether
   they are on the same DNA molecule. Single-cell colonies or long reads settle
   it directly.

**That third point is not a limitation of this analysis. It is the reason Mon
Père used single-cell-derived colonies to make the claim.** We tried to reach a
phasing conclusion with an instrument that does not measure phase.

---

## Where this leaves the TET2 question

Unchanged, and now for a fourth distinct reason.

| attempt | why it did not resolve |
|---|---|
| Watson pooled | effect vanishes when stratified by cohort; 87% of the estimate was composition |
| Watson, power | 36 variants; ceiling 0.79 even if the ramp is real |
| ramp vs wild-type decline | arithmetically the same model in relative fitness |
| **multihit, here** | **needs phase; bulk VAF resolves it for 2% of pairs** |

The multihit explanation remains the most plausible of the three — it is
gene-specific from measured biology and requires no new assumption — and it
remains untested by us.
