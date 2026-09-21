# What this project found, in one read

Ten stages, twelve analyses, two datasets. This is the whole thing in order,
including what was corrected along the way, so nothing has to be reconstructed
from the commit log.

The per-stage documents in [`analysis/`](analysis/) hold the detail. Read this
first.

---

## The question it started with, and the one it ended up answering

**Started:** how strong is the fitness advantage that lets a mutated blood stem
cell outcompete its neighbours?

**Ended:** how much can you trust each way of measuring that advantage?

The change was not a loss of focus. Three published methods disagree by up to a
factor of three on the same gene, so any new estimate is one more number in a
disagreement nobody can adjudicate — while a measurement of the *methods* is
something no dataset can produce and a simulator can.

---

## The results, in the order they were established

### 1. The simulator is correct where it can be checked

Reproduces the exact Moran fixation probability, neutral and under selection,
across population sizes and starting points, all within 3 standard errors.
([`02_validation.py`](analysis/02_validation.py))

**And the Moran approximation is `s/(1+s)`, not `2s`.** The widely quoted `2s`
comes from Wright–Fisher, which has a different offspring variance. Using the
wrong one doubles every fitness estimate.

### 2. Exact beats approximate, and here it is also faster

The multi-clone model first used tau-leaping, which let a one-cell clone die and
give birth in the same slice and inflated survival by 3 standard errors.
Sampling the birth–death transition exactly removed the bias **and** allowed a
time step five times larger. ([`model/population.py`](model/population.py))

### 3. Fitness is identifiable from a VAF spectrum; the mutation rate is not

A grid scored by approximate Bayesian computation puts `s` near 0.13 per year
with every accepted point at the same value, while the mutation rate stays
spread over 91% of the grid. The distance surface shows why: a vertical ridge.
Shape responds to fitness; counts respond to mutation rate, and counts need a
screening denominator the data does not carry. ([`FINDINGS.md`](analysis/FINDINGS.md),
[`04_abc.py`](analysis/04_abc.py))

### 4. Two of the three designs cannot be inverted at all

Stage E found the age regression is non-monotonic in `s`: it rises then falls, so
one observed slope matches several truths. Stage J found the same of the
longitudinal estimator — the median growth rate peaks near `s` = 0.12 and decays
toward zero, because a fitter clone is already saturating when follow-up starts.

The VAF spectrum is monotonic across the whole grid. **It is the only one of the
three that can be fitted to recover a fitness**, which is a sharper reason to
prefer it than any comparison of accuracy.
([`STAGE_J_FABRE_GRID.md`](analysis/STAGE_J_FABRE_GRID.md))

### 5. Study designs recover very different amounts of a known truth

Given a simulated cohort whose fitness is set by construction:

| design | fraction of the truth recovered |
|---|---|
| regression of log(VAF) on the host's age | **18%** |
| longitudinal follow-up over 13 years | **43%** |
| fit to the detected VAF spectrum | **103%** |

The age regression is worse than biased — it is **non-monotonic** in `s`, so it
has no inverse: one observed slope is compatible with several very different
truths. ([`METHOD_BIAS.md`](analysis/METHOD_BIAS.md))

### 6. That ordering holds in 394 real people

Fabre's SardiNIA cohort allows all three designs on the same individuals.

| design | estimate |
|---|---|
| age regression | 0.0024 /yr |
| longitudinal | 0.0455 /yr |
| spectrum fit | 0.110 /yr |

Spectrum ÷ longitudinal: **2.42 observed, 2.40 predicted**, and 2.42 between the
two published papers. Three routes to the same magnitude.

**The part that needs no model at all:** in 770 real trajectories, restricting to
clones already detectable at enrolment drops the measured growth rate from
0.0703 to 0.0455 per year — a 35% loss, non-overlapping intervals.
([`STAGE_F_REAL_COHORT.md`](analysis/STAGE_F_REAL_COHORT.md))

### 7. Clones slow because the marrow fills, not because they tire

Regressing per-clone deceleration on three candidates: the carrier's total
fitness-weighted burden growth gives **R² = 0.956**, the clone's own size 0.271,
its own fitness **0.003**. Deceleration is a property of the host.

Mon Père et al. (2026) derive the same thing analytically — growth rate is
proportional to the difference between a clone's own fitness and the *average
fitness of the population*. ([`OPEN_PROBLEMS.md`](analysis/OPEN_PROBLEMS.md))

### 8. Three independent methods agree once compared like for like

| source | method | DNMT3A, per year |
|---|---|---|
| Watson 2020 | VAF spectrum | 0.112 – 0.160 |
| Mitchell 2022 | phylogenetic | 0.167 – 0.200 |
| this project | VAF spectrum | 0.120 – 0.130 |
| Fabre 2022 | longitudinal | 0.062 |

The phylogenetic estimate shares no machinery with a frequency spectrum and
lands with the spectrum fits. The longitudinal value is the outlier, low, by the
amount stages E and F predict.

### 9. TET2 clones are still growing after 75; DNMT3A clones are not

234 trajectories, each clone split at its own midpoint so it serves as its own
control:

| gene | early half | late half | late 95% interval |
|---|---|---|---|
| DNMT3A | 0.0509 | 0.0096 | −0.0174 to +0.0284 *(includes zero)* |
| TET2 | 0.0711 | 0.0401 | **+0.0133 to +0.0598** |

The between-gene comparison is **borderline** (one-sided P = 0.968, two-sided
interval grazes zero), **three measures were tested**, and the study had
**0.56 power** — computed afterwards, while planning a replication. A borderline
positive at 56% power is the profile of a result that does not replicate. Hold it
loosely. ([`TET2_LONGITUDINAL.md`](analysis/TET2_LONGITUDINAL.md),
[`STAGE_I_LOTHIAN.md`](analysis/STAGE_I_LOTHIAN.md))

---

## What was corrected, and by what

The corrections are more useful than the results, because each names a mistake
that is easy to repeat.

| what was wrong | what caught it | what it cost |
|---|---|---|
| `N` = 50,000 stem cells | reading Watson and Mitchell | fit improved 2.4×; `s` moved 0.10 → 0.13 |
| tau-leaping let dead clones give birth | a validation test against the exact formula | +3σ bias in survival |
| pooling ZinkWGS, which sees nothing below VAF 0.24 | per-cohort detection limits | TET2 estimate had been an implausible 0.346 |
| "Watson is the cross-sectional estimate" | reading their methods | the whole framing of stage E |
| comparing R882H against a gene average | Fabre stating the comparable pair themselves | a ratio of 3.0 that should be 2.42 |
| comparing Mitchell's inferred spectrum against detected clones | opening their own data file | a third method appeared to contradict us; it agrees |
| "the saturation bias grows with `s`" | regressing deceleration on three candidates | right result, wrong mechanism |
| "no prior art for this" | reading the repository of the paper being argued with | there is adjacent prior art |
| a cis/trans threshold 7× smaller than the read noise | computing the noise instead of trusting the null | a false negative that looked clean |
| a bootstrap variable reused by two sections | the number contradicting itself | P = 0.01 reported for a quantity that is 0.90 |
| stage H reported without its power | planning a replication forced the calculation | 0.56, which downgrades a borderline positive |
| reading the caller output instead of the cohort | a sanity check committed with the pre-registration | a spurious failed replication, caught before it was believed |
| banding a ratio whose denominator has no inverse | a diagnostic table the write-up demanded | a 36-minute cluster run that could not have answered its own question |

**The pattern underneath almost all of them:** before comparing two numbers,
establish that they measure the same thing. Hotspot against gene average.
Detected clones against the spectrum they came from. One cohort against six.
Every disagreement this project chased turned out to be that, and never biology.

---

## What the data could not answer, and why

Four questions died, each on a different missing column. Listing them is the
honest inventory of the project's limits.

| question | what was missing |
|---|---|
| the mutation rate | a screening denominator — how many people were tested, not just how many carried something |
| an unbiased age-based estimator | the same denominator; Watson's prevalence-vs-age estimator works and we cannot run it |
| does TET2 fitness rise with age (cross-sectionally) | power — 36 variants, ceiling 0.79 even if the effect were real |
| are two mutations in the same cell | **phase** — bulk VAF resolves it for 2% of pairs; this needs single-cell colonies or long reads |

The recurring lesson: **an aggregated table answers questions about variants.
Every question about people needs the cohorts it was built from.**

---

## A note on terminology

Earlier documents call the slowdown "saturation". The established term is
**clonal interference**, and it is more accurate: a clone is not approaching a
ceiling of its own, it is being overtaken by fitter neighbours competing for the
same niche. Both words appear in the per-stage documents; they refer to the same
mechanism, and the model implements it through the Moran normalisation `W`.

---

## What is open

1. ~~Replication of result 8 in the Lothian cohorts.~~ **Attempted, and it
   cannot be done there.** The CHIP-grade file holds 8 TET2 clones with three or
   more waves against Fabre's 114. Reaching conventional power needs roughly 200
   to 230 — twice Fabre, twenty-five times Lothian. A replication still matters
   more than anything else on this list; it needs a cohort that does not yet
   appear to exist in the open. ([`STAGE_I_LOTHIAN.md`](analysis/STAGE_I_LOTHIAN.md))
2. ~~Fitting the model to trajectories rather than to the spectrum.~~ **Done,
   and it does not work** — not for want of data, but because the trajectory
   estimator is non-monotonic in `s` and therefore has no inverse.
   ([`STAGE_J_FABRE_GRID.md`](analysis/STAGE_J_FABRE_GRID.md))
3. **Phase.** Everything about multihit waits on it, and no amount of bulk
   sequencing supplies it.

---

## Reproducing this

Third-party data is not versioned here. See
[`reference/README.md`](reference/README.md) for provenance and licences.

```bash
uv venv && uv pip install pandas matplotlib scipy

# analyses in order; 01-05 need the Watson table, 10-12 need the Fabre table
for f in analysis/0*.py analysis/1*.py; do .venv/bin/python "$f"; done

# the parameter sweep, which is the only part that wants a cluster
nextflow run pipeline/sweep.nf -profile slurm
```
