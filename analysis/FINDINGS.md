# Exploration findings — Watson et al. 2020

Reproduce with `.venv/bin/python analysis/01_exploration.py`.
Figures in `analysis/figures/`.

---

## 1. The data loss is structural, not random

663 of 1,674 rows have no age. But the loss is not diffuse:

| Cohort | Lost | Total | |
|---|---|---|---|
| Jaiswal2014 | 466 | 466 | **100%** |
| Genovese2014 | 197 | 197 | **100%** |

**The two largest cohorts were excluded entirely.** No other cohort lost a
single row.

This changes the interpretation: we are not working with a 40% smaller sample
of the same universe — we are working with a different universe, missing the
two classic population studies of CHIP. **1,012 variants** remain, and any
conclusion applies to the seven surviving cohorts, not to the set the paper
analysed.

> If Jaiswal and Genovese are ever needed, the ages are in the original
> papers — just not in this aggregated table.

---

## 2. The detection limit varies 300-fold across cohorts

| Cohort | n | smallest VAF | median VAF | median age |
|---|---|---|---|---|
| Young2019 | 158 | 0.0008 | 0.0019 | 61 |
| Young2016 | 26 | 0.0011 | 0.0034 | 58 |
| Acuna2017 | 182 | 0.0016 | 0.0066 | 57 |
| McKerrel2015 | 112 | 0.0079 | 0.0272 | 78 |
| Desai2018 | 31 | 0.0182 | 0.0344 | 70 |
| Coombs2017 | 378 | 0.0192 | 0.0509 | 70 |
| **ZinkWGS** | 125 | **0.2400** | 0.3300 | 71 |

ZinkWGS is whole-genome sequencing — shallow depth, so it only sees clones
that already occupy a quarter of the blood. Young2019 uses deep targeted
sequencing and sees clones 300 times smaller.

**Direct consequence for the model:** comparing a simulated distribution with
an observed one requires applying each cohort's own truncation to the
simulation. Without it, the model is penalised for predicting small clones the
instrument could never have seen.

---

## 3. The growth signal is real, and noisy in the right way

Regression of log(VAF) on age, all 1,012 variants:

```
slope      +0.0387 per year   (standard error 0.0039)
p-value     3.4 x 10^-22
R-squared   0.089
```

Highly significant with almost no explanatory power. This is **not** a
contradiction — it is the signature of a stochastic process.

Age does not determine clone size: it determines the distribution of possible
sizes. Each variant has its own fitness, and drift dominates while the clone
is small. A deterministic model would fit the line and miss everything that
matters; **what we need to reproduce is the cloud, not the line.**

This is the strongest argument for agent-based simulation over a closed-form
expression.

---

## 4. Each driver gene grows at its own rate

| Gene | n | slope (per year) | R-squared | |
|---|---|---|---|---|
| TET2 | 95 | **+0.0969 +/- 0.0141** | 0.338 | significant |
| DNMT3A | 421 | **+0.0468 +/- 0.0063** | 0.116 | significant |
| JAK2 | 45 | +0.0195 +/- 0.0110 | 0.068 | not significant |

TET2 grows roughly **twice as fast** as DNMT3A in this sample, with three
times the R-squared — the age-size relation is far cleaner for it.

JAK2 does not reach significance, but with n=45 that says more about sample
size than about biology.

This sets a model requirement: **fitness is per variant, not global.** A
simulator with a single `s` cannot reproduce three genes growing at different
rates.

---

## 5. R882 confers a measurable advantage

The best-known CHIP hotspot, against the other variants of the same gene:

| Group | n | median VAF | median age |
|---|---|---|---|
| DNMT3A R882 | 96 | **0.0390** | 63 |
| DNMT3A other | 325 | 0.0215 | 63 |

Mann-Whitney U: **p = 5.4 x 10^-5**

**Same median age, clones nearly twice the size.** Because age is controlled
by construction, the difference is attributable to the variant rather than to
exposure time.

This is an excellent validation target: a well-calibrated model should infer a
larger `s` for R882 than for the rest of DNMT3A, and the ratio between them
should be comparable to what is observed here.

---

## 6. What this data does not allow

**Population prevalence of CHIP.** The file carries only detected variants,
with no denominator of how many people were screened per age band. Prevalence
requires that denominator, and it is neither present nor inferable, because
each cohort screened populations of different sizes.

**Individual trajectories.** Each row is a single measurement of one person.
Temporal dynamics are inferred by comparing people of different ages, not by
following the same ones. This is a cross-sectional study and carries all the
limitations of one: cohort effects, survivorship bias, and the impossibility
of observing clones that vanished.

---

## What this decides about the model

| Finding | Requirement it imposes |
|---|---|
| Spread dominates (R-squared 0.09) | stochastic model, not deterministic |
| Detection limit varies 300-fold | truncate the simulation per cohort before comparing |
| Genes grow at different rates | fitness per variant, not one global `s` |
| R882 clones nearly twice the size | quantitative validation target |
| No population denominator | calibrate on the VAF distribution, never on prevalence |
| Cross-sectional design | compare simulated populations by age, not trajectories |

**The calibration observable, stated precisely:** the distribution of VAF
among detected clones, stratified by age and by gene, truncated at the
detection limit of the corresponding cohort.
