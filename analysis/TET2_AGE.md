# The last open problem — TET2 and age-dependent fitness

Reproduce with `.venv/bin/python analysis/08_tet2_age.py`.
Figure `13_tet2_age.png`.

Fabre report that "age was a significant factor specifically for TET2-mutant
clones, which grew faster in older individuals." Our model has no
age-dependent fitness, and its competition term makes every clone slow down
with age, so this runs opposite to everything the simulator does.

---

## 1. The mechanism is not the obvious one

Mouse work attributes the effect mainly to **"the aging-associated reduction in
fitness of aged competitor non-mutant HSCs"** — the TET2 clone does not speed
up, the wild type slows down. In a Moran process fitness is relative, so those
are two different models:

| | what changes | who it affects |
|---|---|---|
| `s_ramp` | the mutation itself gets better with host age | only the classes given a ramp |
| `wt_decline` | the wild type gets worse with host age | **every clone at once** |

Both are now in the model. Median detected VAF at age 70, as a ratio to a world
with no age effect:

| world | DNMT3A-like class | TET2-like class |
|---|---|---|
| TET2 fitness ramps 3%/yr | **×0.96** | ×1.28 |
| wild type decays 1%/yr | **×2.04** | ×2.10 |

**A driver-specific ramp leaves the other gene alone. A decaying wild type drags
it along, and by almost exactly the same factor.** That is a clean
discriminator, and it needs no new sequencing: it asks whether the age effect is
specific to TET2 or shared with every driver.

It also resolves a tension I had not noticed. The mouse mechanism predicts a
**global** acceleration, not a TET2-specific one. The two findings only fit
together if TET2 is differentially *protected* from whatever degrades the wild
type — which is what the inflammation literature reports, TET2-deficient
progenitors responding more to IL-1/IL-6 and being shielded from
inflammation-induced damage.

---

## 2. What the data says, three ways

| test | n TET2 | slope difference | 95% interval | P(TET2 steeper) |
|---|---|---|---|---|
| pooled, floor at VAF 0.0192 | 36 (2 cohorts) | +0.0100 | −0.018 to +0.038 | 0.757 |
| pooled, no floor | 75 (4 cohorts) | **+0.0324** | −0.005 to +0.089 | **0.948** |
| **within cohort, no floor** | 72 (3 cohorts) | **+0.0042** | −0.034 to +0.043 | **0.618** |

Dropping the detection floor nearly triples the apparent effect and takes it to
the edge of significance. **Stratifying by cohort destroys it.**

Per cohort:

| cohort | n TET2 | slope TET2 | slope DNMT3A | difference |
|---|---|---|---|---|
| Acuna2017 | 9 | −0.0132 | +0.0015 | **−0.0148** |
| Coombs2017 | 35 | +0.0217 | +0.0111 | +0.0105 |
| Young2019 | 28 | +0.0416 | +0.0392 | +0.0024 |

Not even consistently positive.

### Why pooling inflates it

The two genes are not drawn from the same cohorts in the same proportions:

| cohort | median age | share of TET2 | share of DNMT3A |
|---|---|---|---|
| Acuna2017 | 56 | 12% | **29%** |
| Coombs2017 | 71 | **47%** | 32% |
| Young2019 | 62 | **37%** | 18% |
| McKerrel2015 | 78 | **0%** | 12% |

TET2 is over-represented in the older cohorts and under-represented in the
youngest. A pooled regression then reads cohort mix as biology. It is the same
mechanism that inflated the detection-floor effect by a factor of 4.1 in stage
E, on a different question.

### The diagnostic that settles which it was

**Low power widens an interval. Confounding moves the point estimate.** Here
the estimate fell from +0.0324 to +0.0042 — an 87% collapse — while the
interval stayed a similar width. That is confounding, not noise.

---

## 3. Could this dataset detect a real ramp at all?

Two simulated worlds built to be deliberately hard to separate: constant
`s = 0.13`, against `s = 0.12` with a 3%/yr ramp, calibrated to produce the
**same median detected VAF at age 70** (agreeing to within 1.6%). Anything that
looks only at age 70 is blind to the difference by construction, so this
measures purely what the age profile adds.

| TET2 variants | P(a real 3%/yr ramp is called correctly) |
|---|---|
| **36** (after the floor) | **0.79** |
| 75 (no floor) | ~0.90 |
| 100 | 0.93 |
| 300 | 0.99 |
| 1000 | 1.00 |

**At n = 36 the ceiling is 0.79, and the real data achieved 0.757.** The data is
performing at very nearly the limit its size allows. The failure to separate is
not a failure of the analysis.

The within-cohort test costs power on top of that: each cohort spans a narrower
age range than the pooled set, so the age leverage per comparison is smaller.
Its P = 0.618 therefore does not prove absence. What it does establish is that
the pooled signal was not robust to the one confounder we already knew about.

---

## 4. What this settles

**For the project.** Stage D's TET2 estimate of 0.130 does not need an
age-dependence correction on the strength of this dataset. The effect Fabre
report may well be real — they had longitudinal data and 697 tracked clones —
but it is not visible here once cohort composition is removed, and the
simulation says it would not be visible even if present.

**For the model.** Age-dependent fitness is implemented in both forms and
validated to behave differently, so it is available when there is data that can
constrain it. Turning it on now would be fitting a parameter the data cannot
identify — the same trap as the mutation rate in stage C.

**What cannot be settled by any amount of this kind of data.** In a
relative-fitness model, "TET2 gets better" and "everyone else gets worse while
TET2 is protected" are the same model. Fitness is a ratio, and the two
descriptions are arithmetically identical. Distinguishing them needs something
the model has no access to: the **absolute** output of wild-type haematopoiesis
with age. That is a measurement, not an inference, and it is exactly why the
mouse transplantation experiments were done — they can hold one side fixed and
vary the other, which an observational human cohort cannot.

The mouse work and Fabre are therefore not in conflict. They describe the same
relative change from opposite sides.

---

## Limitations

- **Three cohorts, 72 TET2 variants.** Acuna2017 contributes 9, and its negative
  difference carries little weight but is not nothing.
- **The estimator is the age slope**, which stage E showed cannot recover `s`.
  It is used here comparatively — two genes measured with the same bent ruler —
  which is valid for ranking and not for magnitude.
- **The ramp shape is a guess.** A 3%/yr linear gain from age 50 is one
  functional form among many; a threshold effect at some inflammatory tipping
  point would leave a different signature and has not been tested.
- **No mortality.** If TET2 carriers with large clones die earlier, the oldest
  age bands are depleted of exactly the observations the test depends on.
