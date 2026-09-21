# Stage F — the design gap, measured in people

Reproduce with `.venv/bin/python analysis/10_real_design_gap.py`.
Figure `14_real_design_gap.png`.

Stage E measured how much each study design under-recovers a known `s`, and
carried one obvious objection: every number came from our own model, so the
biases are real only to the extent the model is.

Fabre's SardiNIA data (figshare `10.6084/m9.figshare.15029118`, CC BY 4.0) has
what the Watson table never had — a person identifier and repeated measurements.
**394 people, 994 clones, two to five timepoints each, median follow-up 12.9
years, median depth 1,140×.** All three designs can therefore be run on the same
individuals.

---

## 0. First, can we read their file correctly?

| gene | n | our median rate/yr | published |
|---|---|---|---|
| TET2 | 210 | 0.0826 | 0.068 |
| DNMT3A | 194 | 0.0483 | 0.062 |
| SF3B1 | 32 | 0.0946 | — |
| ASXL1 | 30 | 0.0975 | — |
| TP53 | 30 | 0.0616 | — |

Same ordering, magnitudes within about a third. **The gap is expected and worth
naming rather than hiding:** Fabre fit a Bayesian model with a fitted
overdispersion term, per-clone random effects and a treatment of the detection
limit. We fit ordinary least squares to log(VAF). The simpler estimator is the
right one here — it is exactly what stage E simulated — but these are our
numbers on their data, not a reproduction of their analysis.

---

## 1. The cleanest result in the project, and it needs no model

Same cohort, same assay, same fitting procedure. The only thing that changes is
whether a clone was already big enough to be detectable at enrolment.

| subset | n | median rate/yr | 95% interval |
|---|---|---|---|
| every clone with a trajectory | 770 | **0.0703** | 0.0611 – 0.0771 |
| only those detectable at baseline (VAF ≥ 0.0192) | 128 | **0.0455** | 0.0224 – 0.0609 |

**Conditioning on detectability costs 35% of the measured growth rate, and the
intervals do not overlap.**

This is stage E's longitudinal bias, in real people, with no simulation
anywhere in the path. A study that enrols the clones it can already see is
measuring the ones that have most nearly finished growing.

---

## 2. The three designs on one cohort

| design | estimate |
|---|---|
| age regression, first timepoint (stage 1 design) | **0.0024** /yr |
| longitudinal, enrolled clones (Fabre design) | **0.0455** /yr |
| spectrum fit, first timepoint, ages 65–75 (Watson design) | **0.110** /yr *(bootstrap 0.101–0.120)* |

Ordering predicted by stage E: **spectrum > longitudinal > age regression.**
Ordering observed: **the same.**

| ratio | predicted | observed |
|---|---|---|
| spectrum / longitudinal | 2.40 | **2.42** |
| longitudinal / age regression | 2.39 | 18.6 |
| spectrum / age regression | 5.72 | 45.0 |

**One ratio matched closely and the other two did not.** Stating that plainly
matters more than the one that did.

**Where it held.** The spectrum-to-longitudinal ratio is the comparison stage E
was built to make, and it is also the pair the published literature disagrees
about — Watson's 15.0% against Fabre's 6.2% for DNMT3A is a ratio of 2.42. Three
routes to roughly 2.4: a simulation with known truth, two published papers, and
now two estimators run on one cohort. The decimals agreeing is luck; the
magnitude arriving three times is not.

**Where it failed.** Stage E predicted the age regression would recover about a
sixth of the spectrum estimate. Here it recovers a fortieth — it is effectively
zero, with a bootstrap interval spanning zero. The direction was right and the
magnitude was badly wrong. Two reasons, neither flattering to the prediction:
this cohort enrols at 55 and runs to 105, so it sits entirely in the regime
where stage E's own curve shows the slope collapsing; and n = 139 makes the
point estimate unstable.

---

## 3. A stage E prediction, tested on data that did not produce it

Stage E found that moving only the detection floor changed the age slope by a
factor of 4.1 in the Watson data. This cohort is sequenced far deeper, so the
same lever can be pulled on data that had no part in making the prediction.

| floor | n | slope | 95% interval |
|---|---|---|---|
| none (VAF > 0) | 799 | **0.0278** | +0.0160 to +0.0399 |
| 0.005 | 452 | 0.0104 | −0.0011 to +0.0220 |
| 0.0192 | 139 | **0.0024** | −0.0117 to +0.0178 |

Same direction, and a factor of 11 across the range. At the floor our grid uses,
the estimate is **not distinguishable from no age effect at all**.

The age regression is not measuring fitness. It is measuring the instrument.

---

## What this settles

**Stage E's central claim survives contact with real data, in the comparison it
was built for.** The spectrum fit and longitudinal follow-up differ by a factor
near 2.4 in the same people, as predicted, and the mechanism — conditioning on
detectability — is demonstrable in those people without any model at all.

**Its quantitative prediction for the age regression does not survive.** The
design fails worse in this cohort than the simulation said it would.

**And the model-free result is the one to keep.** Of everything in stages E and
F, the finding that will still be true if the simulator is wrong is this:
enrolling clones by detectability drops the measured growth rate by 35%, in 770
real trajectories, with non-overlapping intervals.

---

## Limitations

- **One cohort, and an old one.** Entry at 55, median first-phase age 70, up to
  105. The regime where stage E says the designs agree — young carriers, small
  clones — is not represented at all.
- **The spectrum estimate is model-dependent.** It is scored against our grid,
  so the 0.110 inherits our N, our mutation rate range and our Moran
  parameterisation. The longitudinal estimate does not. Half of the 2.42 is
  therefore ours rather than the data's.
- **n = 75 clones** in the 65–75 age band used for the spectrum fit.
- **Our estimator is not theirs.** See section 0.
- **SardiNIA is a founder population** from Sardinia, with its own genetic
  structure. Whether clonal dynamics generalise from it is not something this
  analysis can address.
