# Stage J — a grid built for the cohort, and an arm that cannot be inverted

Reproduce with `.venv/bin/python analysis/14_fabre_grid.py`.
Decision rule fixed in advance:
[`PREREGISTRATION_STAGE_J.md`](PREREGISTRATION_STAGE_J.md), commit `f6af384`.
Figure `16_fabre_grid.png`.

**252 jobs, 3.5 CPU hours, 36 minutes on the three-node cluster.**

---

## The verdict, honoured

**R = 0.40**, against a band of 1.8 to 3.2 committed before the sweep ran.
Outside. Per the pre-registration: **conclusion downgraded, stop.**

| arm | s | 95% interval | distance |
|---|---|---|---|
| VAF spectrum | 0.080 | 0.080 – 0.100 | **0.232** |
| trajectories | 0.200 | 0.060 – 0.220 | 0.766 |

---

## The pre-registration's own reading of this branch was wrong

It said a missed band would mean stage F's 2.42 was partly coincidence. **The
grid says otherwise, and the reason is structural.**

| `s` | median rate, trajectories | median VAF, spectrum |
|---|---|---|
| 0.04 | 0.0152 | 0.0035 |
| 0.08 | 0.0696 | 0.0064 |
| **0.12** | **0.0919** ← peak | 0.0164 |
| 0.18 | 0.0550 | 0.0587 |
| 0.24 | 0.0141 | 0.1326 |
| 0.30 | 0.0032 | 0.2005 |

**The trajectory arm is non-monotonic in `s`. The spectrum arm is monotonic
across the whole grid.**

A fitter clone is closer to saturation when follow-up begins, so past a certain
fitness the measured growth rate *falls* as `s` rises. The map folds.

The observed median of 0.0531 per year therefore matches the grid at **s ≈ 0.066
on the rising branch and s ≈ 0.182 on the falling one**. The fit took the high
branch for trajectories and the low branch for the spectrum, and R is a ratio
between two different branches of a two-valued map.

**R was never a well-defined quantity.** No value of it could have tested what
the pre-registration thought it was testing.

The verdict stands — the band was missed, the result is downgraded, we stop —
but the reason is a flaw in this stage's design, not evidence against stage F.
**Stage F compared a spectrum-fit estimate against a directly measured growth
rate. It never inverted the trajectory distribution, so it does not inherit
this defect.**

---

## What the stage did establish

**A longitudinal growth rate cannot be used as an inverse problem.** This is
stage E's finding about the age regression, reappearing in a different
estimator. It is fine to *measure* a growth rate and report it — stages F and H
do exactly that, and the 35% model-free result depends on nothing else. It is
not a quantity you can fit a model to and read a fitness out of, unless you
already know which branch you are on.

That is now true of two of the three designs this project compared. The spectrum
fit is the only one of the three that is invertible, which is a sharper reason to
prefer it than "it won our test".

**The grid itself is sound and reusable.** 84 points, both summary vectors, built
to Fabre's ages, floor and depth. The spectrum arm fits it well (distance 0.232,
interval 0.080–0.100) and that estimate no longer depends on a grid borrowed from
the Watson cohort — which was the stated limitation this stage set out to close.
It closes, for that arm.

---

## The second pre-registration with a hole in it

Stage I's plan fixed the estimator and forgot to name which file to read. This
one fixed the band and forgot to ask whether the quantity being banded was
invertible.

Both were caught — the first by a sanity check committed alongside the test, the
second by a diagnostic table this write-up demanded. Neither was caught by the
plan itself.

The lesson is narrower than "pre-registration does not work". It is that **a
pre-registered ratio needs a pre-registered check that both of its terms are
identifiable.** That check costs one column of a table and would have prevented
a 36-minute cluster run.

---

## Limitations

- **Both arms come from one model.** This tests the inversion, not the biology —
  the limit stage E stated and stage I's power calculation reinforced. External
  support for the spectrum estimator remains the phylogenetic agreement in stage
  E2, untouched by anything here.
- **Sampling mismatch.** The simulation enrols everyone at 65 and takes five
  evenly spaced draws to 78; the real cohort enrols across a range of ages with
  two to five draws. Same estimator, different sampling.
- **The mutation rate is still unidentifiable**, as stage C established. It is on
  the grid because the simulation needs a value.
