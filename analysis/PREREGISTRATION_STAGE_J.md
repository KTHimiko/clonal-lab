# Pre-registration — stage J, the tailored sweep

**Committed before the sweep was launched.** Stage I showed what a plan is worth
when it fixes the analysis but leaves a hole, so this one is short and its
decision rule is the whole point.

---

## What is being tested

Stage F measured, in 394 real people, a ratio of **2.42** between a VAF-spectrum
fit and thirteen years of longitudinal follow-up. Stage E predicted 2.40 from
simulation, and the two published papers differ by 2.42 on the same gene.

That agreement has a stated weakness: **the numerator came from a grid built for
a different cohort.** Stage F's own limitations say "half of the 2.42 is ours
rather than the data's".

This sweep removes it. The grid is shaped to Fabre's cohort — their baseline age
distribution, their 13-year follow-up, their VAF floor of 0.002, their median
depth of 1,140× — and every grid point emits **two summary vectors from the same
simulated people**: the cross-sectional VAF spectrum, and the per-clone growth
rates a longitudinal study would measure.

Fitting the Fabre cohort both ways, against a grid built for it, is the cleanest
internal check of the project's central result that is available.

---

## The decision rule, fixed now

Let **R** be the ratio between the two estimates our own model produces when
fitted to the Fabre cohort: spectrum arm ÷ trajectory arm.

| outcome | reading | action |
|---|---|---|
| **1.8 ≤ R ≤ 3.2** | the model reproduces the design gap end to end; the borrowed-grid limitation closes | conclusion confirmed — **stop** |
| **R < 1.8 or R > 3.2** | the model does not reproduce the gap the real data shows; stage F's agreement was partly coincidence | conclusion downgraded and recorded as such — **stop** |

The band was chosen before running: wide enough not to demand precision that
three replicates cannot deliver, narrow enough that 1.2 or 4.5 means something is
wrong.

**Both branches end in stopping.** That is the actual purpose of writing this
down. The risk to this project has never been getting an answer wrong; it is
opening a fourth arc instead of closing the third.

---

## The grid

| | |
|---|---|
| `s` | 0.04 to 0.30, step 0.02 (14 values) |
| `mu` | log-spaced, 5×10⁻⁷ to 1.2×10⁻⁵ (6 values) |
| replicates | 3 |
| people per point | 8,000 |
| N | 100,000 |
| total jobs | **252** |

`mu` remains unidentifiable — stage C established that and nothing here changes
it. It is on the grid because the simulation needs a value, and the fit
marginalises over it.

---

## What will not be claimed

That the spectrum arm is correct because it agrees with itself. The grid and the
cohort simulation come from one model, so this tests the **inversion**, not the
biology — the same limit stage E stated and stage I's power calculation
reinforced. The external support for the spectrum estimator remains the
phylogenetic agreement in stage E2, and that is unchanged by this.
