# Stage O — the joint sweep, and the tension that survives it

Reproduce with `.venv/bin/python analysis/16_joint_fit.py`.
Decision rule committed before the sweep: [`PREREGISTRATION_STAGE_O.md`](PREREGISTRATION_STAGE_O.md), `0d4e080`.
Figure `18_joint_fit.png`. **576 points, 144 batched jobs.**

Every earlier run of the hazard fixed the clone population at values fitted with
no mechanism at all. This freed them together.

---

## The identifiability check, run first

**Spectrum monotonic in the mean: 63 of 63 slices.** Adding a hazard did not
fold it. Unlike stage J, whatever this returns is about the model.

---

## The verdict

**0 of 252 points** fall within 25% of the best distance on each arm.
Per the pre-registration: the bracket from stage N stands, and the missing
mechanism is none of the five tried. **Stop.**

---

## Two things the verdict hides, and both matter

### The band was anchored on a degenerate point

`best_tr` is **0.0055**, achieved at mean 0.09, N 100,000 — where `d_cs` is
**1.533** and the simulated median VAF is 0.0241 against 0.0072 observed. That
is a model catastrophically wrong about clone sizes, and requiring "within 25%
of 0.0055" is a bar only that corner can clear.

The rule was well-posed — both arms invert, the monotonicity check passed — but
**badly calibrated**, which is the same class of error as stage J one step
milder. A rule anchored on the best value of each arm assumes those bests are
attained somewhere sensible. Here one of them is not.

### The hazard improves *both* arms once the population is free

| hazard | best d_cs | best d_tr | max declining | median ICC |
|---|---|---|---|---|
| 0.00 | 0.1595 | 0.0224 | 20.7% | **+0.278** |
| 0.02 | 0.1404 | 0.0078 | 25.7% | −0.028 |
| 0.03 | 0.1514 | 0.0055 | 28.5% | −0.069 |
| 0.05 | **0.1278** | 0.0062 | 39.6% | −0.092 |

**Stage N reported the opposite** — that the hazard wrecked the spectrum. It did,
*with the population pinned*. Freed, a larger N and a higher mean absorb the
clones the hazard removes, and the spectrum comes out better than with no
mechanism at all.

That is a correction to stage N's framing, and it is the reason this sweep was
worth running even though the rule says no.

---

## What the best points actually achieve

| mean | shape | N | rate | from | d_cs | d_tr | declining | ICC | P(f\|f) | P(f\|r) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.07 | 4 | 400,000 | 0.05 | 60 | **0.135** | **0.0151** | 29.7% | **−0.040** | **0.54** | **0.29** |
| 0.07 | 4 | 400,000 | 0.03 | 55 | 0.164 | 0.0163 | 21.8% | −0.138 | 0.49 | 0.20 |
| 0.05 | 2 | 200,000 | 0.03 | 55 | 0.151 | 0.0177 | 22.5% | −0.069 | 0.44 | 0.19 |

*Observed: declining 25.1%, ICC 0.129, persistence 0.51 against 0.36.*

**The persistence is reproduced almost exactly** — 0.54 against 0.29, versus
0.51 against 0.36 — and it was never fitted. The declining fraction lands at
29.7% against 25.1%. The spectrum at 0.135 beats the no-mechanism baseline.

**Three of the four, and the fourth is the ICC.**

---

## The finding, which is the tension the ICC exposes

Read the ICC column of both tables together:

| | within-person correlation |
|---|---|
| model, no hazard | **+0.278** |
| model, any hazard tested | **negative** |
| **observed** | **+0.129** |

**The data sits exactly between them.**

The hazard fires independently per clone, so it injects uncorrelated variance
and drives the correlation down. Clonal interference, with the population free
and clones therefore larger, drives it up. Neither alone lands where the data is.

Matching the ICC therefore requires a **weak** hazard — somewhere between 0 and
0.02 — and a weak hazard produces few declining clones. **The two statistics pull
in opposite directions on the same parameter.**

That is the same shape as every constraint this project has found: a mechanism
that buys one observable at the price of another. It is now located on a single
number, which is the sharpest the tension has ever been.

---

## Where five mechanisms leave it

| mechanism | fails on |
|---|---|
| niche structure | 2 points of decline; the VAF ceiling fights the observed tail |
| fluctuating fitness | compounding variance inflates clone sizes 2.4× to 629× |
| blood composition | ICC 0.674 against 0.129 |
| hazard from birth | strips the oldest and largest clones |
| **hazard from carrier age, population free** | **ICC negative against 0.129 observed** |

The specification from stage N survives and gains a clause. The missing
mechanism must produce decline **without preferentially removing old clones**,
**without adding compounding variance**, **without acting on a whole person at
once** — and now **without destroying the modest within-person correlation that
clonal interference generates**.

The last two are nearly contradictory, which is why this is hard: decline must be
clone-specific enough not to correlate clones, yet not so clone-specific that it
washes out the correlation competition creates.

---

## Limitations

- **The rule's calibration**, discussed above. A band anchored on absolute
  distances rather than on each arm's best would have given a different verdict,
  and choosing that after seeing the data is exactly what pre-registration
  exists to prevent. The verdict stands; the flaw is recorded for the next plan.
- **`switch_delta` fixed at 0.20**, on stage N's finding that the rate dominates.
- **Two replicates**, 7,000 people per point.
- **The hazard grid's smallest non-zero rate is 0.02.** The ICC crossing sits
  below it and was not resolved.
