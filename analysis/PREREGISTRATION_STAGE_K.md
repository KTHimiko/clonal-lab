# Pre-registration — stage K, the mixture grid with N free

**Committed before the sweep ran.** Third attempt at writing one of these, and it
carries the correction the second one earned.

---

## What stage J found, and what it left

The trajectory arm could not be fitted, and the diagnosis was not the estimator.
Observed per-clone growth rates span **0.200** from the 10th to the 90th
percentile. A single-`s` model produces at most **0.132** anywhere on its grid.
Real clones differ from each other in fitness; ours did not.

A spot check showed a gamma mixture reaches **0.2015** — essentially exact — and
reproduces the 25% of clones that are shrinking, which a single `s` never did.

**But every mixture wide enough to do that made clones far too large in the VAF
spectrum**, by a factor of 8 to 24. The two observables pull in opposite
directions, monotonically across nine spot-checked points.

`N` was taken from the literature in stage C and never fitted. More stem cells
means the same clone is a smaller share of the blood, which lowers VAF without
touching growth rates — the exact direction needed. So `N` joins the grid.

---

## The question, stated so it can fail

**Is there any (mean, shape, N) at which the model fits both observables of one
cohort at once?**

Not "what is the best fit" — a best fit always exists. Whether a region exists
where neither arm is badly wrong.

---

## The identifiability check this plan owes stage J

Stage J banded a ratio without asking whether its denominator could be inverted.
It could not, so no outcome could have meant anything. **That check is now part
of the plan rather than a diagnosis after the fact.**

Before any fit is reported, two things get printed:

1. **Is the spectrum arm monotonic in `mean` at fixed (shape, N)?** If it folds,
   the same defect applies and the fit is not reported as an estimate.
2. **Is the joint minimum unique?** Counted as local minima on the joint
   distance surface. More than one is reported as such, not resolved by picking
   the lowest.

---

## The decision rule

Distances on each arm are the RMS difference between simulated and observed
quantiles — in log space for VAF, raw for rates.

Let `best_cs` and `best_tr` be the smallest distance any grid point achieves on
each arm separately.

| outcome | reading | action |
|---|---|---|
| **some point is within 25% of `best_cs` AND within 25% of `best_tr`** | the model fits both once `N` is free; report that point and its neighbourhood | **stop** |
| **no point is** | the model cannot fit both even with `N` free; a structural inadequacy, named rather than tuned around | **stop** |

The 25% band is chosen now, before seeing anything, and is deliberately
generous: the question is whether a jointly-acceptable region exists at all, not
where its optimum sits.

**Both branches stop.** As in stage J, that is the point.

---

## The grid

| | |
|---|---|
| mean of the gamma | 0.04 to 0.16, step 0.02 (7) |
| shape | 0.5, 1.0, 2.0, 4.0 (4) — lower is wider |
| **N** | **50,000 · 100,000 · 200,000 · 400,000** |
| replicates | 2 |
| people per point | 6,000 |
| total points | 224, batched into **32 jobs** of 7 means each |

Batched because stage J spent 36 minutes on 3.5 CPU-hours: 252 jobs of 7.5
seconds each, where the scheduler cost more than the work.

---

## Fixed, and therefore a limitation

`mu` stays at 2×10⁻⁶. Stage C established it is unidentifiable from clone-size
shape, and it moves counts rather than quantiles. With `N` varying eightfold the
clone influx `N·mu` varies with it, which is the physically honest choice and
also means `N` and `mu` are not separately identified here. Only `N`'s effect on
VAF *scale* is being tested.
