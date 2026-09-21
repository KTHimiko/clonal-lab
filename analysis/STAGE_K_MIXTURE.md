# Stage K — a fitness distribution, N free, and one failure instead of two

Reproduce with `.venv/bin/python analysis/15_mixture_fit.py`.
Decision rule and identifiability checks fixed in advance:
[`PREREGISTRATION_STAGE_K.md`](PREREGISTRATION_STAGE_K.md), commit `42e75b3`.
Figure `17_mixture_fit.png`.

**224 points in 32 batched jobs.** Stage J spent 36 minutes on 3.5 CPU-hours
running 252 jobs of 7.5 seconds; batching a whole row of means into each job cut
that to minutes.

---

## The identifiability checks came first, as the plan required

Stage J banded a ratio without asking whether its denominator could be inverted.
It could not, and no outcome could have meant anything. This plan made the check
part of the procedure.

| check | result |
|---|---|
| spectrum arm monotonic in the mean, at fixed shape and N | **16 of 16 slices** |
| joint minimum unique | **yes** — 1 point within 5% |

**The test was well-posed.** Whatever it returns, it returns about the model.

---

## The verdict

**0 of 112 grid points fit both arms.** Within 25% of the best distance on each
arm simultaneously: none.

| | mean | shape | N | d_cs | d_tr |
|---|---|---|---|---|---|
| best spectrum | 0.06 | 4.0 | 200,000 | **0.1460** | 0.0426 |
| best trajectory | 0.14 | 4.0 | 400,000 | 1.7881 | **0.0114** |
| best joint | 0.04 | 2.0 | 50,000 | 0.1787 | 0.0291 |

The point that fits the growth rates best is off by **1.79** on the spectrum —
catastrophic. The point that fits the spectrum best is 3.7× the achievable
distance on rates.

Per the pre-registration: **structural inadequacy, named rather than tuned
around. Stop.**

---

## What N did

The first time this project has let N move. It was taken from the literature in
stage C and carried unchanged through eight stages.

| N | best d_cs | best d_tr |
|---|---|---|
| 50,000 | 0.1787 | 0.0232 |
| 100,000 | 0.1672 | 0.0205 |
| 200,000 | **0.1460** | 0.0129 |
| 400,000 | 0.1473 | **0.0114** |

**N helps each arm separately and resolves nothing.** Larger N improves both
best-in-arm distances, and the tension between the arms survives intact — it
just moves to a different place on the grid.

---

## The diagnosis, and it overturns the hypothesis that motivated the stage

The stage was launched on the idea that the missing ingredient was **fitness
heterogeneity**: real clones differ from each other, ours all shared one `s`, and
that was why the observed spread of growth rates (0.200) exceeded anything a
single-`s` model could produce (0.132).

Two measurements say that is not the main story.

**First: there is no excess technical noise.** Fabre's data includes 51 usable
triplets — the same variant, in the same person, at the same draw, sequenced
three times. Comparing their scatter against binomial read sampling:

| | log(VAF) scatter |
|---|---|
| expected from read depth alone | 0.349 |
| observed across replicates | 0.319 |
| **overdispersion** | **0.8× the variance** |

Slightly *less* than binomial, not more. The assay is well behaved, and the
hypothesis that unmodelled technical noise inflated the rate spread is dead.

**Second: measurement noise is a function of clone size, and that ties the two
failures together.** Log-scale noise goes as `1/√(VAF·depth)`, so small clones
are measured imprecisely and large ones precisely:

| VAF | noise on the fitted rate | width it adds |
|---|---|---|
| **0.0072** *(observed median)* | 0.0338 /yr | **0.0867** |
| 0.0149 *(best fit's median)* | 0.0234 /yr | 0.0601 |
| 0.0600 | 0.0114 /yr | 0.0292 |

The model's clones are about **twice as large** as the observed ones. Because
they are larger, their trajectories are measured more precisely, so their fitted
rates cluster more tightly — and the rate distribution comes out too narrow.

**These are not two failures. They are one.** The model cannot produce a
population of clones as small as those observed, and everything else follows.
Roughly 40% of the observed rate spread is measurement noise on small clones,
not biological variation in fitness.

---

## What that leaves

The question is no longer "what is missing from the fitness distribution". It is
**why the model's clones are too big at every parameter combination tried**,
across an eightfold range of N.

Candidates, none tested here:

- **The detection floor.** We admit simulated clones at VAF ≥ 0.002 uniformly.
  The real assay's sensitivity varies with depth per site and per sample, and
  Fabre applied a filtering step we did not reproduce.
- **When clones start.** The model seeds mutations at a constant rate from birth.
  If real driver acquisition is weighted later in life, clones would be younger
  and smaller at the same age.
- **The niche.** The Moran normalisation lets a clone grow until it meets
  competition. Structure in the marrow — clones confined to sub-compartments —
  would cap them earlier.

Each is a different model, not a different parameter. That is what "structural"
means here, and it is why the pre-registration's answer is to stop rather than
keep fitting.

---

## Limitations

- **`mu` fixed at 2×10⁻⁶.** With N varying eightfold the clone influx `N·mu`
  varies with it, so N and `mu` are not separately identified. Only N's effect on
  VAF *scale* was tested.
- **Two replicates per point**, 6,000 people each.
- **51 replicate triplets** is a small basis for the overdispersion estimate,
  though the direction is unambiguous.
- **One cohort.** SardiNIA, entering at 55.
