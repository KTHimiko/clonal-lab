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

## The diagnosis — corrected after the fact, twice

**The first version of this document got the diagnosis wrong, and the grid it
was committed alongside contradicts it.** The wrong version is quoted below,
because deleting it would hide how the error happened.

> *"The model's clones are about twice as large as the observed ones […] These
> are not two failures. They are one."*

### What is actually true about clone sizes

At the point that fits the spectrum best — mean 0.06, shape 4.0, **N = 200,000**
— the VAF distribution matches well at every quantile:

| quantile | observed | simulated | ratio |
|---|---|---|---|
| q10 | 0.0029 | 0.0026 | 0.91 |
| q25 | 0.0045 | 0.0035 | 0.78 |
| q50 | 0.0072 | 0.0061 | 0.85 |
| q75 | 0.0147 | 0.0149 | 1.01 |
| q90 | 0.0401 | 0.0386 | 0.96 |

**The model reproduces the clone size distribution.** The "twice as large" claim
came from a spot check made before the sweep, with N pinned at 100,000. Letting N
move fixed it, and the write-up was not updated to notice.

### What is actually wrong

The failure is at the bottom of the growth-rate distribution, and only there:

| quantile | observed | simulated | difference |
|---|---|---|---|
| q10 | −0.0451 | +0.0187 | **+0.0638** |
| q25 | −0.0002 | +0.0551 | +0.0553 |
| q50 | +0.0531 | +0.0857 | +0.0326 |
| q75 | +0.1006 | +0.1063 | +0.0057 |
| q90 | +0.1548 | +0.1259 | −0.0289 |

The top of the distribution matches. **The model does not produce the clones
that shrink:** 5.6% against 25.1% observed.

### Testing the detection floor, which was the obvious suspect

At a floor of VAF 0.002 the measurement noise on a fitted slope is 0.030 per
year, which is large enough to manufacture apparent decline. Raising the floor
on the observed data — which removes the noisiest clones — tests that directly:

| floor | n | % declining | width | expected noise |
|---|---|---|---|---|
| 0.002 | 518 | 25.1% | 0.1999 | 0.0302 /yr |
| 0.008 | 238 | 23.5% | 0.1387 | 0.0195 /yr |
| **0.030** | 74 | **18.9%** | **0.1138** | 0.0102 /yr |

**Two answers, and they point different ways.**

The *width* discrepancy was largely measurement noise. At a floor of 0.03 the
observed width is 0.1138 against the model's 0.1072 — they nearly agree, and the
gap that motivated this whole stage was mostly an artifact of fitting slopes to
clones measured at three mutant reads.

The *declining fraction* is not. It falls from 25.1% to 18.9% and stops there,
with the noise three times smaller. **About a fifth of real clones are genuinely
shrinking**, and the model produces a quarter of that.

### Why the model cannot do it

It can — at the wrong size. The grid reaches 42.8% declining clones, but only at
mean 0.16 with N = 400,000, where the median VAF is **0.4127**: fifty-seven times
the observed 0.0072.

**In this model, decline is a consequence of size.** A clone shrinks when fitter
neighbours overtake it, which requires the marrow to be filling. A small clone in
a mostly empty marrow has nothing to lose to. Reality has small clones that
shrink anyway, and the model has no mechanism for that.

That is the structural inadequacy, stated correctly: **the model couples decline
to clone size, and the data does not.**

---

## What would decouple them

None of these was tested, and each is a different model rather than a different
parameter:

- **Fluctuating fitness.** A clone's advantage is fixed for life here. If it
  varies — with inflammation, infection, treatment — clones would decline at any
  size.
- **Blood composition.** VAF measures a clone's share of circulating cells, not
  of stem cells. Shifts in lineage output move VAF without moving the clone.
- **Niche structure.** Competition is global here. If clones sit in
  sub-compartments, a small clone can lose its local contest while the marrow at
  large is empty.

---

## Limitations

- **`mu` fixed at 2×10⁻⁶.** With N varying eightfold the clone influx `N·mu`
  varies with it, so N and `mu` are not separately identified. Only N's effect on
  VAF *scale* was tested.
- **Two replicates per point**, 6,000 people each.
- **51 replicate triplets** is a small basis for the overdispersion estimate,
  though the direction is unambiguous.
- **One cohort.** SardiNIA, entering at 55.
