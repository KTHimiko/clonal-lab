# Pre-registration — stage O, the joint sweep

**Committed before the sweep ran.** Fourth of these, and it carries what the
previous three earned: stage I's plan forgot to name a file, stage J's banded a
quantity that had no inverse, stage K's checked identifiability first and was
clean.

---

## What is being tested

Stage N rejected four candidate mechanisms and bounded what was left. A hazard
tied to the **carrier's age** rather than the clone's was the one variation the
bound did not exclude, because it stops preferentially removing old — and
therefore large — clones.

Tested at one parameter point it came closer than anything before it:

| | d spectrum | declining | d rates | ICC | persistence |
|---|---|---|---|---|---|
| no mechanism | 0.139 | 4.4% | 0.0445 | −0.02 | 0.05 / 0.09 |
| **hazard 0.03 from age 55** | **0.177** | **22.0%** | **0.0176** | **0.108** | 0.37 / 0.19 |
| observed | — | 25.1% | — | 0.129 | 0.51 / 0.36 |

All four statistics land in the right region at once, for the first time. The
spectrum misses by 1.27×, against 2.6× to 6× for every mechanism before it.

**And the test had one hand tied.** Every run fixed the clone population at
mean 0.06, shape 4.0, N 200,000 — values fitted with **no mechanism at all**.
With a hazard removing clones, the best-fitting population is a different one.

This sweep frees them together.

---

## The identifiability check, run before any fit is reported

Stage J banded a ratio without asking whether its denominator could be inverted.
That check is part of the plan now:

1. **Is the spectrum arm monotonic in the mean**, at fixed shape, N, rate and
   onset? It was in stage K; adding a hazard could fold it.
2. **Is the accepted region connected**, or does it split? A split region is
   reported as such rather than resolved by picking the lowest point.

---

## The decision rule

Distances are RMS differences between simulated and observed quantiles — log
space for VAF, raw for rates. `best_cs` and `best_tr` are the smallest each arm
achieves anywhere on the grid.

| outcome | action |
|---|---|
| **some point is within 25% of `best_cs` AND within 25% of `best_tr`** | the mechanism fits both arms once the population is free. Report it, then check the three statistics it was **not** fitted to — declining fraction, ICC, persistence — and say plainly whether they follow. **Stop.** |
| **no point is** | four mechanisms and a joint search have failed. The bracket from stage N stands as the result, and the missing mechanism is not any of them. **Stop.** |

Either way, stop. Three previous pre-registrations have all needed that line.

**What a pass would and would not mean.** Fitting the spectrum and the rates is
the fitted part. The declining fraction, the within-person correlation and the
persistence are consequences the sweep does not optimise for — they are the test.
Stage M's mechanism fitted everything it was aimed at and died on exactly such a
statistic.

---

## The grid

| | |
|---|---|
| mean of the fitness gamma | 0.05, 0.07, 0.09, 0.11 |
| shape | 1.0, 2.0, 4.0 |
| N | 100,000 · 200,000 · 400,000 |
| hazard rate | 0 · 0.02 · 0.03 · 0.05 per clone-year |
| hazard onset (carrier age) | 55 · 60 |
| drop | fixed at 0.20 |
| replicates | 2 |
| points | **576**, batched into 144 jobs of four means each |

`switch_delta` is fixed because stage N found the rate, not the depth of the
drop, is what moves the outcome — 0.15 and 0.25 gave nearly the same picture.
That is an assumption, and it is stated rather than hidden.

Rows with rate 0 ignore the onset; the duplicates are kept so every cell has the
same replicate structure, and are collapsed in the analysis.
