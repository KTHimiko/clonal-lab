# Stage H — the TET2 age question, asked of data that can answer it

Reproduce with `.venv/bin/python analysis/12_tet2_longitudinal.py`.

Fabre report TET2 clones growing faster in older people. Four attempts came back
empty, each for a different reason — pooled cohorts hid 87% composition, 36
variants gave a ceiling of 0.79, ramp and wild-type decline are the same model in
relative fitness, and multihit needs phase that bulk VAF cannot resolve.

This one is different because **each clone is its own control**. A clone measured
four or five times over a median of 13.4 years can be split into an early half
and a late half and compared with itself, which removes person, cohort, protocol
and clone identity in one step.

406 clones qualify: **DNMT3A 120, TET2 114**, ages running 68 to 81 at the median.

---

## 1. Every gene decelerates, as stage F predicted

| gene | n | early rate | late rate | deceleration | 95% interval |
|---|---|---|---|---|---|
| DNMT3A | 120 | 0.0509 | 0.0096 | +0.0576 | +0.0334 to +0.0920 |
| TET2 | 114 | 0.0711 | 0.0401 | +0.0514 | +0.0152 to +0.0952 |
| other | 172 | 0.0579 | 0.0172 | +0.0445 | +0.0176 to +0.0780 |

Every interval excludes zero. Clonal interference is not a modelling artifact —
it is visible in every gene group, in the same people, within each clone's own
follow-up.

---

## 2. The comparison Fabre make, in absolute units: nothing

TET2 deceleration minus DNMT3A deceleration: **−0.0062** per year, interval
−0.0525 to +0.0389, P(TET2 decelerates less) = 0.62.

**And that null carries no information.** The spread of per-clone deceleration is
SD 0.164, so with n ≈ 114 nothing smaller than about 0.043 could have been
separated from zero. The observed 0.006 is far inside that. This measure could
not have answered the question however the biology behaved.

---

## 3. The measure the absolute difference was hiding

The two genes do not start from the same place, so comparing their deceleration
in absolute units compares them as if they did.

| gene | early rate | late rate | late 95% interval | retained |
|---|---|---|---|---|
| DNMT3A | 0.0509 | 0.0096 | −0.0174 to +0.0284 | **19%** |
| TET2 | 0.0711 | 0.0401 | **+0.0133 to +0.0598** | **56%** |

**Past a median age of 75, DNMT3A clones are no longer measurably growing.
TET2 clones still are.**

That is the sharpest statement this data supports, and it is about each gene on
its own — a single-gene interval does not depend on a between-group test
surviving.

### The between-gene comparison is borderline, and stays borderline

| half | TET2 − DNMT3A | 95% interval | P(TET2 higher) |
|---|---|---|---|
| early | +0.0202 | −0.0163 to +0.0547 | 0.844 |
| late | +0.0305 | −0.0017 to +0.0644 | **0.968** |

The gap widens in the older half, which is the predicted direction. But the
one-sided P clears 0.95 while the two-sided interval does not clear zero — 3.2%
of the bootstrap sits below it. **Calling this separated would be choosing the
statistic that says so.**

### Three measures were tested and this is the one that worked

| measure | result |
|---|---|
| deceleration, absolute | nothing (−0.0062, P = 0.62) |
| early-half rate | nothing (+0.0202, P = 0.84) |
| late-half rate | borderline (+0.0305, P = 0.97) |

Reporting the third without the first two is the standard way to manufacture a
finding. Three related measures on one dataset, and the best lands at the edge —
roughly what chance produces from time to time.

What keeps it from being noise-mining: **the direction was specified in advance**,
by Fabre and by Mon Père, and all three measures point the same way. A prediction
made before the test is worth more than a p-value found after it. Neither is
worth as much as a replication, and there is not one here.

---

## 4. Controlling for what actually drives deceleration

Stage F found deceleration is a property of the host — the carrier's burden
growth gave R² = 0.956 against 0.003 for the clone's own fitness. TET2 carriers
do differ: burden growth 0.0372 against 0.0226 for DNMT3A carriers.

A regression of deceleration on gene, starting size, carrier burden growth and
starting age gives standardised coefficients:

| term | coefficient |
|---|---|
| carrier burden growth | **+0.117** |
| age at start | −0.069 |
| log VAF at start | −0.048 |
| TET2 vs DNMT3A | −0.032 *(interval −0.167 to +0.092)* |
| **R²** | **0.019** |

Burden growth is again the largest term, and gene the smallest. But R² = 0.019
says almost none of the per-clone variation is explained by any of it — the
per-clone deceleration estimates are dominated by measurement noise at this depth
and number of timepoints, which is the same wall stage G hit.

---

## What this settles

**TET2 clones are still growing after 75 and DNMT3A clones are not**, in 234
trajectories, each clone serving as its own control. That is consistent with
Fabre's claim and is the first positive result in this arc.

**It does not establish that TET2 fitness rises with age.** A clone that is still
growing when others have stopped is equally consistent with the multihit route
(a second driver arriving late), with TET2 being protected from whatever fills
the niche, or with TET2 simply being fitter and therefore further from its
ceiling. Distinguishing those still needs phase.

**And the between-gene test is at the edge.** Treat it as consistent with the
hypothesis rather than as evidence for it.

---

## Limitations

- **One cohort**, SardiNIA, entering at 55. No replication.
- **Three measures tested**, one reached the edge. Flagged rather than hidden.
- **R² = 0.019** on the per-clone model: the estimates are noise-dominated, and
  only the aggregate medians are informative.
- **Splitting a trajectory halves its leverage**, so each half's rate is much
  noisier than the full-trajectory rate. That is the price of using each clone as
  its own control, and it is why nothing per-clone is reportable.
- **Survivorship.** Clones and carriers that dropped out between phases are
  absent, and a clone that vanished is one that stopped growing.
