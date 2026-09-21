# clonal-lab

Agent-based modelling of clonal competition in ageing tissue, calibrated
against public data on clonal haematopoiesis.

As we age, somatic mutations let some blood stem cells outcompete their
neighbours. The resulting clonal expansions — **clonal haematopoiesis** — are
common after 60 and roughly double cardiovascular risk. This project asks how
strong those fitness advantages are, by simulating the competition and fitting
the simulation to what is measured in people.

## Why simulation

Selection alone gives exponential growth. Drift alone gives a random walk.
Together they have no closed-form solution for the distribution of clone
sizes, so the only way to know what the model predicts is to run it — many
thousands of times.

The public data agrees: age explains only 9% of the variance in clone size
(R² = 0.089, p = 3 × 10⁻²²). The signal is real and the spread dominates.
**What must be reproduced is the cloud, not the line.**

## Status

| Stage | Goal | |
|---|---|---|
| Reference survey | inventory public code and data | **done** |
| Exploration | establish what the data can answer | **done** |
| **A** — minimal model | Moran process with selection, validated analytically | **done** |
| **B** — realistic model | multiple clones, continuous mutation influx | **done** |
| **C** — parameter sweep | ABC inference on an HPC cluster | **done** |
| **D** — per variant | fitness per variant class, with bootstrap | **done** |
| **E** — method bias | which study design recovers a known truth | **done** |

## Layout

```
model/        the simulator
analysis/     exploration, validation, figures
reference/    third-party code and data (not versioned)
```

## Running

```bash
uv venv && uv pip install pandas matplotlib scipy
.venv/bin/python analysis/01_exploration.py    # explore the public data
.venv/bin/python analysis/02_validation.py     # validate the single-clone model
.venv/bin/python analysis/03_multiclone.py     # many clones, VAF distributions
nextflow run pipeline/sweep.nf -profile slurm  # the grid, on a cluster
.venv/bin/python analysis/04_abc.py            # score the grid, get a posterior
.venv/bin/python analysis/05_per_variant.py    # fitness per variant class
.venv/bin/python analysis/06_method_bias.py    # what each study design recovers
```

## Key results so far

**The simulator is correct.** It reproduces the exact Moran fixation
probability in the neutral case and under selection, across several population
sizes and starting points — all within 3 standard errors.

**The Moran approximation is s/(1+s), not 2s.** The widely quoted "2s" rule
comes from the Wright-Fisher model, which has a different offspring variance.
Using the wrong one would double every fitness estimate.

**Exact beats approximate, and is faster here.** The multi-clone model first
used tau-leaping, which let a one-cell clone die and give birth in the same
time slice and inflated survival by 3 standard errors. Sampling the
birth-death transition exactly removed the bias and allowed a time step five
times larger.

**Fitness is identifiable from this data; the mutation rate is not.** A 180-point
grid scored by approximate Bayesian computation puts the selection coefficient
at s = 0.13 per year, with every accepted point at the same value. The mutation
rate stays spread over 91% of the grid. The distance surface shows why: a
vertical ridge running across all mutation rates. Shape responds to
fitness; counts respond to mutation rate, and counts cannot be computed without
a screening denominator the data does not carry.

**Fitness is per variant, but this dataset resolves only part of it.** Fitting
each variant class separately gives DNMT3A R882 at 0.130, TET2 at 0.130 and
other DNMT3A variants at 0.120 per year. Bootstrapping the observed variants
shows only one of the three pairwise differences survives the sample size:
TET2 above non-hotspot DNMT3A. The R882 advantage, real in the literature, is
smaller than what 52 variants can resolve.

**Two published methods measure different quantities, and it shows.** Given a
simulated cohort whose fitness is known by construction, a fit to the VAF
spectrum recovers 103% of the truth, following clones for thirteen years
recovers 43%, and regressing log(VAF) on the host's age recovers 18%. The
ordering reproduces the direction of the published disagreement: Watson used
the spectrum, Fabre followed clones. A large clone genuinely grows more slowly
than its fitness, so follow-up measures realised growth rather than fitness
from birth — and the gap widens the fitter the clone is.

**The slope of log(VAF) against age has no inverse.** Across true values from
0.08 to 0.24 it rises and then falls, so one observed slope is compatible with
several very different truths. In the real data, moving only the detection
floor — same cohorts, same variants — changes it by a factor of 4.1. Stage 1's
0.048 and stage D's 0.13 were never in conflict; they are two instruments, one
of which does not measure what its units suggest.

**The estimate was checked against the literature, and a parameter was wrong.**
Reading Watson 2020, Mitchell 2022 and Fabre 2022 showed that two independent
methods put the stem-cell population size at twice our initial guess. Re-running
with the corrected value improved the fit 2.4-fold and moved s from 0.10 to
0.13 per year — into the range of published per-variant estimates.

See [`analysis/METHOD_BIAS.md`](analysis/METHOD_BIAS.md) for the study-design
experiment, [`analysis/FINDINGS.md`](analysis/FINDINGS.md) for the exploration,
[`analysis/LITERATURE_REVIEW.md`](analysis/LITERATURE_REVIEW.md) for what the
papers changed, and [`reference/README.md`](reference/README.md) for data
provenance.

## Data

All data used is open or CC0. The main calibration set is Watson et al. 2020
(*Science*), released under CC0 via Dryad: 1,674 driver variants with allele
frequency, age, gene and cohort, of which 1,012 carry a usable age.

No controlled-access data is used.
