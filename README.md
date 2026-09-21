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
| D — confrontation | compare against observed VAF distributions | planned |

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
at s = 0.10 per year, with every accepted point at the same value. The mutation
rate stays spread over 91% of the grid. The distance surface shows why: a
vertical ridge at s = 0.10 running across all mutation rates. Shape responds to
fitness; counts respond to mutation rate, and counts cannot be computed without
a screening denominator the data does not carry.

**Fitness is per variant, not global.** In the public data TET2 clones grow
about twice as fast as DNMT3A ones, and within DNMT3A the R882 hotspot
produces clones nearly twice the size at the same median age
(p = 5.4 × 10⁻⁵).

See [`analysis/FINDINGS.md`](analysis/FINDINGS.md) for the full exploration
and [`reference/README.md`](reference/README.md) for data provenance.

## Data

All data used is open or CC0. The main calibration set is Watson et al. 2020
(*Science*), released under CC0 via Dryad: 1,674 driver variants with allele
frequency, age, gene and cohort, of which 1,012 carry a usable age.

No controlled-access data is used.
