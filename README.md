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
| **E2** — open problems | checked against the literature and a third method | **done** |
| **E3** — TET2 and age | can this data constrain age-dependent fitness? | **done** |
| **F** — real cohort | the same three designs on 394 real people | **done** |
| **G** — multihit | can bulk VAF separate two hits in one cell from two clones? | **done** |

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
.venv/bin/python analysis/07_open_problems.py  # against a third, independent method
.venv/bin/python analysis/08_tet2_age.py       # can the data see age-dependent fitness?
.venv/bin/python analysis/09_units_check.py    # is s the same quantity everywhere?
.venv/bin/python analysis/10_real_design_gap.py # the three designs on a real cohort
.venv/bin/python analysis/11_multihit.py       # two hits in one cell, or two clones?
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
ordering reproduces the published disagreement at the right size: for DNMT3A,
Watson's spectrum fit gives 15.0% per year and Fabre's follow-up 6.2%, a ratio
of 2.42 against the 2.4 predicted from design bias alone. A large clone
genuinely grows more slowly than its fitness, so follow-up measures realised
growth rather than fitness from birth. Fabre attribute the same slowdown to "an
increasingly competitive oligoclonal landscape" — this model's competition term
in words. A third, independent method agrees: Mitchell's per-clone phylogenetic
estimates put DNMT3A clades at 0.167–0.200 per year, using coalescence patterns
rather than a frequency spectrum. Three detection-conditioned methods land
between 0.11 and 0.20; the longitudinal estimate is the outlier at 0.062.

**The slope of log(VAF) against age has no inverse.** Across true values from
0.08 to 0.24 it rises and then falls, so one observed slope is compatible with
several very different truths. In the real data, moving only the detection
floor — same cohorts, same variants — changes it by a factor of 4.1. Stage 1's
0.048 and stage D's 0.13 were never in conflict; they are two instruments, one
of which does not measure what its units suggest.

**The design gap is real, and measurable without a model.** Fabre's SardiNIA
cohort — 394 people, 994 clones, up to five timepoints each — allows all three
designs on the same individuals. Enrolling clones by detectability drops the
measured growth rate from 0.0703 to 0.0455 per year across 770 real
trajectories, with non-overlapping intervals: a 35% loss, no simulation
involved. The spectrum fit and longitudinal follow-up differ by a factor of 2.42
in these people, against 2.40 predicted from simulation and 2.42 between the two
published papers. The age regression does worse than predicted — effectively
zero, with an interval spanning zero — and moving only the detection floor
changes it by a factor of 11.

**The TET2 age effect is not visible here, and could not have been.** Fabre
report TET2 clones growing faster in older people. Pooled across cohorts this
dataset appears to agree — the gap against DNMT3A reaches P = 0.948 once the
detection floor is dropped — but stratifying by cohort collapses the estimate by
87%, because TET2 is over-represented in the older cohorts and absent from one
entirely. Simulating two worlds calibrated to be identical at age 70 shows that
36 variants could call a real 3%/yr ramp correctly only 79% of the time, and the
data achieved 0.757. **It is performing at the ceiling its size allows.** About
300 variants would be needed. A deeper point survives the arithmetic: in a
relative-fitness model, "TET2 improves" and "everyone else degrades while TET2 is
protected" are the same model, and separating them needs an absolute measurement
of wild-type output that no observational cohort provides.

**The estimate was checked against the literature, and a parameter was wrong.**
Reading Watson 2020, Mitchell 2022 and Fabre 2022 showed that two independent
methods put the stem-cell population size at twice our initial guess. Re-running
with the corrected value improved the fit 2.4-fold and moved s from 0.10 to
0.13 per year — into the range of published per-variant estimates.

See [`analysis/METHOD_BIAS.md`](analysis/METHOD_BIAS.md) for the study-design
experiment and
[`analysis/METHOD_BIAS_LITERATURE.md`](analysis/METHOD_BIAS_LITERATURE.md) and
[`analysis/OPEN_PROBLEMS.md`](analysis/OPEN_PROBLEMS.md) for what checking it
against the papers corrected,
[`analysis/TET2_AGE.md`](analysis/TET2_AGE.md) for the age-dependence question,
[`analysis/STAGE_F_REAL_COHORT.md`](analysis/STAGE_F_REAL_COHORT.md) for the
same designs run on a real cohort,
[`analysis/MULTIHIT.md`](analysis/MULTIHIT.md) for why bulk sequencing cannot
answer a phasing question,
[`analysis/LITERATURE_E_RESULTS.md`](analysis/LITERATURE_E_RESULTS.md) for what
the papers say about each of those results, [`analysis/FINDINGS.md`](analysis/FINDINGS.md) for the exploration,
[`analysis/LITERATURE_REVIEW.md`](analysis/LITERATURE_REVIEW.md) for what the
papers changed, and [`reference/README.md`](reference/README.md) for data
provenance.

## Data

All data used is open or CC0. The main calibration set is Watson et al. 2020
(*Science*), released under CC0 via Dryad: 1,674 driver variants with allele
frequency, age, gene and cohort, of which 1,012 carry a usable age.

A second calibration set was added in stage F: Fabre et al. 2022 (*Nature*),
released under CC BY 4.0 via figshare — 4,287 measurements of 994 clones in 394
individuals, with person-level identifiers and two to five timepoints each.

No controlled-access data is used. The raw sequencing behind Fabre 2022 sits in
the EGA under controlled access and was not requested or used; only the open
derived variant table.
