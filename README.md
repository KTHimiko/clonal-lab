# clonal-lab

Agent-based modelling of clonal competition in ageing tissue, calibrated
against public data on clonal haematopoiesis.

As we age, somatic mutations let some blood stem cells outcompete their
neighbours. The resulting clonal expansions — **clonal haematopoiesis** — are
common after 60 and roughly double cardiovascular risk.

The project set out to ask how strong those fitness advantages are. It ended up
asking something more useful: **how much to trust each way of measuring them.**
Published estimates for the same gene differ by up to a factor of three, and no
dataset can say which is right — but a simulator whose answer is known by
construction can, and a cohort measured repeatedly can check it.

**Start with [`SYNTHESIS.md`](SYNTHESIS.md)**: every result, every correction,
and every question the data could not answer, in one read.

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
| **H** — TET2 longitudinal | each clone as its own control, over 13 years | **done** |
| **I** — replication | pre-registered, then abandoned for a stated reason | **done** |
| **J** — tailored grid | two estimators on one cohort; one has no inverse | **done** |
| **K** — mixture, N free | can any model version fit both observables? | **done** |
| **L** — mechanisms | niche structure and fluctuating fitness, both tested | **done** |
| **M** — composition | fits the target, refuted by a statistic it predicts | **done** |
| **N** — clone-specific switch | closest yet; the bracket closes from both sides | **done** |

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
.venv/bin/python analysis/12_tet2_longitudinal.py  # TET2 and age, within clone
.venv/bin/python analysis/13_lothian_replication.py # why the replication cannot be run
.venv/bin/python analysis/14_fabre_grid.py     # two estimators, one cohort, one grid
.venv/bin/python analysis/15_mixture_fit.py    # a fitness distribution, with N free
```
## What it found

**The full account is in [`SYNTHESIS.md`](SYNTHESIS.md)** — every result, every
correction, and every question the data could not answer, in one read. The
per-stage documents in [`analysis/`](analysis/) hold the detail.

The short version:

**Study designs do not measure the same thing.** Against a simulated cohort
whose fitness is known by construction, a fit to the VAF spectrum recovers 103%
of the truth, thirteen years of follow-up recovers 43%, and regressing log(VAF)
on the host's age recovers 18% — and that last one is non-monotonic in `s`, so it
has no inverse at all.

**The ordering holds in 394 real people.** Fabre's cohort allows all three
designs on the same individuals: 0.110, 0.0455 and 0.0024 per year. Spectrum ÷
longitudinal is 2.42 observed against 2.40 predicted, and 2.42 between the two
published papers.

**One result needs no model.** In 770 real trajectories, restricting to clones
already detectable at enrolment drops the measured growth rate from 0.0703 to
0.0455 per year — a 35% loss, with non-overlapping intervals. A study that
enrols the clones it can already see is measuring the ones that have most nearly
finished growing.

**Clones slow because the marrow fills.** Per-clone deceleration is explained by
the carrier's total clonal burden growth (R² = 0.956), not by the clone's own
fitness (R² = 0.003). Mon Père et al. derive the same relation analytically.

**Three independent methods agree once compared like for like.** Spectrum fits
and phylogenetic reconstruction both put DNMT3A between 0.11 and 0.20 per year;
the longitudinal estimate is the outlier at 0.062, low by the predicted amount.

**TET2 clones are still growing after 75; DNMT3A clones are not.** 234
trajectories, each clone its own control. The between-gene comparison is
borderline, three measures were tested, and the study had **0.56 power** —
computed afterwards while planning a replication, which is a retrospective
downgrade of the result rather than a footnote to it. The replication was
pre-registered and then abandoned: the candidate cohort holds 8 TET2 clones
against Fabre's 114.

**The model cannot fit both observables of one cohort.** Letting fitness be a
distribution rather than a number, and letting the stem-cell population range
over eightfold, no combination of the two reproduces both the sizes of clones and
their growth rates. The diagnosis is single: the model's clones come out about
twice as large as the observed ones, and because measurement noise scales as
`1/√(VAF·depth)`, larger clones are measured more precisely and their rate
distribution comes out too narrow. One failure, not two — and 51 technical
replicate triplets rule out the obvious alternative, showing the assay is if
anything *less* noisy than binomial.

**And four questions could not be answered**, each on a different missing column:
the mutation rate and the one unbiased age-based estimator both need a screening
denominator; the cross-sectional TET2 test had a power ceiling of 0.79; and
multihit needs phase, which bulk sequencing does not measure.

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
