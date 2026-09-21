# References — what exists and how to get it

A survey of what is publicly available on clonal dynamics in haematopoiesis,
done before writing any model. Knowing the shape of the real data avoids
building a simulator that talks to nothing.

> [!IMPORTANT]
> None of this is our code. This directory holds third-party material for
> reading and reference, and is **not versioned** with the project.

---

## Obtained

### Mitchell et al. 2022 — code

**Clonal dynamics of haematopoiesis across the human lifespan**, *Nature*
607:343–350. Maps clonal dynamics across the human lifespan and shows the
regime change after roughly age 70, with a few clones dominating.

- Repository: <https://github.com/emily-mitchell/normal_haematopoiesis>
- Cloned into `code/normal_haematopoiesis` — 27 MB, R
- Open-access article: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9177428/>

> [!CAUTION]
> **The repository declares no licence.** With no licence, the legal default
> is all rights reserved: you may read and learn from it, but **not** reuse
> the code in your own project. If something is useful, reimplement it from
> the method described in the paper rather than copying.

**What it teaches, which is the valuable part:**

| Directory | Contents |
|---|---|
| `6_population_modelling` | estimating the stem-cell population size |
| `7_phylofit` | fitting dynamics from phylogenies |
| `8_driver_modelling` | **ABC to infer the fitness effect of driver mutations** |
| `10_simulations_for_figures` | the simulations behind the paper's figures |
| `11_LOY_simulations` | loss of chromosome Y |

The core method is **ABC** (Approximate Bayesian Computation): simulate
thousands of times while varying parameters, compute summary statistics for
each run, and keep the ones that land close to the observed values. The
distribution of the kept parameters is the estimate.

It is massively parallel — each simulation is independent of the others — and
therefore exactly the kind of load a cluster exists for.

They submit with `bsub`, the LSF scheduler at the Sanger Institute. We use
SLURM. The method is scheduler-agnostic; only the submission line changes.

Their simulation engine is the `Rsimpop` R package.

---

## Data obtained — inventory

### Watson 2020 (CC0) — `data/watson2020/`

Eight folders matching the paper's figures, 196 CSVs in total. The file that
matters is in `Maximum_likelihood_estimations/`:

**`all_studies_trimmed_all_genes.csv`** — the project's calibration target.

```
VAF,age,variant,gene,study
0.0076,55,R404*,ASXL1,Acuna2017
```

| | |
|---|---|
| Rows | 1,674 variants |
| Columns | `VAF`, `age`, `variant`, `gene`, `study` |
| Cohorts | Jaiswal2014 (466), Coombs2017 (378), Genovese2014 (197), Acuna2017 (182), Young2019 (158), ZinkWGS (125), McKerrel2015 (112), Desai2018 (31), Young2016 (26) |
| Main genes | DNMT3A (779), TET2 (126), JAK2 (86), ASXL1 (73), TP53 (64), SF3B1 (63), SRSF2 (44), CBL (35) |
| VAF range | 0.0008 to 0.9091 |
| Age range | 5 to 98 years |

> [!WARNING]
> **663 of the 1,674 rows carry `noagedata` instead of an age.** That leaves
> **1,012 usable variants**. Filtering this is the first step of any analysis —
> and missing it would produce a silently wrong age-VAF curve.

Other relevant folders in the same collection:

| Folder | Useful for |
|---|---|
| `Mutation_rate_calculations` | mutation rates by trinucleotide context — **model input** |
| `Maximum_likelihood_estimations` | the authors' own fitness estimates — **ground truth for our inference** |
| `Age_prevalence_of_DNMT3A_R882H_and_R882C_variants` | raw cohorts: McKerrel (112 rows with VAF and age) and Coombs (1,591 clinical samples) |
| `Estimating_fitness_effects_of_infrequently_mutated_sites` | per-site fitness effects |

### Mitchell 2022 (CC BY 4.0) — `data/mitchell2022/`

The full Mendeley package holds **53 files and 5.15 GB uncompressed**, almost
all of it mutation sets for dN/dS analysis — not our problem right now.

Only the useful part was extracted: **`Summary_cut.csv`** (364 KB, 3,592 rows)
with `donor_id`, `age`, `colony_ID`, `cell_type`, `sample_type`, `timepoint`.

---

## Provenance — where these came from

Both require a browser: programmatic attempts hit authentication or bot
protection. Recorded here so the effort is not repeated if a re-download is
ever needed.

### Watson et al. 2020 — data

**The evolutionary dynamics and fitness landscape of clonal hematopoiesis**,
*Science* 367:1449–1454. Estimated the fitness effect of driver variants from
data on roughly 50,000 people, and bounded the number of haematopoietic stem
cells.

- Data: <https://doi.org/10.5061/dryad.83bk3j9mw> — **CC0 licence**
  (public domain, unrestricted use)
- Downloaded: `doi_10_5061_dryad_83bk3j9mw__v20200327.zip`, 566 KB, containing
  eight nested zips
- Open PDF:
  <https://web.stanford.edu/group/dsfisher/papers/pdf/watson_et_al_2020.pdf>

*The Dryad API answers `Unauthorized, must have current bearer token` and the
direct endpoint returns 403. Only the page's "Download dataset" button works.*

### Mitchell et al. 2022 — data matrices

- Mendeley Data: <https://data.mendeley.com/datasets/np54zjkvxr/1>
- Downloaded: `np54zjkvxr-1.zip`, 4.6 GB
- Raw sequencing data: EGA, **controlled access** (EGAD00001007851) — out of
  scope, requires institutional approval

---

## Background reading

Suggested order, from phenomenon to method:

1. **Jaiswal et al. 2014**, *NEJM* — establishes CHIP as a phenomenon and its
   association with mortality and cardiovascular risk.
2. **Genovese et al. 2014**, *NEJM* — the same finding, independently.
3. **Watson et al. 2020**, *Science* — quantifies fitness from allele
   frequency. The calibration target for our model.
4. **Mitchell et al. 2022**, *Nature* — dynamics across the lifespan, and the
   ABC methodology we want to reproduce.

---

## What this settles for the project

**The calibration observable** is the distribution of variant allele frequency
as a function of age. It is what the public data carries and what the model
must predict.

**The inference method** is ABC: there is no closed-form likelihood for an
agent-based model, so comparison happens through summary statistics.

**The computational load** is the sweep: thousands of independent simulations.
That is where the cluster comes in, and why it exists in this project.

**The data constraint holds:** everything we need is open or CC0. The only
controlled material (EGA) is raw sequencing, which we would not use anyway.
