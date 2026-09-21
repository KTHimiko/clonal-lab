# What the literature says about each stage E result

A second pass through the literature, one section per result the stage produced.
The first pass ([`METHOD_BIAS_LITERATURE.md`](METHOD_BIAS_LITERATURE.md))
checked the claims stage E made about *other people's* numbers. This one checks
the claims it made about *mechanism*.

Three results are corroborated, one from an unexpectedly close direction. One
result gains a competing explanation that neither we nor the paper we were
arguing with had considered.

---

## 1. Deceleration is driven by the population's average fitness

**Our finding** ([`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md)): regressing per-clone
deceleration on three candidates in a mixed-fitness cohort gave the carrier's
total fitness-weighted burden growth R² = 0.956, the clone's own size 0.271, and
its own fitness **0.003**. We described this as "the Moran normalisation
surfacing as the only term that matters."

**Mon Père, Terenzi & Werner (2026)**, *Cancer Discovery* 16(9):1765–1779,
derive the same thing analytically:

> "the rate at which a clone grows […] is proportional to the difference between
> its own innate fitness and the **average fitness in the HSC population**"

and state the consequence we measured:

> "a clone with positive innate fitness will always outcompete WT HSCs and
> therefore grows in the absence of other clones, [but] a polyclonal stem cell
> population containing fitter clones can force its decline."

**This is the strongest external result in the whole stage.** An independent
group, working from trajectory data rather than simulation, derived the
mechanism our regression recovered numerically. Our `W` term *is* their average
population fitness. They also record the observation that motivated our problem
2 — "variants typically saturate below fixation, many variant trajectories
decrease late in life."

They name it **clonal interference**, which is the established term. We had been
calling it saturation, which is narrower and slightly wrong: saturation suggests
a ceiling the clone approaches, whereas the mechanism is being overtaken by
fitter neighbours. Adopt their term.

---

## 2. Detection thresholds bias what can be concluded

**Our finding** ([`METHOD_BIAS.md`](METHOD_BIAS.md)): moving only the detection
floor, on the same cohorts and the same variants, changed the age slope by a
factor of 4.1. The truncation removes small clones, and small clones are what
make the young end of the line low.

**Independent confirmation in a different currency.** A 2025 medRxiv analysis of
CHIP detection finds that calling CHIP from genome or exome sequencing rather
than targeted deep sequencing captures only about **9% of the true disease
hazard**, with power for a hazard ratio of 1.5 / 2.0 / 2.5 / 3.0 running
**8% / 28% / 53% / 70%**.

Same mechanism, different endpoint: theirs is clinical risk, ours is a fitness
estimate. Both are dominated by which clones the instrument could see. It is
worth knowing that the effect is large enough to matter in the clinical
literature, not only in inference.

---

## 3. Independent methods agree once compared like for like

**Our finding** ([`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md)): three
detection-conditioned methods land between 0.11 and 0.20 for DNMT3A; the
longitudinal estimate at 0.062 is the outlier.

Two pieces of external support:

- **Williams et al. (2023)**, *Estimating single cell clonal dynamics in human
  blood using coalescent theory* — the methodology behind the `phylofit`
  estimates we used. They report phylofit correlating with ABC-based estimates
  at **r = 0.96**. The third method in our comparison is itself cross-validated
  against a fourth.
- **Mon Père et al. (2026)** run Bayesian inference on longitudinal trajectories
  and on single-cell genetic heterogeneity independently and report that
  "**inferences on both data converge**."

So agreement across method families is not a coincidence of our pairing; it is
reproducible when the quantities are matched.

---

## 4. Detected clones are not the spectrum — a third instance

This kept being the error underneath the other errors. It now has a third
independent illustration.

| paper | underlying spectrum | detected / expanded |
|---|---|---|
| Mitchell 2022 | gamma, mode 5–10% / yr | 46 largest clades, 10–30% / yr |
| Watson 2020 | >90% of DNMT3A nonsynonymous effectively neutral | top-20 variants, 11–23% / yr |
| **Mon Père 2026** | **innate fitness η = 0.07 (+0.04 / −0.03) per year** | **elevated-fitness clones 0.26–0.49, one at 1.03** |

Three papers, three pairs of numbers that look like disagreements and are not.
Every time this project compared across that line, it got a wrong answer — the
R882-versus-gene-average pair, and the Mitchell-gamma-versus-spectrum-fit pair.

---

## 5. TET2 and age — a third mechanism, which we did not model

This is the substantive change.

Stage E modelled two candidates for the TET2 age effect: the mutation's own
fitness ramping, and the wild type declining. **Mon Père et al. supply a third,
and it fits the observation better than either:**

> "*DNMT3A* variants were enriched for single-hit clones, whereas ***TET2***,
> *ASXL1*, *JAK2*, *SF3B1* and *SRSF2* showed enrichment for **multihit
> evolution**"

> "The fittest clones emerge predominantly later in life in accordance with a
> **multistep evolutionary process**"

> "higher fitness trajectories had a higher likelihood of being detected
> alongside at least one other driver (62% of highest-fitness trajectories vs
> 42%)"

**A TET2 clone that acquires a second driver later in life accelerates. Under a
single-mutation model that is indistinguishable from "TET2 fitness rises with
age."** And it is gene-specific in exactly the way the observation is: DNMT3A,
enriched for single-hit clones, would show no such acceleration.

Why this beats our two candidates:

| mechanism | predicts TET2-specific? | needs an extra assumption? |
|---|---|---|
| wild type declines | **no** — lifts every clone (we measured ×2.04 for the untouched class) | — |
| TET2 fitness ramps | yes, by construction | yes: the ramp itself, with no mechanism |
| **multihit** | **yes, from measured gene-level biology** | **no** |

The multihit explanation requires nothing new: second drivers arrive at a rate
we could estimate, and TET2's enrichment for them is measured, not posited.

**What it changes in our conclusion.** The negative result stands —
[`TET2_AGE.md`](TET2_AGE.md) showed the apparent effect does not survive
stratifying by cohort and could not have been detected at n = 36 regardless. But
the framing was wrong in one respect: we presented the choice as between two
mechanisms that a relative-fitness model cannot distinguish, and concluded the
question needs an absolute measurement of wild-type output. **The multihit route
is distinguishable without that**, because it predicts something observable:
accelerating TET2 clones should carry a second driver and non-accelerating ones
should not — a question about variant co-occurrence rather than about absolute
haematopoietic output.

**But not in this dataset, and it is worth being exact about why.** The Watson
table has five columns — `VAF, age, variant, gene, study` — and **no person
identifier**. Two rows cannot be assigned to the same individual, so
co-occurrence is not computable. Using `(study, age)` as a substitute key fails
on inspection: 192 of the 479 pairs hold more than one variant, one of them 24,
and that group spans ASXL1, CSF3R, CUX1 and six different DNMT3A variants —
obviously several people who happen to share an age in a cohort of hundreds, not
one person with nine drivers.

This is the third time the same structural loss has bitten. Stage C found the
mutation rate unidentifiable because the table carries no screening denominator.
Stage E found the one working age-based estimator unusable for the same reason.
Now the multihit test fails because person-level linkage was aggregated away.
**An aggregated table answers questions about variants; every question about
people needs the cohorts it was built from.**

---

## 6. Still no prior art for the design-bias measurement itself

Searched again with different terms. Nothing found that simulates a cohort with
known `s` and measures what each study design recovers. The closest remains
BESTish (2025), which unifies cross-sectional and longitudinal data in one
Bayesian framework rather than measuring the bias of each, and which explicitly
does not analyse how detection thresholds bias estimates.

"Not found" is weaker than "does not exist", and this is a large literature.

---

## What this changes in the project

**Done here**

1. **Terminology.** "Saturation" becomes **clonal interference** where the
   mechanism is being overtaken by fitter clones. The established term is more
   accurate than ours.
2. **`TET2_AGE.md` gains the multihit alternative**, and the claim that the
   question requires an absolute measurement is narrowed — it is true for the
   ramp-versus-wild-type pair, and not true for multihit.

**Next, in order of expected value**

3. **Test the multihit hypothesis — which needs data we do not have.** The
   question is whether TET2 variants co-occur with a second driver more often in
   older carriers. It cannot be asked of the Watson table, which has no person
   identifier. It can be asked of any cohort that publishes per-individual
   variant lists, and Mon Père's own analysis is built on exactly that, so the
   route is to their data rather than to a new computation on ours.
4. **Implement second hits in the model.** A clone that acquires another driver
   gains fitness. The per-clone fitness matrix already supports it; what is
   missing is a second mutation process acting on existing clones rather than on
   wild type.
5. **Compare our inferred fitness against η = 0.07 per year** — Mon Père's
   underlying innate-fitness mean — being careful, this time, to compare it
   against an underlying spectrum and not against our detection-conditioned
   estimates.

---

## Sources added in this pass

- Mon Père NV, Terenzi F, Werner B (2026). *The Evolution of Polyclonal
  Competition in Aging Hematopoiesis*. **Cancer Discovery** 16(9):1765–1779.
  [publisher](https://aacrjournals.org/cancerdiscovery/article/16/9/1765/787658/The-Evolution-of-Polyclonal-Competition-in-Aging) ·
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC13531000/) ·
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/42084533/)
- Williams N, Lee J, Mitchell E, et al. (2023). *Estimating single cell clonal
  dynamics in human blood using coalescent theory*. **bioRxiv**.
  [preprint](https://www.biorxiv.org/content/10.1101/2023.02.24.529817v2.full)
  — the method behind the `phylofit` estimates used in `OPEN_PROBLEMS.md`.
- *Detection of clonal hematopoiesis of indeterminate potential via genome or
  exome sequencing profoundly underestimates disease associations* (2025).
  **medRxiv**.
  [preprint](https://www.medrxiv.org/content/10.1101/2025.08.11.25333294v1) ·
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12948416/)
- *BESTish: A Diffusion-Approximation Framework for Inferring Selection and
  Mutation in Clonal Hematopoiesis* (2025).
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC12879653/)
- *HSC Niche Dynamics in Regeneration, Pre-malignancy, and Cancer: Insights From
  Mathematical Modeling* (2023). **Stem Cells** 41(3):260.
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC10020982/)
- *Impact of the clonal architecture and probable biallelic inactivation of TET2
  mutations in myelodysplastic neoplasms* (2026). **Blood Cancer Journal**.
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC13221467/)
  — biallelic TET2 accumulating with age, consistent with the multihit route.
