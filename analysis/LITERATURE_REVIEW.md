# What the literature says, and what it changes here

Three papers read in full or in detail, chosen because each answers a question
our model had left open.

| Paper | Method | Answers |
|---|---|---|
| Watson et al. 2020, *Science* | cross-sectional VAF spectrum | fitness per variant, mutation rates |
| Mitchell et al. 2022, *Nature* | phylogenetic reconstruction + ABC | stem-cell numbers, dynamics after 70 |
| Fabre et al. 2022, *Nature* | longitudinal, 697 clones over 13 years | whether fitness is constant over life |

---

## 1. Our population size was wrong by a factor of two

We chose N = 50,000 stem cells with a one-year division time, giving
Nτ = 50,000. That was a guess.

**Two independent methods converge on Nτ ≈ 100,000.**

Watson infers it from the VAF spectrum: Nτ ≈ 100,000 ± 30,000 years. Mitchell
infers it from phylogenies built on 3,579 whole-genome-sequenced colonies,
plus telomere shortening at 30 bp per year — and lands on the same plateau of
about 100,000 HSC-years, with 20,000 to 200,000 long-term stem cells and
τ ≈ 1 year.

Agreement between a frequency-spectrum method and a phylogenetic one is
strong evidence: they share almost no assumptions.

**This matters because Nτ sets the strength of drift.** Halving it doubles the
randomness, and the model then needs a different fitness to reproduce the same
observed spread. Our s = 0.10 was conditional on a weaker drift than reality.

*Action taken: the sweep was re-run with N = 100,000, and it changed the
answer in the direction the literature predicts.*

| | N = 50,000 | N = 100,000 |
|---|---|---|
| inferred s | 0.101 / yr | **0.129 / yr** |
| distance to observed | 0.177 | **0.075** |
| simulated median VAF | 0.0460 | 0.0559 |
| observed median VAF | 0.0509 | 0.0509 |

The fit improved 2.4-fold, and s moved from just below Watson's per-variant
range into it — 12.9% sits among their top twenty (11.2 to 23.1%) and near
DNMT3A R882H at 14.8%.

A correction that both improves the fit and moves the estimate toward
independently published values is the shape a correct one has.

---

## 2. Cross-sectional and longitudinal methods disagree

This is the most consequential finding, and it is not about our code.

| Gene / variant | Watson 2020 (cross-sectional) | Fabre 2022 (longitudinal) |
|---|---|---|
| DNMT3A | R882H at 14.8% / yr | ~5% / yr overall |
| TET2 | — | 6.8% / yr |
| SRSF2 | P95R at 23.1% / yr | P95H **> 50%** / yr |
| splicing genes | 15–23% / yr | 5.4% / yr |

The two families of estimate are not measuring quite the same thing — Watson
reports per-variant hotspots, Fabre reports gene-level averages across 697
clones tracked in the same people for a median of 13 years — but the gaps run
in both directions and are large.

Mitchell, from phylogenies, reports a gamma-distributed spectrum with median
s of 5–10% per year and a heavy tail above 10%, with the 46 largest clones at
10–30%.

**Where our estimate sits.** Our pooled s = 0.10 per year is at the top of
Mitchell's median band, above Fabre's gene-level rates, and below Watson's
hotspot values. For a single coefficient fitted to a detection-truncated
mixture, that is a defensible place to be — but the spread between published
methods is itself the honest headline: **fitness estimates in this field carry
method-dependent uncertainty comparable to the effect being measured.**

---

## 3. Fitness is not constant over life

Our model assumes a clone's fitness `s` is fixed from birth. So does Watson's.

Fabre's longitudinal data contradicts it:

- **DNMT3A clones expand preferentially early in life and slow in old age**,
  in what the authors call "an increasingly competitive oligoclonal
  landscape".
- **Splicing-gene mutations drive expansion only later in life.**
- Clones carrying the *same* mutation differ in growth rate by about
  ±5% per year — which is most of the effect for a slow driver.

So fitness depends on age, on which other clones are present, and on
something individual beyond the mutation itself.

The reassuring part: 92.4% of clones still expanded at a stable exponential
rate over the observation window. Constant fitness is a serviceable
approximation within a decade; it is not one across a lifetime.

---

## 4. After 70, the marrow is a different place

Mitchell's phylogenies show that after age 70 haematopoiesis reorganises:

- **12 to 18 independent expanded clones** per person, against one or fewer in
  younger adults
- those clones produce **30 to 60% of all blood**
- clonal diversity falls off a cliff

This is the regime where our model's competition term stops being a formality.
While clones are small, competition is negligible and a Moran process and a
branching process agree. Once a dozen clones hold half the marrow, they are
competing with each other, not just with wild-type cells — and the two model
families part company.

It also explains Fabre's slowdown of DNMT3A in old age: a clone growing into
a marrow that is already half-claimed grows slower, not because its mutation
changed but because its neighbourhood did.

---

## What this changes in the project

**Done now**

1. **N corrected to 100,000**, consistent with Nτ ≈ 100,000 and τ ≈ 1 year
   from two independent methods. Sweep re-run.

**Next, in order of expected impact**

2. **Fitness per variant** (stage D). Now with two independent targets to
   check against: Watson's per-variant table and Fabre's gene-level rates.
   Landing between them would be a meaningful result; landing outside both
   would mean something is wrong.

3. **Age-dependent fitness.** Fabre's finding is directly implementable: let
   `s` decline as the mutant fraction of the marrow rises. Our model already
   computes that fraction every time slice — the competition term W — so the
   hook exists.

4. **Compare model families.** Re-run the same inference under a branching
   process. Given how much published estimates disagree by method, measuring
   our own method-dependence is not optional.

**Not worth doing**

Chasing a tighter posterior for the mutation rate on this dataset. Stage C
showed it is unidentifiable without a screening denominator, and no amount of
extra simulation fixes missing data.

---

## Sources

- Watson CJ, Papula AL, Poon GYP, Wong WH, Young AL, Druley TE, Fisher DS,
  Blundell JR (2020). *The evolutionary dynamics and fitness landscape of
  clonal hematopoiesis*. **Science** 367:1449–1454.
  [publisher](https://www.science.org/doi/10.1126/science.aay9333) ·
  [open PDF](https://web.stanford.edu/group/dsfisher/papers/pdf/watson_et_al_2020.pdf)
- Mitchell E, Spencer Chapman M, Williams N, et al. (2022). *Clonal dynamics
  of haematopoiesis across the human lifespan*. **Nature** 606:343–350.
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9177428/) ·
  [publisher](https://www.nature.com/articles/s41586-022-04786-y) ·
  [code](https://github.com/emily-mitchell/normal_haematopoiesis)
- Fabre MA, de Almeida JG, Fiorillo E, et al. (2022). *The longitudinal
  dynamics and natural history of clonal haematopoiesis*. **Nature**
  606:335–342.
  [publisher](https://www.nature.com/articles/s41586-022-04785-z) ·
  [open access](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9177423/) ·
  [preprint](https://www.biorxiv.org/content/10.1101/2021.08.12.455048v1)
- Gene-specific fitness from longitudinal data, **Nature Medicine** 2022 —
  <https://www.nature.com/articles/s41591-022-01883-3>
