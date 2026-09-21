# Checking our estimate against the published values

Our stage-C result was `s = 0.10 per year`, inferred by ABC over a 180-point
grid. This document checks it against Watson et al. 2020, the paper our
calibration data comes from.

Source: Watson, Papula, Poon, Wong, Young, Druley, Fisher & Blundell (2020),
*The evolutionary dynamics and fitness landscape of clonal hematopoiesis*,
**Science** 367:1449–1454. Read from the open PDF hosted by the Fisher lab
at Stanford.

---

## Their fitness estimates

Figure 2A of the paper lists inferred fitness effects for the twenty most
commonly observed variants, in **percent per year** — the same unit we used.

| Variant | s (% per year) |
|---|---|
| SRSF2 P95R | 23.1 |
| SF3B1 K700E | 22.9 |
| **DNMT3A R882C** | **18.7** |
| SRSF2 P95H | 16.1 |
| DNMT3A R729W | 16.0 |
| DNMT3A R326C | 15.8 |
| **DNMT3A R882H** | **14.8** |
| JAK2 V617F | 14.6 |
| IDH2 R140Q | 14.5 |
| DNMT3A R736H | 14.1 |
| … | … |
| DNMT3A P904L | 11.2 |

They also classify the whole spectrum into three bands: **low** (s ≤ 4%),
**moderate** (4–10%) and **high** (s ≥ 10%) per year, and report that roughly
40% of variants confer moderate-to-high fitness.

## Where our number lands

**Our pooled estimate of 10% per year sits at the boundary between their
moderate and high classes, and below every one of the top twenty variants
(which run 11.2% to 23.1%).**

That is the right place for it, and the reason is instructive.

We fitted a **single s shared by every clone**. The observed set we calibrated
against is truncated at VAF ≥ 0.0192, which preferentially retains variants
fit enough to have grown large — so our pooled estimate reflects the detected
subset, enriched for fitness, rather than the full spectrum. Meanwhile the
paper reports that in DNMT3A, TET2 and ASXL1 over 90% of nonsynonymous
variants fall in the **low** fitness class.

So a pooled 10% for the detected subset is coherent with a distribution whose
bulk is below 4% and whose visible tail runs 11–23%.

## The mutation rate, unexpectedly

Our μ posterior was wide — 6.6 × 10⁻⁷ to 1 × 10⁻⁵ per cell per year, spanning
91% of the grid — with a **median of 3.4 × 10⁻⁶**.

The paper states that mutations conferring s > 4% per year arise at a rate of
about **4 × 10⁻⁶ per cell per year**.

Our weakly-identified median lands within a factor of 1.2 of their estimate.
That agreement is partly luck given how flat our posterior is, but it does say
the grid was centred on the right order of magnitude rather than by accident.

---

## Two discrepancies worth stating

### 1. They used a branching process; we used a Moran process

The paper builds "a simple stochastic branching model of HSC dynamics". Ours
is a Moran process with a fixed population size.

The two are not interchangeable. A branching process lets the population grow
or shrink; a Moran process holds it constant by construction. They agree in
the regime where clones are small relative to N — which is most of life — and
diverge as a clone approaches taking over the marrow.

This also revives the factor-of-two caution from stage A: fixation and escape
probabilities differ between model families, and a fitness estimate is only
meaningful alongside the model that produced it.

### 2. Our N is a choice; theirs is a product they could not separate

We set N = 50,000 stem cells and implicitly took one year between divisions.

The paper is explicit that population-genetic analysis **cannot separate N
from τ** (the time between symmetric self-renewal divisions) — only their
product. They infer **Nτ ≈ 100,000 ± 30,000 years**, and from it bound the
stem-cell count between **25,000 and 1.3 million**.

Our N = 50,000 with τ = 1 year gives Nτ = 50,000: inside their bounds on N,
but **half their inferred Nτ**.

That matters, because Nτ is what actually sets the strength of drift. Holding
Nτ at 100,000 instead would change how much randomness the model produces, and
therefore the fitness needed to explain the observed spread. Our s = 0.10
should be read as *conditional on Nτ = 50,000*, not as a free-standing number.

---

## What this verification changes

**It does not invalidate the result.** The pipeline works, the inference runs,
and the number lands where the literature says it should for the subset we can
see.

**It does sharpen what the number means.** `s = 0.10 per year` is a pooled
estimate for detected clones, under a Moran model with Nτ = 50,000. Every one
of those qualifiers is load-bearing.

**It names the next corrections**, in order of how much they would change the
answer:

1. **Fitness per variant instead of a single pooled s** — stage D. The paper's
   Figure 2A is a ready-made target: fit DNMT3A R882H separately and check
   against 14.8%.
2. **Treat Nτ as the parameter**, not N, and sweep it over their inferred
   range rather than fixing it.
3. **Compare model families**: rerun the same inference under a branching
   process and see how much of the estimate is the model rather than the data.

---

## Sources

- Watson et al. 2020, *Science* — [publisher](https://www.science.org/doi/10.1126/science.aay9333)
  · [open PDF](https://web.stanford.edu/group/dsfisher/papers/pdf/watson_et_al_2020.pdf)
- Dryad dataset (CC0) — <https://doi.org/10.5061/dryad.83bk3j9mw>
- Fabre et al. 2022, *Nature*, longitudinal dynamics —
  <https://www.nature.com/articles/s41586-022-04785-z>
- Mitchell et al. 2022, *Nature*, dynamics across the lifespan —
  <https://www.nature.com/articles/s41586-022-04786-y>
- Review of gene-specific fitness from longitudinal data, *Nature Medicine*
  2022 — <https://www.nature.com/articles/s41591-022-01883-3>
