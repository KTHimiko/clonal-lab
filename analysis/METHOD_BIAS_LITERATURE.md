# Checking stage E against the literature

Stage E measured how much each study design under-recovers a known `s`. Every
number came from our own simulator, so the obvious failure mode is that the
biases are real only inside the model. This document checks the four claims
that depend on facts outside it.

Stage C's literature check found a parameter wrong by a factor of two and
improved the fit 2.4-fold. That precedent is why this is not optional.

---

## C1 — Is Watson's estimator really a fit to the VAF spectrum?

**Confirmed, verbatim.** The whole reframing of stage E rests on this, and it
is stated plainly in the paper:

> "Because the characteristic maximum VAF, f, depends on the fitness effect, s,
> by estimating f from the VAF spectrum, we can infer a variant's fitness."
> — Watson et al. 2020, *Science* 367:1449–1454

No regression against age appears anywhere in their inference. Calling Watson
"the cross-sectional estimate" and contrasting it with Fabre as though the axis
were sampling design is wrong, and stage D said exactly that. Corrected.

### The unexpected part: Watson has an age-based cross-check, and it works

They do use age — as an independent validation, with a different estimator than
the one we tested:

> "the age prevalence of these variants does increase linearly with age […] the
> rate of this increase is consistent with a fitness effect of s ≈ 14% per year,
> which is in agreement with estimates inferred from the VAF distribution."

That is **prevalence** against age, predicted to rise at rate `2Ntμs` — not
log(VAF) against age. 14% against 15% from the spectrum: two estimators, one
answer.

So stage E's conclusion needs narrowing. *Age carries usable information about
fitness.* What fails is one specific way of extracting it — the slope of
log(VAF) among detected clones. The prevalence slope is a different functional
and, on Watson's own data, an unbiased one.

**And we cannot use it.** Prevalence needs a screening denominator: how many
people were tested per age band, not just how many carried a variant. Stage C
established this dataset has none, which is why the mutation rate came out
unidentifiable. The same missing denominator blocks the one age-based estimator
that would have worked. Two findings, one cause.

---

## C2 — Is the saturation mechanism real, or an artifact of the Moran term?

**Supported, and it is the authors' own explanation.**

Fabre attributes their deceleration to "an increasingly competitive oligoclonal
landscape" — clones growing into a marrow that is already claimed. That is the
Moran normalisation in words: the model's competition term was not invented to
rescue this result, it is the standard account of the phenomenon.

Independently, the BESTish diffusion-approximation framework (2025) derives
mean-field VAF trajectories that are **logistic**, not exponential. Saturation
falls out of the analysis rather than being imposed.

### But two observations do not fit

**Fast drivers decelerate least, and the model predicts the opposite.** Fabre
reports deceleration "most marked for clones bearing mutations in DNMT3A, BRCC3
and TP53", while "fast-growing clones harbouring U2AF1, SRSF2 P95H, PTPN11 or
IDH1 mutations" showed "almost no deceleration."

Our model says the saturation bias **grows with `s`**, because a fitter clone
reaches the saturating regime sooner — that is exactly why our longitudinal
estimate collapses from 43% to 4% of the truth as `s` goes from 0.08 to 0.24. A
clone growing at 50% per year from a detectable size should saturate within the
observation window and show dramatic deceleration. Fabre sees almost none.

Candidate explanations, none verified: those clones were caught earlier in
their trajectory and were still small through the window; the marrow tolerates
a larger mutant fraction before competition bites than the Moran normalisation
implies; or their follow-up ended before saturation. **This is the sharpest
empirical challenge to the mechanism in stage E, and it is unresolved.**

**TET2 grows faster in older individuals.** Fabre finds "age was a significant
factor specifically for TET2-mutant clones, which grew faster in older
individuals." Our model has no age-dependent fitness at all, and this runs
opposite to the deceleration we do model. Any TET2 number from our pipeline
carries a bias we cannot currently estimate.

---

## C3 — Does a third, independent method break the tie?

**Yes, and it does not side with us.** This is the open problem.

Mitchell et al. 2022 reconstruct clonal histories from phylogenies built out of
3,579 single-cell-derived colonies — a method that shares no machinery with
either VAF spectra or serial VAF measurement. They report growth rates "ranging
from 5% (DNMT3A and TP53) to more than 50% per year (SRSF2 P95H)."

| Method | DNMT3A | Family |
|---|---|---|
| Watson 2020 | ~15% / yr | VAF spectrum |
| Fabre 2022 | 6.2% / yr | longitudinal |
| Mitchell 2022 | ~5% / yr | phylogenetic |

Two of the three land near 5–6%. Stage E argued the spectrum fit recovers the
truth and follow-up under-recovers; the phylogenetic estimate sides with
follow-up.

**Our experiment cannot adjudicate this**, and it is important to say why:
stage E validated the spectrum fit against a grid generated by our own model.
That checks the inversion, not the biology (section 5 of `METHOD_BIAS.md` says
so). An upward bias present in real data but absent from our simulator would be
invisible to the entire experiment by construction.

Two reconciliations are available and we cannot choose between them:

1. **Mixture, not bias.** Watson's headline values are per-variant hotspots —
   the fit tail. Watson themselves report that in DNMT3A over 90% of
   nonsynonymous variants are effectively neutral. A gene-level average over
   that distribution *should* be far below the hotspot value, so Mitchell's 5%
   and Watson's 15% may be measuring different subsets rather than disagreeing.
2. **Shared bias.** Mitchell's cohort is elderly, and a single exponential
   fitted over a trajectory that has been saturating returns an average dragged
   below the early-life rate — the same effect stage E measured in the
   longitudinal design.

Distinguishing them means restricting all three methods to the same variant
class, which the published summaries do not permit. **Flagged as the first
thing to settle in any future stage.**

---

## C4 — Is the published gap the size stage E predicts?

**Yes, on the one comparison that is like for like — and stage E was quoting
the wrong pair.**

Stage D and the first draft of stage E used DNMT3A R882H at 0.148 (Watson)
against "~5%" (Fabre), a ratio of 3.0. That pair mixes a per-variant hotspot
with a gene-level average, which is precisely the mixture problem above.

Fabre state the clean pair themselves, in their Discussion:

> "expansion of clones was substantially faster in younger individuals (15.0%
> per year) compared with older individuals (6.2% per year)"

Both numbers are DNMT3A gene-level; the 15.0% is Watson's spectrum estimate,
the 6.2% is Fabre's own longitudinal measurement.

| | |
|---|---|
| Published ratio, DNMT3A, spectrum ÷ longitudinal | 15.0 / 6.2 = **2.42** |
| Ratio stage E predicts from design bias alone | 103% / 43% = **2.4** |

The agreement is closer than the experiment deserves — three seeds, one
mutation rate, a single follow-up window — and should be read as consistency,
not as a fourth-significant-figure match. But it is the right size, on the
right pair, in the right direction.

**Fabre reach the same conclusion by a different route.** They read the gap as
DNMT3A clones genuinely growing faster in early life, not as either method
being wrong. Stage E reproduces that gap in a simulation where `s` is constant
by construction, with competition as the only mechanism. The two readings agree
on the mechanism and differ on what to call it: a clone whose growth slows
because the niche filled has not lost fitness, and calling the late-window rate
its "fitness" is the error.

---

## C5 — Has anyone already run this experiment?

**Not found, which is weaker than "nobody has."**

Searches for simulation-with-known-ground-truth comparisons of CH study designs
returned nothing directly equivalent. The closest work, BESTish (2025), builds
one Bayesian framework that ingests both cross-sectional and longitudinal data
rather than measuring the bias of each. It explicitly does not analyse how
detection thresholds bias parameter estimates.

The clone-age problem is, however, stated in the literature as a known
limitation rather than a discovery of ours:

> "Inferring fitness from a single timepoint creates additional uncertainty
> about whether a mutation has arisen recently and has grown rapidly […] or
> arose a long time ago and has grown slowly."

Stage E's contribution is putting a number on it, not noticing it.

---

## What this changes

**Corrected in `METHOD_BIAS.md`**

1. The comparison pair is DNMT3A 15.0 vs 6.2 (ratio 2.42), not R882H 0.148 vs
   0.050 (ratio 3.0). The old pair was not like for like.
2. "Regression on host age" is too broad a claim. What fails is the **slope of
   log(VAF) among detected clones**. Watson's prevalence-against-age estimator
   uses the same sampling design and is unbiased on their data.
3. The saturation mechanism now carries a named counter-observation: fast
   drivers should decelerate most and are reported to decelerate least.

**Open, in priority order**

1. **Reconcile with the phylogenetic estimates.** Requires all three methods on
   one variant class. Until then, stage E's claim that the spectrum fit is the
   unbiased one rests on self-consistency, and a third method disagrees.
2. **Explain the fast-driver discrepancy**, or find that our competition term
   saturates too early. This is directly testable: fit our own simulated fast
   drivers over a 13-year window and check whether they decelerate as sharply
   as the model implies.
3. **Age-dependent fitness for TET2**, in the direction Fabre reports.

---

## Sources

- Watson CJ, Papula AL, Poon GYP, Wong WH, Young AL, Druley TE, Fisher DS,
  Blundell JR (2020). *The evolutionary dynamics and fitness landscape of
  clonal hematopoiesis*. **Science** 367:1449–1454.
  [publisher](https://www.science.org/doi/10.1126/science.aay9333) ·
  [open PDF](https://web.stanford.edu/group/dsfisher/papers/pdf/watson_et_al_2020.pdf)
- Fabre MA, de Almeida JG, Fiorillo E, et al. (2022). *The longitudinal
  dynamics and natural history of clonal haematopoiesis*. **Nature**
  606:335–342.
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9177423/)
- Mitchell E, Spencer Chapman M, Williams N, et al. (2022). *Clonal dynamics of
  haematopoiesis across the human lifespan*. **Nature** 606:343–350.
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9177428/)
- Robertson NA, Latorre-Crespo E, Terradas-Terradas M, et al. (2022).
  *Longitudinal dynamics of clonal hematopoiesis identifies gene-specific
  fitness effects*. **Nature Medicine** 28:1439–1446.
  [publisher](https://www.nature.com/articles/s41591-022-01883-3) ·
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9307482/)
- *BESTish: A Diffusion-Approximation Framework for Inferring Selection and
  Mutation in Clonal Hematopoiesis* (2025).
  [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC12879653/)
