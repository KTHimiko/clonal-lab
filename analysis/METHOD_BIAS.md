# Stage E — which measurement method recovers the truth?

Reproduce with `.venv/bin/python analysis/06_method_bias.py`.
Figures `10_method_recovery.png` and `11_bias_decomposition.png`.
Checked against the literature in
[`METHOD_BIAS_LITERATURE.md`](METHOD_BIAS_LITERATURE.md), which corrected two
claims below and left one open problem.

---

## The question

Stage D compared our estimate with two published ones and found they disagree
on the same gene. For DNMT3A, Watson 2020 gives s = 0.150 per year and Fabre
2022 gives 0.062 — a factor of 2.4. Both are careful work on real cohorts. They
differ in **design**: Watson reads one blood sample from many people of
different ages, Fabre follows the same people for a median of thirteen years.

With real data that disagreement cannot be settled. Nobody knows the true
fitness of a clone in a living person, so measuring harder does not help.

**In a simulation there is an arbiter, because we set `s` ourselves.** Build a
cohort whose fitness is known by construction, apply each measurement design to
it, and compare what comes back against the number we typed in. Any gap is bias
attributable to the design.

This is what the model is for. It is not a better estimate of `s` — it is the
only way to put an error bar on a *method*.

---

## 1. The three designs, named precisely

| Design | What it does | Who uses it |
|---|---|---|
| **Regression on host age** | one draw per person, slope of log(VAF) on the host's age | stage 1 of this project; the first thing anyone reaches for |
| **Longitudinal follow-up** | same people at 65 and 78, one growth rate per clone | Fabre 2022 |
| **VAF spectrum fit** | match the shape of the detected VAF distribution against a simulated grid | stages C and D; the family Watson 2020 belongs to |

**This corrects a shortcut from stage D.** Calling Watson "the cross-sectional
estimate" is true about the sampling and misleading about the estimator. A
single-timepoint study can be read with a slope against age *or* with a
likelihood fit to the spectrum, and those are not the same instrument. Watson
used the second. Only the second is being compared with Fabre here.

---

## 2. What each design recovers

Five true values, three seeds each, 6,000 simulated people, N = 100,000 stem
cells, detection limit VAF ≥ 0.0192, sequencing depth 1,000×.

| true s | regression on age | longitudinal | spectrum fit |
|---|---|---|---|
| 0.08 | 0.019 (−77%) | 0.064 (−19%) | 0.080 (0%) |
| 0.12 | 0.033 (−72%) | 0.087 (−27%) | 0.120 (0%) |
| 0.16 | 0.033 (−80%) | 0.067 (−58%) | 0.162 (+1%) |
| 0.20 | 0.024 (−88%) | 0.029 (−85%) | 0.217 (+8%) |
| 0.24 | 0.017 (−93%) | 0.010 (−96%) | 0.251 (+5%) |
| **mean recovery** | **18%** | **43%** | **103%** |

Three things in that table matter more than the numbers.

**The ordering reproduces the published disagreement, at the right size.**
Fabre state the clean comparison in their own discussion: DNMT3A at **15.0%**
per year from Watson's spectrum fit against **6.2%** per year from their
follow-up — both gene-level, both DNMT3A. Ratio **2.42**. From design bias
alone this experiment predicts **2.4**.

> An earlier draft used R882H at 0.148 against Fabre's ~0.050, a ratio of 3.0.
> That pair is not like for like — a per-variant hotspot against a gene-level
> average — and Watson report that over 90% of nonsynonymous DNMT3A variants are
> effectively neutral, so the gene average belongs well below the hotspot for
> reasons unrelated to method. Corrected.

The agreement is closer than three seeds and one mutation rate deserve; read it
as consistency, not as a match. **Fabre reach the same place by another route**,
attributing the slowdown to "an increasingly competitive oligoclonal landscape"
— this model's competition term in words.

**The longitudinal bias grows with fitness.** A clone cannot keep growing
exponentially once it owns a large share of the marrow: the cells it competes
against are increasingly its own. The model contains this through the Moran
normalisation, so a large clone's measured growth is genuinely slower than its
fitness. A study that enrols clones big enough to detect, then watches them for
thirteen years, spends much of that window in the saturating regime — and a
fitter clone gets there sooner. At s = 0.24, thirteen years of follow-up
recover almost nothing.

> This is not an error in Fabre's measurement. The growth they report is real.
> It is an error to read that growth as a fitness that would apply to a clone
> starting from one cell.

**The regression on age is not merely biased — it is non-monotonic.** It
returns 0.019 at s = 0.08, rises to 0.033 in the middle, and falls back to
0.017 at s = 0.24. A non-monotonic map has no inverse. One observed slope is
compatible with several very different truths, and no amount of extra data
fixes that.

**But age is not the problem, and the claim stops here.** Watson run an
age-based cross-check of their own: the prevalence of a variant at a fixed
detection threshold, predicted to rise linearly at rate `2Ntμs`. It returns 14%
per year against the 15% their spectrum fit gives. Same sampling design, same
ages, a different functional of them, no detectable bias. What fails is this
particular extraction — the slope of log(VAF) among **detected** clones — and
it fails because it conditions on clearing the threshold. Prevalence counts the
people who carry a variant instead, and escapes.

We cannot borrow it: prevalence needs a screening denominator, and stage C
established this dataset has none. The same missing column that made the
mutation rate unidentifiable blocks the one age-based estimator that works.

---

## 3. Where the bias comes from

One cohort, s = 0.14. Each row switches off **one** artifact a real study
cannot switch off, from the same baseline.

**Regression on host age**

| condition | estimate | error |
|---|---|---|
| as a real study sees it | 0.035 | −75% |
| clone age known | 0.060 | −57% |
| perfect reads (no noise) | 0.034 | −76% |
| deepest real limit (0.0008) | 0.047 | −67% |
| clone age + deep limit | 0.092 | −34% |
| all three off at once | 0.097 | −31% |

**Longitudinal follow-up**

| condition | estimate | error |
|---|---|---|
| as a real study sees it | 0.084 | −40% |
| perfect reads (no noise) | 0.086 | −38% |
| small clones only (VAF < 0.05) | 0.106 | −24% |
| both off at once | 0.110 | −21% |

**The dominant artifact in the regression is that a clone's age is not its
host's age.** A clone in a 70-year-old may be eight years old. The regressor is
therefore the wrong clock, and the error is not noise — it systematically
flattens the line. Handing the model the true clone age nearly doubles the
estimate on its own.

**Read noise contributes essentially nothing** at depth 1,000, in either
design. This was worth testing rather than assuming: selecting clones on a
noisy measurement and then measuring their growth is a textbook route to
regression toward the mean, and it turns out to be negligible here because the
median is robust and the noise at VAF 0.02 with 1,000 reads is small relative
to thirteen years of growth.

**For the longitudinal design the dominant artifact is saturation.**
Restricting to clones below VAF 0.05 recovers a quarter of the missing signal.

> This decomposition works **only** because `s` is held fixed across every clone
> in the simulation. The same stratification on real data means something else:
> size correlates with fitness there, so small clones are also the less fit
> ones. Published work finds micro-CH (VAF 0.5–2%) growing *slower* than CHIP
> (VAF ≥ 2%), the opposite direction — which is what confounding by fitness
> looks like, not evidence against saturation.

**Nothing brings either design all the way back.** Even with every artifact
off, the regression returns 0.097 against a truth of 0.14. The residual is
truncation that cannot be removed: no instrument detects a one-cell clone, and
the clones missing from the young end of the line are exactly the ones that
would make it steep.

---

## 4. The same bias, visible in the real data

Take the six cohorts stage D kept and move **only** the detection floor.

| | slope | n |
|---|---|---|
| no floor (down to VAF 0.0008) | 0.0479 | 887 |
| floor at VAF 0.0192 | 0.0117 | 524 |

The same cohorts, the same variants, a factor of **4.1** between the two
numbers. Nothing about the cells changed — only where the instrument stops
seeing.

The simulation reproduces the sign, and the deep-sequencing end almost exactly:
with a floor of 0.0008 it returned 0.047, against 0.0479 observed. But raising
the floor moves it by a factor of 1.3, not 4.1. **The model does not account
for the whole effect**, and the missing part is identifiable rather than
mysterious:

| | share of the pool |
|---|---|
| Coombs2017, before the floor | 43% |
| Coombs2017, after the floor | 72% |
| Young2019 variants surviving the floor | 158 → 9 |

The simulation is one homogeneous cohort, so raising its floor removes clones
and nothing else. The real pooled regression is six cohorts with different
sensitivities **and** different age distributions. Young2019 sequences deeply
and skews young; Coombs2017 has a limit sitting right at the floor and a median
age of 70. Applying the floor deletes the young, deep-sequenced end of the pool
and leaves a regression run almost within a single cohort. The age lever is
largely gone before the fit starts.

So the 4.1 is two effects stacked: truncation, which the model reproduces at
about 1.3, and a composition shift the model does not contain because it was
never given more than one cohort. Both are artifacts of measurement. Neither is
biology.

---

## 5. What the spectrum-fit result does and does not establish

The grid was produced by this same model, so recovering the input is a check
that the inversion is unbiased and that five quantiles carry enough information
to pin `s` down. **Both could have failed and neither did.**

It is **not** independent evidence that the model describes real
haematopoiesis. No self-consistency check can be. What it rules out is a
specific failure we could not otherwise exclude: that stage D's estimates were
an artifact of the fitting procedure rather than a property of the data.

---

## 6. What this changes

**Stage 1 and stage D were never in conflict.** The regression slope of 0.0479
and the spectrum estimate of ~0.13 are the same population read with two
instruments, one of which does not measure what its units suggest.

**A published `s` is not comparable across designs without a correction.**
Stage D's verdict column — "between the two methods" — treated Watson and Fabre
as two noisy readings of one quantity. They are readings of different
quantities: fitness from birth, and realised growth of an already-large clone
over a specific window at a specific age.

**The comparison to make in future work** is against Watson and other spectrum
fits directly, and against Fabre only after either restricting to small clones
or propagating our model's own saturation through the same thirteen-year
window.

---

## Limitations

- **One model, one truth.** Every number here assumes the simulator's biology:
  constant `s` over life, constant N, no cell death from disease, no
  competition structure beyond the Moran normalisation. A bias measured against
  a wrong model is a bias in the wrong units.
- **No mortality.** People who die are not removed, so the oldest simulated age
  bands contain individuals a real cohort would have lost — plausibly those
  with the largest clones, given the cardiovascular association.
- **Age-independent fitness.** Fabre reports fitness declining with age. If
  real, it would add a bias to the longitudinal design that this experiment
  cannot see, because the model has no such effect to recover.
- **One cohort per simulation.** As section 4 shows, this is exactly what keeps
  the model from reproducing the full real-data effect.

## The two open problems the literature check raised — both now addressed

Worked through in [`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md). In short:

**The third method does not disagree.** Mitchell's per-clone phylogenetic
estimates, which ship with their code, put DNMT3A clades at 0.167–0.200 per year
— above Watson's spectrum fit, not below it. The "~5% per year" recorded in the
literature check is their *inferred underlying spectrum* over every driver that
arises, including the ones that never expand; the measured expanded clades run
10–30%. Comparing a detection-conditioned estimate against the population it was
drawn from was the error, and it is the same error as the R882-versus-gene
-average pair, one layer down. With the comparison corrected, three independent
detection-conditioned methods land between 0.11 and 0.20 and the longitudinal
estimate is the outlier at 0.062. **This is the external corroboration section 5
says this experiment cannot supply for itself.**

**The saturation mechanism was misstated here.** This document says the bias
"grows with `s`". Regressing deceleration on three candidates in a mixed-fitness
cohort says otherwise: the carrier's total fitness-weighted burden growth gives
R² = 0.956, the clone's own size 0.271, and its own fitness **0.003**. A clone
slows by the amount its marrow filled, whoever filled it. Held at fixed size,
deceleration is non-monotonic in fitness — it peaks while the clone is still
expanding and collapses to zero once it has fixed.

That does not fully reconcile with Fabre, and the residue is specific: their
fast drivers show little deceleration *while still growing fast*, and this model
only produces little deceleration once a clone has stopped. Settling it needs
the baseline VAF distribution of their fast-driver clones.

The 43% figure is unaffected — it was measured with one `s` per cohort, where
size and fitness move together as they do in any real detected sample.
