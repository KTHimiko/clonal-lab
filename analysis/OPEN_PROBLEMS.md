# Stage E, part 2 — settling the two open problems

Reproduce with `.venv/bin/python analysis/07_open_problems.py`.
Figure `12_open_problems.png`.

The literature check on stage E
([`METHOD_BIAS_LITERATURE.md`](METHOD_BIAS_LITERATURE.md)) left two problems.
One is resolved. The other is narrowed to a single checkable question, and the
mechanism underneath it turned out to be different from what stage E claimed.

---

## Problem 1 — the phylogenetic method does not disagree

### What the check recorded

Mitchell 2022 reconstructs clonal histories from phylogenies, machinery shared
with neither VAF spectra nor serial sampling. The check recorded them as putting
DNMT3A near 5% per year, siding with longitudinal follow-up against the spectrum
fit — the one result that would have undercut stage E from outside.

### What their own data says

Their per-clone estimates ship with the paper's code, which this project already
holds. 49 expanded clades across 7 donors aged 29–81:

| | s per year |
|---|---|
| median | **0.192** |
| range | 0.106 – 0.355 |
| quartiles | 0.166 / 0.223 |

The nine clades carrying a named DNMT3A driver: median **0.188**, range
0.122–0.200.

### The error was comparing two different quantities

Both numbers are in the paper, and both are correct:

| | what it is | value |
|---|---|---|
| observed expanded clades | measured, the 46 largest clones | **10–30% / yr** |
| underlying gamma spectrum | inferred, over *every* driver that arises, including those that never expand | 5–10% / yr |

The first is conditioned on detection. The second is the distribution those
detections are drawn **from**, and most of what it contains never becomes
visible. Watson's per-variant values and Fabre's followed clones are *both*
conditioned on detection, so comparing either against the underlying spectrum
compares a detected subset against the population it came from — and a detected
subset of a fitness-skewed distribution is fitness-skewed by construction.

This is the same mistake as the R882-versus-gene-average pair from the
literature check, one layer down: two correctly reported numbers, put side by
side when they measure different things.

### The comparison, like for like

| source | method | variant class | s per year |
|---|---|---|---|
| Watson 2020 | VAF spectrum | DNMT3A top-20 variants | 0.112 – 0.160 |
| Mitchell 2022 | phylogenetic | DNMT3A clades (n=9) | 0.167 – 0.200 |
| this project | VAF spectrum | DNMT3A R882 / other | 0.120 – 0.130 |
| Fabre 2022 | longitudinal | DNMT3A, gene level | 0.062 |

**The three detection-conditioned methods sit between 0.11 and 0.20. The
longitudinal one is the outlier, and it is low** — the direction stage E
predicted and roughly the size it predicted.

**This is the independent corroboration stage E could not provide for itself.**
Section 5 of `METHOD_BIAS.md` conceded that recovering our own input from our own
grid checks the inversion and not the biology. The phylogenetic estimate is not
ours, does not use a frequency spectrum, and lands where the spectrum fits land.

### The caveat that does not go away

A phylofit clade enters the analysis *because it expanded*. That is a fitness
filter at least as strong as a VAF threshold. These estimates are the right
comparator for Watson and for us, and the wrong one for Mitchell's own gamma.

---

## Problem 2 — the mechanism was wrong, the conflict is narrowed

### The problem

Stage E wrote that the saturation bias **grows with `s`**. Fabre report the
opposite ordering in real data: deceleration most marked in DNMT3A, TP53 and
BRCC3, almost none in fast drivers like U2AF1 and SRSF2 P95H.

### The experiment

Fabre compare slow and fast drivers **within the same people**, so the
simulation has to as well. The model now draws each clone's fitness from a
heavy-tailed spectrum — `s ∈ {0.05, 0.08, 0.12, 0.20, 0.40}` with weights
`{0.42, 0.27, 0.18, 0.10, 0.03}`, mirroring Mitchell's gamma and Watson's
finding that most nonsynonymous DNMT3A variants are effectively neutral.

Every clone detected at 65 is followed to 78, with its growth rate measured
twice: 65→71.5 and 71.5→78. Deceleration is the first minus the second.

### Deceleration is a property of the host, not the driver

Three candidates, each regressed on deceleration alone:

| predictor | R² |
|---|---|
| **the person's total fitness-weighted burden growth** | **0.956** |
| the clone's own VAF at 65 | 0.271 |
| the clone's own fitness `s` | **0.003** |

A clone slows by almost exactly the amount the fitness-weighted clonal burden of
*its marrow* grew during the window — whether the clone caused that growth or a
neighbour did. That is the Moran normalisation surfacing as the only term that
matters, and it is not what stage E claimed.

### Deceleration is non-monotonic in fitness, once size is held

| VAF at 65 | s=0.08 | s=0.12 | s=0.20 | s=0.40 |
|---|---|---|---|---|
| 0.019–0.04 | 0.006 | 0.014 | 0.042 | **0.183** |
| 0.04–0.08 | — | 0.019 | 0.050 | 0.176 |
| 0.08–0.16 | — | 0.023 | 0.056 | 0.134 |
| 0.16+ | — | — | 0.020 | **0.001** |

In the smallest band the fastest clones decelerate **most** — they are about to
take over the marrow and are generating the burden growth themselves. In the
largest band the same clones decelerate **least** — they already took it over,
and a clone that has finished expanding has nothing left to decelerate.

### What this does to the conflict with Fabre

**It narrows it to one checkable question rather than settling it.**

This model produces "almost no deceleration" for a fast clone caught *after* it
has fixed — the flat tail on the right of the table. It produces the opposite
for a fast clone caught while still expanding. The arithmetic favours the first
state: at `s = 0.4` a surviving clone needs roughly twenty years to reach VAF
0.02 from one cell and about five more to approach fixation, so the window in
which it is both detectable and still growing is short.

**But Fabre report their fast drivers still growing at high rates**, which is
the second state, not the first. Those two facts do not sit together, and
nothing measured here makes them. What would settle it is the baseline VAF
distribution of their fast-driver clones, which the published summaries do not
give.

### The prediction this makes

Testable on an existing longitudinal cohort without new sequencing:
**deceleration should be predicted by the carrier's total clonal burden growth,
not by which gene is mutated.** Match clones on baseline VAF and on carrier
burden, and the gene-level differences should largely disappear. If they do not,
this model's competition term is wrong.

### What stands

The fraction of its own fitness a clone is realising at the moment follow-up
starts, measured directly and for the first time here in a mixed cohort:

| true s | fraction of fitness realised at age 65 |
|---|---|
| 0.08 | 87% |
| 0.12 | 81% |
| 0.20 | 29% |
| 0.40 | ~0% |

Monotonic, severe, and the reason a longitudinal design recovers 43% of the
truth rather than 100%. **Stage E's number survives; its explanation has been
replaced by a better one.**

---

## What changed in the model

`simulate_cohort` now accepts a sequence for `s` with optional `s_weights`, so
each new clone draws its own fitness. A single float still behaves exactly as
before. Snapshots gained `fitness_at`, parallel to `sizes` and `births`, because
with a mixture a recycled slot changes meaning and the final fitness matrix is
not valid for earlier snapshots.

This is a better model regardless of stage E: real marrow contains a spectrum of
drivers competing with each other, and the previous version could only represent
a cohort in which every clone was identical.
