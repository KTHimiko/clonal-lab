# Stage M — the third mechanism fits the target and is refuted by the data

Stage L left one candidate standing for the decline the model cannot produce:
**blood composition.** VAF measures a clone's share of *circulating* cells, not
of stem cells, so a shift in lineage output moves the measurement without moving
the clone.

It was the only candidate not ruled out by stage L's argument, because it does
not add variance to real growth — the effect perturbs the measurement rather
than compounding through it.

**It works, and then it dies on a statistic it also predicts.**

---

## 1. It fits, where nothing else did

Implemented as a per-person Ornstein–Uhlenbeck factor in log space, applied to
every clone in that person at that draw, mean-corrected so it shifts spread and
trend but not the average.

| mechanism | d spectrum | declining | d rates |
|---|---|---|---|
| none | 0.147 | 5.8% | 0.0431 |
| niche structure (12) | 0.376 | 8.1% | 0.0368 |
| fluctuating fitness (σ 0.20, τ 2) | **0.856** | 24.0% | 0.0156 |
| **composition (σ 0.75, τ 8)** | **0.141** | **21.2%** | **0.0163** |

Declining clones reach 21.2% against 25.1% observed. The rate fit ties the best
this project has produced. **And the spectrum is untouched** — 0.141 against
0.139 with no mechanism at all, with VAF q90 at 0.0412 against 0.0401 observed.

Exactly as predicted: because the effect does not compound, it buys decline
without paying in clone size.

---

## 2. The strength it needs is enormous

At σ = 0.75 the factor drifts by a standard deviation of **0.072 per year** over
a thirteen-year window. On its own that produces a rate spread of 0.185, against
the 0.200 actually observed.

**For this mechanism to explain the decline, nearly the entire observed spread
of growth rates must be measurement artifact rather than biology.**

That is a strong claim, and it makes a falsifiable prediction: if a person's
clones are all being moved by one shared factor, **they must move together**.

---

## 3. The prediction, tested

Intraclass correlation of growth rates among clones within the same person.

| σ | τ | ICC predicted |
|---|---|---|
| 0.00 | 8 | −0.020 |
| 0.30 | 8 | 0.209 |
| 0.50 | 8 | 0.530 |
| **0.75** | **8** | **0.674** |
| 0.75 | 50 | 0.369 |

**Observed in Fabre's data: 0.129** (396 clones in 145 people; permutation null
0.040 ± 0.059).

**The strength needed to explain the decline predicts an ICC of 0.674. The data
says 0.129 — off by a factor of five.**

Rejected. The mechanism reproduces the statistics it was fitted to and fails the
independent one it also predicts, which is the only kind of test worth running.

---

## 4. What the same test does establish

Read the table the other way. At σ = 0 the model predicts ICC = −0.020; the data
says 0.129, about 1.5 null standard deviations above. Interpolating, the observed
correlation is consistent with **σ ≈ 0.2**.

So there is evidence for a composition effect. It is simply far too small to be
the answer: at σ = 0.2 the declining fraction is around 7%, against 25% observed.

**Two statements, both true:** a modest composition effect is supported by the
within-person correlation, and it explains almost none of the decline.

And it is worth noting what the plain model does *not* explain: it produces an
ICC of −0.02 while the data shows 0.129. Our model contains clonal interference,
and in this regime that interference produces no within-person correlation at
all — clones are too small to interact. **The observed correlation is real and
unexplained by the model as it stands.**

---

## 5. All three candidates are now closed

| mechanism | why it fails |
|---|---|
| niche structure | 2 points of decline; a clone confined to a compartment cannot exceed VAF `1/(2·niches)`, so the ceiling and the requirement fight |
| fluctuating fitness | hits the decline target; variance compounds multiplicatively and inflates clone sizes 2.4× to 629× |
| blood composition | hits the decline target and preserves sizes; predicts ICC 0.674 against 0.129 observed |

Each was rejected by a different argument, which is the useful part. Together
they bound what the missing mechanism can be:

> It must **not add variance to real growth** (or sizes inflate).
> It must **not require nearby competitors** (clones are too small and too few).
> It must **not act on all of a person's clones at once** (or the ICC rises).

That leaves something **clone-specific, episodic, and uncorrelated within a
host**: a clone acquiring a deleterious second lesion, immune clearance of
particular clones, differentiation or exhaustion of a single lineage. None of
these is in the model, and none was tested here.

---

## Limitations

- **One parameter point** for the clone population (mean 0.06, shape 4.0,
  N 200,000), and one seed per cell of every table.
- **The ICC comparison is not confounder-free.** Clonal interference also
  predicts within-person correlation. In this model it produces none, which is
  why the comparison is usable here — but in a model where clones were larger it
  would not be.
- **The OU process is one shape** of composition drift. A process with occasional
  large shifts rather than smooth excursions would produce a different ICC for
  the same decline, and was not tried.
