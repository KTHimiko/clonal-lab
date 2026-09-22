# Stage L — two mechanisms tested, both fail the same way

Reproduce with the model parameters `niches` and `s_volatility`, both added
here. Each reduces exactly to the previous model at its default, verified
bit-for-bit.

Stage K's corrected diagnosis: **the model couples decline to clone size, and
the data does not.** About a fifth of real clones shrink, at sizes where the
model's clones have nothing to lose to. Three candidates were listed. Two are
now tested.

All runs are fixed at the point that fits the spectrum best — mean 0.06,
shape 4.0, N = 200,000 — so only the new mechanism moves.

Targets: VAF q50 **0.0072**, q90 **0.0401**; declining **25.1%**; width 0.1999.

---

## 1. Niche structure — competition made local

The stem cell pool is divided into compartments; a clone competes only against
the cells in its own. The idea was that a clone could lose a local contest while
the marrow at large is nearly empty.

| niches | VAF q90 | d spectrum | declining | d rates |
|---|---|---|---|---|
| 1 | 0.0360 | **0.147** | 5.8% | 0.0431 |
| 4 | 0.0298 | 0.198 | 7.2% | 0.0380 |
| 12 | **0.0202** | **0.376** | **8.1%** | 0.0368 |

**Rejected.** Declining clones move from 5.8% to 8.1% — two points against the
nineteen that are missing — while the spectrum fit degrades 2.6×.

**Why, and why smaller compartments cannot save it.** With `mu` = 2×10⁻⁶ and
N = 200,000 the influx is 0.4 clones per person per year: about 28 arrivals over
seventy years, most of which die young. Split across 12 compartments that is
roughly **two clones per niche, nearly all microscopic**. Decline needs something
to lose to, and there is nothing there.

Raising the compartment count does not help, because a clone confined to one
cannot exceed VAF `1/(2·niches)` — at 12 compartments the ceiling is 0.042, and
the observed q90 is 0.0401. **The two requirements fight each other by
construction.**

---

## 2. Fluctuating fitness — a clone's advantage varies over life

Each clone's fitness follows an Ornstein–Uhlenbeck excursion around its own
mean, with volatility σ and correlation time τ. A clone can have a bad decade at
any size, whatever its neighbours are doing.

| σ | τ | VAF q90 | d spectrum | declining | d rates |
|---|---|---|---|---|---|
| 0.00 | — | 0.0360 | **0.147** | 5.8% | 0.0431 |
| 0.10 | 2 | 0.0614 | 0.269 | 11.8% | 0.0258 |
| **0.20** | **2** | **0.1724** | **0.856** | **24.0%** | **0.0156** |
| 0.35 | 8 | 0.4939 | 2.335 | 39.0% | 0.0591 |

**It hits the target.** At σ = 0.20 with a two-year correlation time the
declining fraction is 24.0% against 25.1% observed, and the rate distribution is
the **best fit this project has produced** — d = 0.0156, against 0.0431 for the
fixed-fitness model.

**And it destroys the size distribution.** VAF q90 goes from 0.0360 to 0.1724,
four times the observed 0.0401. The spectrum fit degrades **5.8×**.

---

## 3. Why both fail the same way, and what that rules out

Growth is multiplicative: a clone's size goes as the exponential of its
accumulated growth. For a symmetric fluctuation of variance `V` accumulated over
the window,

```
E[size] = exp(mean) × exp(V/2)
```

**Variance inflates the mean.** Measured on the actual parameters:

| σ | τ | accumulated variance | inflation of mean size |
|---|---|---|---|
| 0.10 | 2 yr | 0.44 | ×1.2 |
| 0.20 | 2 yr | 1.76 | ×2.4 |
| 0.20 | 8 yr | 4.21 | ×8.2 |
| 0.35 | 8 yr | 12.89 | ×629 |

And there is no room for it. The observed q90 is 0.0401; the fixed-fitness model
already produces 0.0360. **The size distribution is at its ceiling before any
new variance is added.**

That generalises past the two mechanisms tested:

> **Any mechanism that adds variance to a clone's real growth will inflate clone
> sizes, and the sizes have no headroom. Decline cannot be bought that way.**

Niche structure fails for the same reason in disguise: to make local competition
bite, clones must be large relative to their compartment, which means either
large clones or a low ceiling — and both break the spectrum.

---

## 4. What survives

| candidate | status |
|---|---|
| ~~niche structure~~ | **tested, rejected** — 2 points of decline, 2.6× worse spectrum |
| ~~fluctuating fitness~~ | **tested, rejected** — hits the decline target, 5.8× worse spectrum |
| **blood composition** | **untested, and now the only one that can work** |

The third candidate is the only one that does not touch real growth. VAF measures
a clone's share of *circulating cells*, not of stem cells. If lineage output
shifts — myeloid against lymphoid, with age or with illness — the measured VAF
moves while the clone does not.

That produces apparent decline **without adding any variance to growth**, which
section 3 shows is the binding constraint. It is the only one of the three that
is not ruled out by the argument that rules out the other two.

---

## Limitations

- **One parameter point** (mean 0.06, shape 4.0, N 200,000) for every run. The
  mechanisms were tested where the spectrum fits best, not re-optimised jointly.
  A joint search might find a compensating region, though section 3's argument
  suggests the ceiling is structural rather than local.
- **One seed per cell** of both tables.
- **The OU process is one choice** among many for how fitness might vary. A
  process with occasional large drops rather than symmetric excursions would add
  less variance for the same decline, and was not tried.
