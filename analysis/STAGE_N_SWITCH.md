# Stage N — the fourth mechanism, the closest, and the bracket it closes

Stages L and M rejected three candidates, each by a different argument, and the
three rejections together bounded what was left: the missing mechanism must be
**clone-specific, episodic, and uncorrelated within a host.**

A permanent drop in a clone's own fitness, arriving at a random time, is what a
deleterious second lesion, immune clearance, and lineage exhaustion all look
like in this model. It is the one thing that satisfies all three constraints,
and it escapes stage L's inflation argument because it is **asymmetric** — a
one-way drop cannot compound upward.

**It gets three of four checks and fails the fourth, in the opposite direction
from everything before it.**

---

## 1. The signature that motivated it, found in the data first

If decline is a switch rather than a fluctuation, a clone that has started
falling keeps falling. Measured on 406 Fabre clones with four or more draws:

| | rises next | falls next |
|---|---|---|
| rose first | 183 | 105 |
| **fell first** | 58 | **60** |

`P(fall | fell) = 0.51` against `P(fall | rose) = 0.36`, with a permutation null
of 0.40 ± 0.04 and **P(observed > null) = 0.996**.

**Decline persists.** That is a switch, not a wobble — and unlike the
within-person correlation that killed blood composition, this one clears its
null with room to spare.

---

## 2. What the mechanism reproduces

| rate | drop | d spectrum | declining | d rates | P(f\|f) | P(f\|r) |
|---|---|---|---|---|---|---|
| 0.00 | — | **0.139** | 4.4% | 0.0445 | 0.05 | 0.09 |
| 0.01 | 0.15 | 0.156 | 11.3% | 0.0315 | 0.17 | 0.13 |
| 0.02 | 0.15 | 0.269 | 17.9% | 0.0220 | 0.27 | 0.21 |
| **0.03** | **0.20** | **0.447** | **23.7%** | **0.0147** | **0.42** | **0.22** |
| 0.05 | 0.15 | 0.530 | 36.1% | 0.0240 | 0.44 | 0.32 |

*Targets: declining 25.1%, persistence 0.51 against 0.36, d spectrum ≈ 0.139.*

At a hazard of 0.03 per clone-year with a drop of 0.20:

- **declining 23.7%** against 25.1% observed
- **persistence 0.42 against 0.22** — the gap is 0.20, against 0.15 observed,
  and in the right direction
- **d rates 0.0147** — the best fit this project has produced, against 0.0445
  for the model with no mechanism

Three of four, and the persistence was not calibrated for: it was measured in
the data first and the mechanism was chosen because of it.

---

## 3. Where it fails, and why the failure is structural

**The spectrum degrades 3.2×**, from 0.139 to 0.447, with VAF q90 falling from
0.0404 to 0.0237 against 0.0401 observed. The clones become too small.

The hazard is per year of clone life, so **the clones most likely to have
switched are the ones that have existed longest — which are the largest.** The
mechanism preferentially removes exactly the tail the spectrum needs.

And the trade-off is monotone: every step toward the observed decline is a step
away from the observed sizes. At a hazard of 0.01 the spectrum survives (0.156)
and the decline reaches only 11.3%, half the target.

---

## 4. A correction to how the first attempt was judged

The first run of this test reported **no persistence at all** — 0.21 against
0.22 — and would have rejected the mechanism outright.

That comparison was rigged, by me. The simulation required a clone to be present
at all five draws, while the observed analysis required three or more draws above
the floor. **A clone that shrinks out of sight still contributes to the observed
set, with a short trajectory and a negative slope, and was excluded from the
simulated one** — removing precisely the decliners the mechanism creates.

With the selection matched, persistence appears. The conclusion changed because
the comparison was fixed, not because the model was.

This is the third bug in this session that lived in the comparison rather than
in the model — after a bootstrap variable reused across sections, and an age
vector read from the wrong scope. None of them crashed. **All three returned
plausible numbers**, which is what makes that layer dangerous: a wrong model
usually fails loudly, and a wrong comparison quietly agrees with you.

---

## 5. The bracket, which is the actual result

| mechanism | effect on the size distribution |
|---|---|
| fluctuating fitness | clones **too large** — q90 0.172 against 0.040 |
| niche structure | ceiling `1/(2·niches)` fights the observed tail |
| **switch (this)** | clones **too small** — q90 0.024 against 0.040 |
| blood composition | **size preserved**, killed by ICC 0.674 against 0.129 |

Four mechanisms, and the size distribution is now pinned from both sides.
**Anything strong enough to produce 25% declining clones pushes the spectrum out
of range in one direction or the other — except composition, which preserves it
and is refuted by the within-person correlation instead.**

That stops being a list of failed candidates and becomes a specification. The
missing mechanism must:

> produce decline **without preferentially removing old clones** (or the tail goes);
> **without adding compounding variance** (or the sizes inflate);
> **without acting on a whole person at once** (or the correlation rises).

None of the four satisfies all three. Something that does would be a real
contribution, and this project does not have it.

---

## Limitations

- **One parameter point** for the clone population (mean 0.06, shape 4.0,
  N 200,000) throughout, one seed per row.
- **The hazard is constant per clone-year.** A hazard tied to something other
  than clone age — to the carrier's age, or to the number of divisions — would
  not preferentially remove the largest clones and was not tried. That is the
  most obvious thing to test next and is not tested here.
- **The observed persistence is modest.** 0.51 against 0.36 clears its null, but
  part of it is simply a low-fitness clone being consistently low, which the
  model already produces (0.05 against 0.09 at zero hazard — weak and inverted).
