# Pre-registration — stage I, the Lothian cohort

**Written and committed before the Lothian growth rates were computed.** The
design facts below (sample sizes, depths, wave structure) were inspected first,
because they are properties of the study rather than of its outcome. No rate,
no comparison and no figure from this cohort has been produced at the time of
writing.

The reason for the formality: stage H tested three related measures and reported
the one that reached the edge of significance. That was disclosed, but
disclosure does not undo it. Fixing the analysis in advance does.

---

## 1. What stage H found, and how much to trust it

234 Fabre trajectories, each clone split at its own midpoint:

| gene | late-half rate | 95% interval |
|---|---|---|
| DNMT3A | 0.0096 | −0.0174 to +0.0284 |
| TET2 | 0.0401 | +0.0133 to +0.0598 |

Difference **+0.0305**, one-sided P = 0.968, two-sided interval grazing zero.

**Its power, computed here for the first time, was 0.56.** Resampling stage H's
own empirical distributions at its own sample sizes, a one-sided test at α = 0.05
detects that difference only 56% of the time.

That changes how the stage H result should be read. A borderline positive from a
study with 56% power is the profile of a finding that does not replicate — not
evidence that it is wrong, but a reason to hold it loosely. **This is recorded
as a retrospective correction to stage H's own write-up.**

---

## 2. Why this cohort cannot settle it

Robertson 2022, Lothian Birth Cohorts, GEO **GSE178936**, open access.

| | |
|---|---|
| participants | 85 (LBC1921 42, LBC1936 43) |
| waves per clone | mostly 3 or 4 |
| wave spacing | LBC1921 ≈ 3 yr, LBC1936 ≈ 4 yr |
| depth | median 1,127× overall; 1,968× for TET2 |
| usable clones (≥3 waves) | **TET2 49, DNMT3A 129** |
| LBC1921 only, entirely past 75 | TET2 23, DNMT3A 69 |

Power for the stage H test, from stage H's own spread at these sample sizes:

| configuration | power |
|---|---|
| stage H itself (114/120) | 0.56 |
| Lothian ≥3 waves (49/129) | 0.35 |
| the same, with the extra noise of two-point halves | **0.25** |
| LBC1921 alone (23/69) | lower still |

**A significance test here would be uninformative in both directions.** A null
carries a 75% chance of missing a real effect. A positive, at this power, is more
likely to be noise than confirmation.

---

## 3. What will be done instead, fixed in advance

**Not a hypothesis test. An estimate with an interval.**

A small study cannot declare significance, but it can contribute an effect
estimate that is either consistent or inconsistent with +0.0305. That is the
useful thing it can do, and it is the only thing that will be claimed.

**The single pre-specified quantity:** the difference in median late-half growth
rate, TET2 minus DNMT3A, in clones with at least three waves.

**Fixed analysis choices**, identical to stage H except where the data forces a
change, and each change named:

| choice | value | why |
|---|---|---|
| VAF floor | 0.002 | as stage H |
| minimum timepoints | **3** (stage H used 4) | only 56 clones have 4, and they exist in one cohort only |
| split | at each clone's own median wave, middle point shared | as stage H |
| rate fit | weighted least squares on log(VAF), weights from read depth | as stage H |
| time axis | wave × 3 years (LBC1921), wave × 4 years (LBC1936) | ages are in dbGAP under controlled access and will not be requested |
| interval | 4,000-sample bootstrap of the median difference | as stage H |

**What will be reported:** the point estimate, its interval, and whether stage
H's +0.0305 falls inside it.

**What will not be claimed:** that the effect is confirmed, or that it is
refuted. Neither conclusion is available at this power, and reaching for either
would be the same error stage H's three-measure disclosure was meant to avoid.

**A prediction, so this is falsifiable:** if stage H's effect is real and the
same size, the Lothian point estimate should be positive, and stage H's +0.0305
should fall within the Lothian interval. If the point estimate is negative, that
is evidence against — weak evidence, but evidence.

---

## 4. The known weaknesses of this test, listed before seeing it

- **Ages are approximated from wave number.** A constant scaling error inflates
  or deflates both genes' rates together, so the *sign* of the difference is
  robust and its *magnitude* is not.
- **Two-point halves.** With three waves, each half is two measurements, so each
  per-clone rate is a single difference. The aggregate median is still
  meaningful; nothing per-clone will be reported.
- **A different gene panel.** 75 genes by targeted capture, against Fabre's
  panel. Which variants are called is not identical.
- **A different population.** Scottish birth cohorts against a Sardinian founder
  population.
- **Older on entry.** Most of the Lothian window sits past the age where stage H
  says the genes diverge, so there is less contrast between the halves than in
  Fabre — which reduces the expected effect size on top of reducing power.
