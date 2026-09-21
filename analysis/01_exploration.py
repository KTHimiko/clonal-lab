#!/usr/bin/env python3
"""
Exploration of the aggregated dataset from Watson et al. 2020.

Goal: understand what the data can answer BEFORE writing the model. A
simulator calibrated against an observable the data does not measure is
wasted work.

Usage:  .venv/bin/python analysis/01_exploration.py
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
CSV = (ROOT / "reference/data/watson2020/Maximum_likelihood_estimations"
            / "Maximum likelihood estimations - data files"
            / "all_studies_trimmed_all_genes.csv")
FIGS = ROOT / "analysis/figures"
FIGS.mkdir(parents=True, exist_ok=True)

# palette: light surface, recessive ink, at most three series in a scatter
SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "axes.labelsize": 10,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def clean_axes(ax, grid="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(axis=grid, which="major")
    ax.grid(axis=("x" if grid == "y" else "y"), which="major", visible=False)

def section(t):
    print(f"\n{'-'*72}\n{t}\n{'-'*72}")


# ------------------------------------------------------- load and clean
section("1. LOAD AND CLEAN")

raw = pd.read_csv(CSV, lineterminator="\r")
raw.columns = [c.strip() for c in raw.columns]
print(f"rows read: {len(raw)}")

# age arrives as text: 663 rows carry 'noagedata'. Converting with
# errors='coerce' turns the invalid entries into NaN instead of blowing up.
raw["age"] = pd.to_numeric(raw["age"], errors="coerce")
raw["VAF"] = pd.to_numeric(raw["VAF"], errors="coerce")

df = raw.dropna(subset=["age", "VAF"]).copy()
print(f"dropped for missing age: {len(raw) - len(df)}")
print(f"usable: {len(df)}")
print(f"\nage:  {df.age.min():.0f} to {df.age.max():.0f} years "
      f"(median {df.age.median():.0f})")
print(f"VAF:  {df.VAF.min():.4f} to {df.VAF.max():.4f} "
      f"(median {df.VAF.median():.4f})")

# who was dropped matters: if the missing rows come from a single study, the
# loss is systematic, not random.
lost = raw[raw.age.isna()]
if len(lost):
    print("\ncohorts lost to missing age:")
    for study, n in lost.study.value_counts().items():
        tot = (raw.study == study).sum()
        print(f"  {study:<16} {n:>4} of {tot:>4}  ({100*n/tot:.0f}%)")


# ------------------------------------------- the bias that must be faced
section("2. DETECTION LIMIT BY COHORT")
print("Each study sequenced at a different depth, so each one sees clones")
print("only above a different minimum size. Comparing VAF across cohorts")
print("without accounting for this mixes biology with instrument.\n")

lim = (df.groupby("study")
         .agg(n=("VAF", "size"), vaf_min=("VAF", "min"),
              vaf_p05=("VAF", lambda s: s.quantile(0.05)),
              vaf_median=("VAF", "median"), age_median=("age", "median"))
         .sort_values("vaf_min"))
print(lim.to_string(float_format=lambda x: f"{x:.4f}"))


# ------------------------------------------------- the central relation
section("3. VAF AGAINST AGE")
print("Under exponential clone growth, VAF ~ exp(s*t): log(VAF) should rise")
print("linearly with age, and the slope carries the fitness effect s.\n")

x, y = df.age.values, np.log(df.VAF.values)
slope, intercept, r, p, err = stats.linregress(x, y)
print("regression of log(VAF) on age, all genes:")
print(f"  slope:    {slope:+.5f} per year  (standard error {err:.5f})")
print(f"  p-value:  {p:.3e}")
print(f"  R-squared: {r**2:.4f}")
print("\nA low R-squared is expected and informative: age alone does not")
print("determine clone size. Each variant has its own fitness, and chance")
print("dominates while a clone is small. That is exactly what a stochastic")
print("model must reproduce — a spread, not a curve.")

print("\nslope by gene (genes with n >= 40 only):")
rows = []
for gene, g in df.groupby("gene"):
    if len(g) < 40:
        continue
    s, _, rr, pp, ee = stats.linregress(g.age.values, np.log(g.VAF.values))
    rows.append((gene, len(g), s, ee, pp, rr**2))
for gene, n, s, ee, pp, r2 in sorted(rows, key=lambda t: -t[2]):
    sig = "*" if pp < 0.05 else " "
    print(f"  {gene:<8} n={n:>4}  slope {s:+.5f} +/- {ee:.5f}{sig}  R2={r2:.3f}")
print("  (* p < 0.05)")


# ------------------------------------------- the concrete biological question
section("4. DNMT3A R882 AGAINST THE REST OF DNMT3A")
print("R882 is the best-known CHIP hotspot. If it confers a larger advantage,")
print("its clones should be bigger at the same age.\n")

d3 = df[df.gene == "DNMT3A"].copy()
d3["group"] = np.where(d3.variant.str.startswith("R882"), "R882", "other")
for grp, g in d3.groupby("group"):
    print(f"  {grp:<6} n={len(g):>4}  median VAF {g.VAF.median():.4f}  "
          f"median age {g.age.median():.0f}")

a = d3[d3.group == "R882"].VAF.values
b = d3[d3.group == "other"].VAF.values
u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
print(f"\n  Mann-Whitney U: p = {pu:.3e}")
print("  (non-parametric on purpose: VAF distributions are heavily skewed,")
print("   so comparing means would mislead)")


# ------------------------------------------------------------------ figures
section("5. FIGURES")

# Fig 1 - the central relation. Single series: no legend, the title names it.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.scatter(df.age, df.VAF, s=14, alpha=0.45, color=S1,
           edgecolors="none", rasterized=True)
xx = np.linspace(df.age.min(), df.age.max(), 100)
ax.plot(xx, np.exp(intercept + slope*xx), color=INK, linewidth=2)
ax.set_yscale("log")
ax.set_xlabel("age (years)")
ax.set_ylabel("variant allele frequency (VAF)")
ax.set_title("Detected clones grow with age, but spread dominates")
ax.text(0.98, 0.04, f"n = {len(df)} variants\nslope {slope:+.4f}/year   R2 = {r**2:.3f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=INK2)
clean_axes(ax, grid="both")
ax.grid(axis="x", visible=True)
fig.tight_layout(); fig.savefig(FIGS/"01_vaf_by_age.png"); plt.close(fig)
print("  01_vaf_by_age.png")

# Fig 2 - the instrumental bias. Two marks, legend present.
fig, ax = plt.subplots(figsize=(7.2, 4.0))
o = lim.sort_values("vaf_min")
pos = np.arange(len(o))
ax.hlines(pos, o.vaf_min, o.vaf_median, color=BASE, linewidth=2)
ax.scatter(o.vaf_min, pos, s=64, color=S1, zorder=3, label="smallest VAF detected")
ax.scatter(o.vaf_median, pos, s=64, color=S2, zorder=3, label="median VAF")
ax.set_xscale("log")
ax.set_yticks(pos); ax.set_yticklabels([f"{i}  (n={int(n)})" for i, n in zip(o.index, o.n)])
ax.set_xlabel("variant allele frequency (VAF, log scale)")
ax.set_title("Each cohort sees clones above a different size")
ax.legend(loc="lower right")
clean_axes(ax, grid="x")
fig.tight_layout(); fig.savefig(FIGS/"02_detection_limit.png"); plt.close(fig)
print("  02_detection_limit.png")

# Fig 3 - two series: legend present.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for grp, colour in (("R882", S2), ("other", S1)):
    v = np.sort(d3[d3.group == grp].VAF.values)
    ax.step(v, np.arange(1, len(v)+1)/len(v), where="post",
            color=colour, linewidth=2, label=f"DNMT3A {grp} (n={len(v)})")
ax.set_xscale("log")
ax.set_xlabel("variant allele frequency (VAF, log scale)")
ax.set_ylabel("cumulative proportion")
ax.set_title("R882 produces larger clones than other DNMT3A variants"
             if pu < 0.05 else "R882 clones are not larger in this dataset")
ax.text(0.03, 0.95, f"Mann-Whitney p = {pu:.1e}", transform=ax.transAxes,
        va="top", fontsize=9, color=INK2)
ax.legend(loc="lower right")
clean_axes(ax, grid="both")
fig.tight_layout(); fig.savefig(FIGS/"03_dnmt3a_r882.png"); plt.close(fig)
print("  03_dnmt3a_r882.png")

# Fig 4 - small multiples: avoids painting nine cohorts in nine colours.
top = df.gene.value_counts().head(6).index.tolist()
fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.0), sharex=True, sharey=True)
for ax, gene in zip(axes.ravel(), top):
    g = df[df.gene == gene]
    ax.scatter(g.age, g.VAF, s=12, alpha=0.5, color=S1, edgecolors="none")
    if len(g) >= 10:
        s, i2, rr, pp, _ = stats.linregress(g.age.values, np.log(g.VAF.values))
        xs = np.linspace(g.age.min(), g.age.max(), 50)
        ax.plot(xs, np.exp(i2 + s*xs), color=INK, linewidth=1.6)
        ax.set_title(f"{gene}   n={len(g)}   {s:+.4f}/yr", fontsize=10)
    ax.set_yscale("log")
    clean_axes(ax, grid="both")
for ax in axes[1]:
    ax.set_xlabel("age (years)")
for ax in axes[:, 0]:
    ax.set_ylabel("VAF")
fig.suptitle("Each driver gene grows at its own rate",
             x=0.012, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGS/"04_by_gene.png"); plt.close(fig)
print("  04_by_gene.png")


# ------------------------------------------------- what the data allows
section("6. WHAT THIS DATA ALLOWS - AND WHAT IT DOES NOT")
print("""ALLOWS
  - VAF distribution conditional on detection, by age and by gene
  - comparing growth rates between driver genes
  - estimating fitness from the slope of log(VAF) on age

DOES NOT ALLOW
  - population prevalence of CHIP: the file carries only the DETECTED
    variants, with no denominator of how many people were screened per age
    band. Prevalence needs that denominator, and it is not here.
  - individual trajectories: each row is a single measurement of one person.
    Temporal dynamics are inferred by comparing people of different ages, not
    by following the same ones.

CONSEQUENCE FOR THE MODEL
  The observable to reproduce is the DISTRIBUTION of VAF by age, truncated at
  each cohort's detection limit — not a mean curve. The simulator must
  generate realistic spread, and the comparison against observation must
  apply the same truncation, or it compares apples with oranges.""")
