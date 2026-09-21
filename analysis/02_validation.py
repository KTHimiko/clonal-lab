#!/usr/bin/env python3
"""
Stage A — validate the simulator against the analytical expectation.

A simulator that fails to reproduce the case where a formula exists does not
deserve trust where none does. This is the only stage of the project where we
know the right answer in advance, and therefore the only chance to check.

Usage:  .venv/bin/python analysis/02_validation.py
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from model.moran import (fixation_probability, simulate_fixation_probability,
                         trajectories)

FIGS = ROOT / "analysis/figures"
FIGS.mkdir(parents=True, exist_ok=True)

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.size": 10,
    "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.titlepad": 14, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.grid": True, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 9, "legend.labelcolor": INK2,
    "figure.dpi": 130,
})

def clean(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

def section(t):
    print(f"\n{'-'*74}\n{t}\n{'-'*74}")


# =================================================== test 1: the neutral case
section("TEST 1 - NEUTRAL CASE (s = 0)")
print("With no advantage at all, a clone's chance of taking over must equal")
print("the share it already holds: i0/N. It is the simplest case and the one")
print("that catches the crudest implementation bugs.\n")

print(f"{'N':>6}  {'i0':>4}  {'expected':>10}  {'simulated':>10}  {'+/-':>8}  {'sigmas':>8}")
failures = 0
for N, i0 in [(100, 1), (100, 10), (100, 50), (500, 1), (500, 25)]:
    exp = fixation_probability(N, 0.0, i0)
    sim, err, _ = simulate_fixation_probability(N, 0.0, i0, replicates=10_000, seed=42)
    z = abs(sim - exp) / err if err > 0 else 0
    mark = "ok" if z < 3 else "FAILED"
    failures += z >= 3
    print(f"{N:>6}  {i0:>4}  {exp:>10.5f}  {sim:>10.5f}  {err:>8.5f}  {z:>6.2f}  {mark}")

print("\n'sigmas' is the gap between simulated and expected measured in standard")
print("errors. Below 3 it is indistinguishable from sampling noise, which is")
print("what a correct implementation should produce.")


# ============================================== test 2: with a fitness effect
section("TEST 2 - WITH A SELECTIVE ADVANTAGE")
print("Now the full formula:  P_fix = (1 - (1/r)^i0) / (1 - (1/r)^N),")
print("with r = 1+s. This is where selection and drift genuinely compete.\n")

N = 200
print(f"N = {N}, starting from a single mutant cell\n")
print(f"{'s':>7}  {'expected':>10}  {'simulated':>10}  {'+/-':>8}  {'sigmas':>8}")
svals = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
sim_vals, err_vals = [], []
for s in svals:
    exp = fixation_probability(N, s, 1)
    sim, err, _ = simulate_fixation_probability(N, s, 1, replicates=10_000, seed=7)
    sim_vals.append(sim); err_vals.append(err)
    z = abs(sim - exp) / err if err > 0 else 0
    mark = "ok" if z < 3 else "FAILED"
    failures += z >= 3
    print(f"{s:>7.2f}  {exp:>10.5f}  {sim:>10.5f}  {err:>8.5f}  {z:>6.2f}  {mark}")


# ================================= test 3: the rule of thumb, actually checked
section("TEST 3 - THE LARGE-POPULATION APPROXIMATION")
print("For large N and small s the formula collapses into something memorable.")
print("Worth checking which value is right, because two of them circulate.\n")

N = 100_000
print(f"N = {N:,}\n")
print(f"{'s':>7}  {'exact':>10}  {'s/(1+s)':>10}  {'2s':>10}")
for s in [0.001, 0.005, 0.01, 0.05, 0.1]:
    exact = fixation_probability(N, s, 1)
    print(f"{s:>7.3f}  {exact:>10.6f}  {s/(1+s):>10.6f}  {2*s:>10.6f}")

print("\nThe correct approximation for the MORAN process is s/(1+s), i.e. ~s.")
print("The '2s' rule found in many texts comes from the Wright-Fisher model,")
print("which has a different offspring variance. Both models describe the same")
print("biology and disagree by a factor of two — using the wrong one would")
print("double our fitness estimates.")


# ===================================================================== figures
section("FIGURES")

# Fig 1 - simulated against analytical. Two series: legend present.
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ss = np.linspace(0, 0.55, 200)
ax.plot(ss, [fixation_probability(200, s, 1) for s in ss],
        color=S1, linewidth=2, label="exact formula")
ax.errorbar(svals, sim_vals, yerr=np.array(err_vals)*1.96, fmt="o",
            color=S2, markersize=7, capsize=4, linewidth=1.6,
            label="simulation (10k replicates, 95% CI)")
ax.set_xlabel("selection coefficient  s")
ax.set_ylabel("fixation probability")
ax.set_title("The simulator reproduces the exact formula")
ax.text(0.97, 0.06, "N = 200, starting from 1 mutant cell",
        transform=ax.transAxes, ha="right", fontsize=9, color=INK2)
ax.legend(loc="upper left")
clean(ax)
fig.tight_layout(); fig.savefig(FIGS/"05_fixation_validation.png"); plt.close(fig)
print("  05_fixation_validation.png")

# Fig 2 - trajectories: why one simulation says nothing.
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
for ax, s, title in [(axes[0], 0.0, "no advantage (s = 0)"),
                     (axes[1], 0.1, "with advantage (s = 0.1)")]:
    N = 200
    t, m = trajectories(N, s, i0=1, n=40, steps=150_000, seed=1000)
    fixed = int((m[:, -1] >= N).sum())
    for k in range(m.shape[0]):
        won = m[k, -1] >= N
        ax.plot(t, m[k]/N, color=(S2 if won else S1),
                alpha=(0.85 if won else 0.30),
                linewidth=(1.4 if won else 1.0))
    ax.set_title(f"{title}   —   {fixed} of 40 fixed", fontsize=10)
    ax.set_xlabel("time (generations)")
    clean(ax)
axes[0].set_ylabel("fraction of the population")
fig.suptitle("Same parameters, different fates: drift decides",
             x=0.012, ha="left", fontsize=12, fontweight="bold", color=INK)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIGS/"06_trajectories.png"); plt.close(fig)
print("  06_trajectories.png")


# ===================================================================== verdict
section("VERDICT")
if failures == 0:
    print("All tests passed within 3 sigma.")
    print("\nThe simulator reproduces the exact formula in the neutral case and")
    print("under selection, across several population sizes and starting")
    print("points. It is fit to be used where no formula exists.")
else:
    print(f"WARNING: {failures} test(s) beyond 3 sigma. Do not proceed without investigating.")
    sys.exit(1)
