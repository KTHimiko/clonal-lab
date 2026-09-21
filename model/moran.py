"""
Moran process with selection — the minimal model of clonal competition.

THE RULE, at each time step:
  1. one cell is chosen to reproduce, with probability proportional to
     fitness (mutants weigh r = 1+s, wild-type cells weigh 1)
  2. one cell is chosen to die, uniformly at random
  3. the population returns to its original size

Step 1 is selection. Step 2 is drift. Both forces emerge from the same rule,
without being programmed separately.

TWO IMPLEMENTATION DECISIONS THAT DECIDE FEASIBILITY

1. We do not track individual cells. With only two types, the state of the
   system is a single number: `i`, how many cells are mutant. The model
   becomes a birth-death chain — exact, and with no cost proportional to N.

2. We do not run one replicate at a time. All replicates advance together in
   a numpy array, so the number of Python-level operations is the number of
   STEPS, not steps x replicates. In practice this is worth two orders of
   magnitude: the per-replicate loop never finished; the vectorised version
   takes seconds.
"""

import numpy as np


def fixation_probability(N, s, i0=1):
    """
    Probability that the clone takes over the whole population, in closed form.

    Classic result for the Moran process:

                    1 - (1/r)^i0
        P_fix  =  ----------------       with r = 1 + s
                    1 - (1/r)^N

    In the neutral case (s = 0) the limit is simply i0/N: with no advantage,
    the chance of winning is the share you already hold.

    This is the formula the simulator must reproduce.
    """
    if s == 0:
        return i0 / N
    x = 1.0 / (1.0 + s)
    return (1.0 - x**i0) / (1.0 - x**N)


def _vectorised_step(i, N, s, rng):
    """
    Advance one step of the embedded chain, for every active replicate.

    THE EMBEDDED CHAIN
    Most Moran steps change nothing — a mutant replaces another mutant.
    Rather than drawing those wasted steps, we draw only the DIRECTION of the
    next change, conditional on a change occurring:

        P(up | changed) = p_up / (p_up + p_down)

    Identical in distribution, far faster. The price is losing real time,
    which is irrelevant when the question is "does it fix or not".
    """
    r = 1.0 + s
    x = i.astype(np.float64)
    weight = r * x + (N - x)
    p_up = (r * x / weight) * ((N - x) / N)
    p_down = ((N - x) / weight) * (x / N)
    p = p_up / (p_up + p_down)
    return np.where(rng.random(x.size) < p, 1, -1)


def simulate_fixation_probability(N, s, i0=1, replicates=10_000, seed=None,
                                  max_steps=50_000_000):
    """
    Estimate the fixation probability by running many realisations in parallel.

    Returns (estimate, standard_error, steps_to_absorb_all).

    The standard error comes from the binomial — sqrt(p(1-p)/n). It is how
    much the estimate wobbles purely because of the finite replicate count,
    and it decides whether a departure from the formula is real or just
    sampling noise.
    """
    rng = np.random.default_rng(seed)
    i = np.full(replicates, i0, dtype=np.int64)
    active = (i > 0) & (i < N)
    steps = 0

    while active.any():
        idx = np.flatnonzero(active)
        i[idx] += _vectorised_step(i[idx], N, s, rng)
        active = (i > 0) & (i < N)
        steps += 1
        if steps > max_steps:
            raise RuntimeError(f"no absorption after {steps} steps (N={N}, s={s})")

    p = (i >= N).mean()
    return p, np.sqrt(max(p * (1 - p), 1e-12) / replicates), steps


def trajectories(N, s, i0=1, n=40, steps=200_000, seed=None):
    """
    Run `n` realisations with explicit time, recording each path.

    Unlike the function above, every step here is a real Moran step —
    including the ones that change nothing — because drawing a trajectory
    needs real time.

    Returns (times, matrix) with matrix of shape (n, len(times)). Replicates
    already absorbed stay frozen at their final value, which is the correct
    behaviour: extinct stays extinct.

    Time is in generations; N Moran steps is roughly one generation.
    """
    rng = np.random.default_rng(seed)
    r = 1.0 + s
    i = np.full(n, i0, dtype=np.int64)
    history = [i.copy()]

    for _ in range(steps):
        active = (i > 0) & (i < N)
        if not active.any():
            break
        idx = np.flatnonzero(active)
        x = i[idx].astype(np.float64)
        weight = r * x + (N - x)
        p_up = (r * x / weight) * ((N - x) / N)
        p_down = ((N - x) / weight) * (x / N)
        u = rng.random(x.size)
        delta = np.where(u < p_up, 1, np.where(u < p_up + p_down, -1, 0))
        i[idx] += delta
        history.append(i.copy())

    m = np.array(history).T                   # (n, steps)
    t = np.arange(m.shape[1]) / N
    return t, m
