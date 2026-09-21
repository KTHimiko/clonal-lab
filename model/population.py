"""
Multi-clone Moran process with a continuous influx of mutations.

Stage A simulated one clone against the wild type, where an exact formula
exists. This module drops that restriction: clones arise over time, coexist,
and compete with each other as well as with the wild type. No closed form
survives — which is the whole reason for simulating.

WHAT IT PRODUCES
For each simulated person, the sizes of every clone alive at a given age.
Divided by 2N, those become variant allele frequencies, directly comparable
with the observed data.

THREE DESIGN DECISIONS

1. TAU-LEAPING INSTEAD OF SINGLE STEPS
   A single Moran step changes the population by one cell. With N = 50,000
   stem cells over 80 years that is four million steps per person, and we
   need thousands of people. Instead we advance by a time slice `dt` and draw
   how many births and deaths happened in it from a Poisson distribution.
   Exact in the limit of small dt, and thousands of times faster.

2. A PEOPLE-BY-CLONES MATRIX
   Every simulated person advances in the same numpy operation. The state is
   a matrix: one row per person, one column per clone slot. Empty slots hold
   size zero and cost nothing.

3. COMPETITION IS EXPLICIT, NOT ASSUMED AWAY
   It is tempting to let each clone grow independently, which is accurate
   while clones are small. We keep the full Moran normalisation instead, so
   the model stays correct when a clone becomes large — exactly the regime
   that matters after age 70.
"""

import numpy as np


def simulate_cohort(N=50_000, s=0.10, mu=2e-6, years=80.0, people=1_000,
                    dt=0.1, max_clones=64, seed=None, record_ages=None,
                    initial_clones=0, s_weights=None):
    """
    Simulate a cohort of people from birth to `years`.

    Parameters
    ----------
    N : int
        Number of haematopoietic stem cells. Published estimates span
        50,000 to 200,000.
    s : float or sequence of float
        Selection coefficient of a driver mutation, per year. A sequence makes
        every new clone draw its own value from that set, which is what a real
        marrow looks like: one person carries a slow DNMT3A clone and a fast
        splicing-gene clone at the same time, and they compete with each other.
        A single float keeps every clone identical, which is what isolating one
        effect at a time requires.
    s_weights : sequence of float, optional
        Probabilities for the values in `s`. Uniform when omitted.
    mu : float
        Driver mutation rate per cell per year. With N cells, new clones
        appear at rate N*mu per year.
    years : float
        Age to simulate up to.
    people : int
        How many independent people to simulate.
    dt : float
        Time slice, in years. Smaller is more accurate and slower.
    max_clones : int
        Slots per person. A person who exceeds this is rare; the excess is
        counted and reported rather than silently dropped.
    record_ages : sequence of float, optional
        Ages at which to snapshot clone sizes. Defaults to the final age.
    initial_clones : int
        Seed every person with this many clones of size 1 at birth. Used for
        validation: with `mu=0` and one seeded clone the model reduces to the
        single-clone case of stage A, where the exact answer is known.

    Returns
    -------
    dict with
        'ages'      the ages snapshotted
        'sizes'     list of (people, max_clones) arrays, one per age
        'births'    list of (people, max_clones) arrays, one per age, holding
                    the age at which each live clone arose (NaN for empty
                    slots). Two things need it. First, a clone's age is not
                    its host's age, and telling those apart is the whole
                    subject of stage E. Second, slots are recycled: a clone
                    that dies frees its slot for a different clone, so
                    following 'the clone in slot 7' across two snapshots
                    without checking the birth time silently splices two
                    different clones into one fake trajectory.
        'fitness_at' list of (people, max_clones) arrays, one per age, holding
                    each live clone's own s. With a mixture of fitnesses a slot
                    changes meaning when it is recycled, so the final matrix is
                    not valid for earlier snapshots.
        'fitness'   (people, max_clones) array of each slot's s at the END
        'overflow'  how many mutation events found no free slot
    """
    rng = np.random.default_rng(seed)
    record_ages = np.atleast_1d(record_ages if record_ages is not None else years)

    s_values = np.atleast_1d(np.asarray(s, dtype=np.float64))
    if s_weights is not None:
        s_weights = np.asarray(s_weights, dtype=np.float64)
        s_weights = s_weights / s_weights.sum()

    def draw_s(size):
        if s_values.size == 1:
            return np.full(size, s_values[0])
        return rng.choice(s_values, size=size, p=s_weights)

    sizes = np.zeros((people, max_clones), dtype=np.float64)
    fitness = np.zeros((people, max_clones), dtype=np.float64)
    births = np.full((people, max_clones), np.nan, dtype=np.float64)
    overflow = 0

    if initial_clones:
        sizes[:, :initial_clones] = 1.0
        fitness[:, :initial_clones] = draw_s((people, initial_clones))
        births[:, :initial_clones] = 0.0

    snapshots, birth_snaps, fitness_snaps = [], [], []
    snap_at = list(np.sort(record_ages))
    steps = int(round(years / dt))
    new_clone_rate = N * mu * dt          # expected new clones per person per slice

    for step in range(1, steps + 1):
        t = step * dt

        # --- competition: the Moran weights over the whole population -------
        # W is the total reproductive weight. Wild-type cells weigh 1 each,
        # clone cells weigh (1+s) each. A clone's share of W is its chance of
        # being the one that reproduces.
        mutant_total = sizes.sum(axis=1, keepdims=True)
        wild = N - mutant_total
        W = (sizes * (1.0 + fitness)).sum(axis=1, keepdims=True) + wild

        # --- exact birth-death transition over the slice -------------------
        # A clone of c cells is c independent lineages. Over dt years, with
        # per-cell birth rate lam and death rate 1, the linear birth-death
        # process has a KNOWN transition distribution — we sample from it
        # directly instead of approximating with Poisson counts.
        #
        # WHY NOT TAU-LEAPING
        # Drawing births and deaths from the size at the start of the slice
        # lets a one-cell clone die and give birth in the same slice, so it
        # resurrects. That inflated survival by 3 standard errors at dt=0.1.
        # Applying deaths first fixes the resurrection and creates the
        # opposite error, ten times larger, by killing clones that would have
        # reproduced before dying. Neither ordering is right, because in the
        # real process the events interleave.
        #
        # THE EXACT FORM
        # Each lineage goes extinct with probability alpha; each survivor
        # leaves a geometrically distributed number of cells. So the size
        # after the slice is:
        #     survivors  ~ Binomial(c, 1 - alpha)
        #     size       = survivors + NegativeBinomial(survivors, 1 - beta)
        # which is exact, and needs no small dt to be correct.
        lam = N * (1.0 + fitness) / W          # per-cell birth rate, per year
        r = lam - 1.0                          # net growth rate
        E = np.exp(np.clip(r * dt, -50, 50))

        den = lam * E - 1.0
        near_neutral = np.abs(r) < 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha = np.where(near_neutral, lam * dt / (1.0 + lam * dt),
                             (E - 1.0) / den)
            beta = np.where(near_neutral, alpha, lam * (E - 1.0) / den)
        alpha = np.clip(np.nan_to_num(alpha, nan=1.0), 0.0, 1.0 - 1e-12)
        beta = np.clip(np.nan_to_num(beta, nan=0.0), 0.0, 1.0 - 1e-12)

        alive = sizes > 0
        counts = sizes.astype(np.int64)
        survivors = np.zeros_like(counts)
        survivors[alive] = rng.binomial(counts[alive], (1.0 - alpha)[alive])

        grew = survivors > 0
        extra = np.zeros_like(counts)
        extra[grew] = rng.negative_binomial(survivors[grew], (1.0 - beta)[grew])
        sizes = (survivors + extra).astype(np.float64)
        np.clip(sizes, 0, N, out=sizes)

        # a clone that hit zero is gone: free its slot
        extinct = sizes <= 0
        sizes[extinct] = 0.0
        fitness[extinct] = 0.0
        births[extinct] = np.nan

        # --- new mutations --------------------------------------------------
        n_new = rng.poisson(new_clone_rate, size=people)
        for k in np.flatnonzero(n_new):
            free = np.flatnonzero(sizes[k] == 0)
            take = min(n_new[k], free.size)
            overflow += int(n_new[k] - take)
            if take:
                slots = free[:take]
                sizes[k, slots] = 1.0
                fitness[k, slots] = draw_s(slots.size)
                births[k, slots] = t

        # --- snapshot --------------------------------------------------------
        while snap_at and t >= snap_at[0] - 1e-9:
            snapshots.append(sizes.copy())
            birth_snaps.append(births.copy())
            fitness_snaps.append(fitness.copy())
            snap_at.pop(0)

    while snap_at:                                   # ages beyond the run
        snapshots.append(sizes.copy())
        birth_snaps.append(births.copy())
        fitness_snaps.append(fitness.copy())
        snap_at.pop(0)

    return {"ages": np.sort(record_ages), "sizes": snapshots,
            "births": birth_snaps, "fitness_at": fitness_snaps,
            "fitness": fitness, "overflow": overflow}


def sizes_to_vaf(sizes, N):
    """
    Convert clone sizes to variant allele frequency.

    A clone of `c` cells out of `N` stem cells makes up c/N of the blood. The
    mutation sits on one of two chromosome copies, so the fraction of mutated
    reads is half that:

        VAF = c / (2N)

    This factor of two is the most common source of confusion in the field.
    """
    return sizes / (2.0 * N)


def detected(vaf, limit):
    """
    Keep only the clones an instrument with this detection limit would see.

    Every real cohort has one, and it varies 300-fold across the studies in
    our data. Comparing a simulation against observation without applying the
    same truncation penalises the model for predicting clones nobody could
    have measured.
    """
    v = vaf[vaf >= limit]
    return v[v > 0]
