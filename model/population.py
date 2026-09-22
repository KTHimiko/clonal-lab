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
                    initial_clones=0, s_weights=None,
                    s_ramp=None, wt_decline=0.0, age_effects_from=50.0,
                    niches=1, s_volatility=0.0, s_tau=5.0,
                    switch_rate=0.0, switch_delta=0.0):
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
    s_ramp : float or sequence of float, optional
        Fractional gain in a clone's OWN fitness per year of host age beyond
        `age_effects_from`. One value per entry of `s` when a sequence. This is
        the "the mutation gets better" mechanism.
    wt_decline : float
        Exponential decay rate of the WILD-TYPE reproductive weight per year of
        host age beyond `age_effects_from`. This is the "everyone else gets
        worse" mechanism, and it is not the same thing: fitness in a Moran
        process is relative, so degrading the wild type lifts every clone at
        once, while `s_ramp` lifts only the classes it is given for. Mouse work
        attributes the TET2 age effect mainly to this one — aged non-mutant
        stem cells losing ground rather than the mutant gaining it.
    age_effects_from : float
        Host age at which both effects start. Before it nothing changes, so a
        run with either parameter still reduces exactly to the constant model
        over the early decades.
    switch_rate : float
        Hazard, per clone per year, of a permanent drop in that clone's own
        fitness. Zero reproduces the previous model exactly.

        Stage M closed the last of three candidate mechanisms and the three
        rejections together bound what is left: it cannot add variance to real
        growth (sizes inflate), cannot need nearby competitors (clones are too
        small), and cannot act on all of one person's clones at once (the
        within-person correlation rises). A clone-specific, one-way, episodic
        event satisfies all three, and it is what a deleterious second lesion,
        immune clearance, or lineage exhaustion all look like in this model.

        It escapes stage L's inflation argument because it is ASYMMETRIC. The
        variance that inflated clone sizes came from symmetric fluctuation,
        whose upswings compound. A drop that only goes down cannot inflate
        anything.
    switch_delta : float
        How much fitness falls when a clone switches. Subtracted, so a clone at
        s = 0.06 with a delta of 0.10 becomes s = -0.04 and shrinks.
    s_volatility : float
        Standard deviation of a clone's fitness fluctuation around its own
        mean. Zero reproduces the fixed-fitness model exactly.

        Stage K found the data demanding decline that is NOT a consequence of
        size: about a fifth of real clones shrink, at sizes where the model's
        clones have nothing to lose to. Niche structure was tested and moved
        that fraction by two points while breaking the size distribution. A
        fitness that varies over time decouples the two directly — a clone can
        have a bad decade at any size, whatever its neighbours are doing.
    s_tau : float
        Correlation time of that fluctuation, in years. Short means the
        fluctuation averages out within a follow-up window and changes little;
        long means a clone can spend a whole study period disadvantaged.
    niches : int
        Number of independent compartments the stem cell pool is divided into.
        Competition is resolved WITHIN a compartment: a clone born in one
        competes only against the N/niches cells there.

        With `niches=1` this is exactly the global Moran process and every
        earlier result is unchanged. Above 1 it decouples decline from global
        size, which is what stage K found the data demanding — a clone can lose
        its local contest while the marrow at large is nearly empty.

        The ceiling it imposes is real and worth knowing: a clone confined to
        one compartment cannot exceed VAF 1/(2*niches), so the observed q90 of
        0.04 needs `niches` no larger than about 12.
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

    ramp_values = (np.zeros_like(s_values) if s_ramp is None
                   else np.broadcast_to(np.atleast_1d(
                       np.asarray(s_ramp, dtype=np.float64)),
                       s_values.shape).copy())

    def draw_s(size):
        """Draw a fitness and the ramp that belongs to it, as a pair."""
        if s_values.size == 1:
            idx = np.zeros(size, dtype=np.int64)
        else:
            idx = rng.choice(s_values.size, size=size, p=s_weights)
        return s_values[idx], ramp_values[idx]

    sizes = np.zeros((people, max_clones), dtype=np.float64)
    fitness = np.zeros((people, max_clones), dtype=np.float64)
    # which compartment each clone was born into; -1 marks an empty slot
    niche = np.full((people, max_clones), -1, dtype=np.int64)
    # the current fitness excursion of each clone, an Ornstein-Uhlenbeck process
    excursion = np.zeros((people, max_clones), dtype=np.float64)
    rho = float(np.exp(-dt / s_tau)) if s_volatility > 0 else 0.0
    switched = np.zeros((people, max_clones), dtype=bool)
    p_switch = 1.0 - float(np.exp(-switch_rate * dt)) if switch_rate > 0 else 0.0
    niche_capacity = N / float(niches)
    person_of = np.repeat(np.arange(people), max_clones).reshape(people, max_clones)
    births = np.full((people, max_clones), np.nan, dtype=np.float64)
    ramp = np.zeros((people, max_clones), dtype=np.float64)
    overflow = 0

    if initial_clones:
        sizes[:, :initial_clones] = 1.0
        fitness[:, :initial_clones], ramp[:, :initial_clones] = \
            draw_s((people, initial_clones))
        births[:, :initial_clones] = 0.0
        niche[:, :initial_clones] = rng.integers(0, niches,
                                                 (people, initial_clones))

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
        # --- age effects -----------------------------------------------------
        # Both act only after `age_effects_from`, and they are distinguishable:
        # `s_ramp` raises the fitness of the classes it was given for, while
        # `wt_decline` lowers the wild type and therefore lifts EVERY clone,
        # including ones with no ramp at all. That difference is the whole
        # reason for implementing them separately.
        # A clone that switches keeps the lower fitness for good, so the drop
        # is applied to `fitness` itself rather than to the per-step value.
        if p_switch > 0:
            eligible = (sizes > 0) & ~switched
            hit = eligible & (rng.random(sizes.shape) < p_switch)
            if hit.any():
                fitness[hit] -= switch_delta
                switched[hit] = True

        elapsed = max(t - age_effects_from, 0.0)
        if s_volatility > 0:
            # z_new = rho*z + sqrt(1-rho^2)*sigma*noise keeps the marginal
            # standard deviation at sigma whatever dt is, so the amount of
            # fluctuation does not depend on the integration step.
            excursion = (rho * excursion
                         + np.sqrt(1.0 - rho * rho) * s_volatility
                         * rng.standard_normal(excursion.shape))
            fitness_now = fitness * (1.0 + ramp * elapsed) + excursion
            # a clone cannot have a negative birth rate
            np.clip(fitness_now, -0.9, None, out=fitness_now)
        else:
            fitness_now = fitness * (1.0 + ramp * elapsed)
        wt_weight = np.exp(-wt_decline * elapsed)

        if niches == 1:
            mutant_total = sizes.sum(axis=1, keepdims=True)
            wild = N - mutant_total
            W = ((sizes * (1.0 + fitness_now)).sum(axis=1, keepdims=True)
                 + wild * wt_weight)
            W_here = W
            capacity = N
        else:
            # Totals per (person, compartment), accumulated into a flat array so
            # the whole cohort is reduced in one pass.
            live = niche >= 0
            flat = person_of * niches + np.where(live, niche, 0)
            cells = np.zeros(people * niches)
            weight = np.zeros(people * niches)
            np.add.at(cells, flat[live], sizes[live])
            np.add.at(weight, flat[live], (sizes * (1.0 + fitness_now))[live])
            wild_n = np.maximum(niche_capacity - cells, 0.0)
            W_flat = weight + wild_n * wt_weight
            W_here = np.where(live, W_flat[flat], 1.0)
            capacity = niche_capacity

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
        lam = capacity * (1.0 + fitness_now) / W_here   # per-cell birth rate
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
        np.clip(sizes, 0, capacity, out=sizes)

        # a clone that hit zero is gone: free its slot
        extinct = sizes <= 0
        sizes[extinct] = 0.0
        fitness[extinct] = 0.0
        ramp[extinct] = 0.0
        niche[extinct] = -1
        excursion[extinct] = 0.0
        switched[extinct] = False
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
                fitness[k, slots], ramp[k, slots] = draw_s(slots.size)
                births[k, slots] = t
                niche[k, slots] = rng.integers(0, niches, slots.size)

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
