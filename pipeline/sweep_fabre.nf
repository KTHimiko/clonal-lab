#!/usr/bin/env nextflow

/*
 * Stage J — a grid shaped to the Fabre cohort, emitting two summary vectors.
 *
 * Stage C's grid was built for the Watson table and stage F had to score a
 * different cohort against it, which is the limitation this closes. Every point
 * here uses Fabre's baseline age spread, their 13-year follow-up, their VAF
 * floor and their sequencing depth, and records BOTH the cross-sectional
 * spectrum and the per-clone growth rates from the same simulated people.
 *
 * Decision rule fixed in advance: analysis/PREREGISTRATION_STAGE_J.md
 *
 *   nextflow run pipeline/sweep_fabre.nf -profile slurm
 */

// Explicit lists, not computed ranges: Nextflow 26 rejects bare statements at
// script level and trips over integer arithmetic inside closures. Listing the
// values also puts the whole grid in plain sight.
params.s_list     = '0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30'
params.mu_list    = '5e-07,9.5e-07,1.8e-06,3.4e-06,6.4e-06,1.2e-05'
params.replicates = 3
params.people     = 8000
params.N          = 100000
params.outdir     = 'results_fabre'
params.code       = "${projectDir}/.."

process SIMULATE {
    tag "s=${String.format('%.2f', s)} mu=${String.format('%.1e', mu)} r=${rep}"

    // Measured, not guessed: the worst grid point peaks at 336 MB, and since
    // this session's cgroup work the limit is enforced rather than advisory.
    // The profile default of 512 MB would leave no headroom, and an OOM
    // halfway through 252 jobs costs more than the reserved memory.
    memory '768 MB'
    cpus 1
    time '20m'

    input:
    tuple val(s), val(mu), val(rep)

    output:
    path "point.csv"

    script:
    """
    python3 ${params.code}/scripts/sweep_fabre.py \\
        --s ${s} --mu ${mu} --N ${params.N} \\
        --people ${params.people} --seed ${rep} \\
        --out point.csv
    """
}

process COLLECT {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'part_*.csv'

    output:
    path 'sweep_fabre.csv'

    script:
    """
    head -1 \$(ls part_*.csv | head -1) > sweep_fabre.csv
    for f in part_*.csv; do tail -n +2 \$f >> sweep_fabre.csv; done
    """
}

workflow {
    def s_grid  = params.s_list.toString().split(',').collect  { v -> v.toDouble() }
    def mu_grid = params.mu_list.toString().split(',').collect { v -> v.toDouble() }

    // `as int` is not optional: a parameter from the command line arrives as a
    // String, and 1.."3" builds a range over character codes. It silently ran
    // replicates numbered 49 and 50 the first time this was written.
    def reps = (1..(params.replicates as int)).toList()

    Channel.fromList([s_grid, mu_grid, reps].combinations())
    | SIMULATE
    | collect
    | COLLECT
}
