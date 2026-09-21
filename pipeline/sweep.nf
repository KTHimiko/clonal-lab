#!/usr/bin/env nextflow

/*
 * Stage C — parameter sweep for Approximate Bayesian Computation.
 *
 * Each point of the (s, mu) grid is an independent simulation of a whole
 * cohort. Nothing is shared between points, which is what makes this the
 * textbook case for a cluster: the work divides with no coordination.
 *
 * The pipeline deliberately stops at summary statistics. Scoring them
 * against an observed cohort is cheap and happens off the cluster, so the
 * same grid can later be re-scored against a different cohort or a different
 * distance function without re-running a single simulation.
 *
 *   nextflow run pipeline/sweep.nf -profile slurm
 */

// Grids are explicit lists rather than computed ranges. Nextflow 26 rejects
// bare statements at script level and its stricter parser trips over integer
// arithmetic inside closures; listing the values avoids both, and has the
// side benefit of putting the whole grid in plain sight.
// Fitness is linear from 0.02 to 0.40; mutation rate is log-spaced, because
// we have a rough idea of the magnitude of one and not of the other.
params.s_list     = '0.02,0.0471,0.0743,0.1014,0.1286,0.1557,0.1829,0.21,0.2371,0.2643,0.2914,0.3186,0.3457,0.3729,0.4'
params.mu_list    = '5e-07,6.565e-07,8.62e-07,1.132e-06,1.486e-06,1.951e-06,2.562e-06,3.364e-06,4.417e-06,5.8e-06,7.616e-06,1e-05'
params.replicates = 2
params.people     = 5000
params.N          = 50000
params.ages       = '50,60,70,80'
params.limit      = 0.0192
params.outdir     = 'results'
params.code       = "${projectDir}/.."

process SIMULATE {
    tag "s=${String.format('%.3f', s)} mu=${String.format('%.2e', mu)} r=${rep}"

    input:
    tuple val(s), val(mu), val(rep)

    output:
    path "point.csv"

    script:
    """
    python3 ${params.code}/scripts/sweep_point.py \\
        --s ${s} --mu ${mu} --N ${params.N} \\
        --people ${params.people} --ages '${params.ages}' \\
        --limit ${params.limit} --seed ${rep} \\
        --out point.csv
    """
}

process COLLECT {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'part_*.csv'

    output:
    path 'sweep.csv'

    script:
    """
    head -1 \$(ls part_*.csv | head -1) > sweep.csv
    for f in part_*.csv; do tail -n +2 \$f >> sweep.csv; done
    """
}

workflow {
    def s_grid  = params.s_list.toString().split(',').collect  { v -> v.toDouble() }
    def mu_grid = params.mu_list.toString().split(',').collect { v -> v.toDouble() }

    // `as int` is not optional: a parameter given on the command line arrives
    // as a String, and `1.."1"` builds a range over character codes rather
    // than numbers — it silently ran replicates numbered 48 and 49.
    def reps = (1..(params.replicates as int)).toList()

    Channel.fromList([s_grid, mu_grid, reps].combinations())
    | SIMULATE
    | collect
    | COLLECT
}
