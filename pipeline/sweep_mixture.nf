#!/usr/bin/env nextflow

/*
 * Stage K — a fitness DISTRIBUTION, with N free.
 *
 * Decision rule and identifiability checks fixed in advance:
 * analysis/PREREGISTRATION_STAGE_K.md
 *
 * Batched deliberately. Stage J ran 252 jobs of 7.5 seconds and spent 36
 * minutes on 3.5 CPU-hours, because submitting and polling a job cost more
 * than running it. Each job here sweeps a whole row of means.
 *
 *   nextflow run pipeline/sweep_mixture.nf -profile slurm
 */

params.mean_list  = '0.04,0.06,0.08,0.10,0.12,0.14,0.16'
params.shape_list = '0.5,1.0,2.0,4.0'
params.n_list     = '50000,100000,200000,400000'
params.replicates = 2
params.people     = 6000
params.mu         = 2e-6
params.outdir     = 'results_mixture'
params.code       = "${projectDir}/.."

process SIMULATE {
    tag "N=${N} shape=${shape} r=${rep}"

    // Measured on the worst point of the previous sweep at 336 MB; since this
    // session's cgroup work the limit is enforced, not advisory.
    memory '768 MB'
    cpus 1
    time '30m'

    input:
    tuple val(N), val(shape), val(rep)

    output:
    path "batch.csv"

    script:
    """
    python3 ${params.code}/scripts/sweep_mixture.py \\
        --mean-list '${params.mean_list}' \\
        --shape ${shape} --N ${N} --mu ${params.mu} \\
        --people ${params.people} --seed ${rep} \\
        --out batch.csv
    """
}

process COLLECT {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'part_*.csv'

    output:
    path 'sweep_mixture.csv'

    script:
    """
    head -1 \$(ls part_*.csv | head -1) > sweep_mixture.csv
    for f in part_*.csv; do tail -n +2 \$f >> sweep_mixture.csv; done
    """
}

workflow {
    def n_grid     = params.n_list.toString().split(',').collect     { v -> v.toInteger() }
    def shape_grid = params.shape_list.toString().split(',').collect { v -> v.toDouble() }

    // `as int`: a parameter from the command line arrives as a String, and
    // 1.."2" builds a range over character codes.
    def reps = (1..(params.replicates as int)).toList()

    Channel.fromList([n_grid, shape_grid, reps].combinations())
    | SIMULATE
    | collect
    | COLLECT
}
