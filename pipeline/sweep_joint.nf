#!/usr/bin/env nextflow

/*
 * Stage O — the clone population and the hazard, swept together.
 *
 * Every earlier run of the hazard fixed mean, shape and N at values fitted with
 * no mechanism at all, which understates it: a hazard that removes clones
 * changes which population fits best.
 *
 * Plan and decision rule: analysis/PREREGISTRATION_STAGE_O.md
 *
 *   nextflow run pipeline/sweep_joint.nf -profile slurm
 */

params.mean_list  = '0.05,0.07,0.09,0.11'
params.shape_list = '1.0,2.0,4.0'
params.n_list     = '100000,200000,400000'
params.rate_list  = '0.0,0.02,0.03,0.05'
params.from_list  = '55,60'
params.replicates = 2
params.people     = 7000
params.outdir     = 'results_joint'
params.code       = "${projectDir}/.."

process SIMULATE {
    tag "N=${N} sh=${shape} r=${rate} f=${frm} rep=${rep}"

    // 336 MB measured on the largest earlier point; enforced since this
    // session's cgroup work, so the headroom is deliberate.
    memory '768 MB'
    cpus 1
    time '30m'

    input:
    tuple val(N), val(shape), val(rate), val(frm), val(rep)

    output:
    path "batch.csv"

    script:
    """
    python3 ${params.code}/scripts/sweep_joint.py \\
        --mean-list '${params.mean_list}' \\
        --shape ${shape} --N ${N} --rate ${rate} --from-age ${frm} \\
        --people ${params.people} --seed ${rep} \\
        --out batch.csv
    """
}

process COLLECT {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'part_*.csv'

    output:
    path 'sweep_joint.csv'

    script:
    """
    head -1 \$(ls part_*.csv | head -1) > sweep_joint.csv
    for f in part_*.csv; do tail -n +2 \$f >> sweep_joint.csv; done
    """
}

workflow {
    def n_grid     = params.n_list.toString().split(',').collect     { v -> v.toInteger() }
    def shape_grid = params.shape_list.toString().split(',').collect { v -> v.toDouble() }
    def rate_grid  = params.rate_list.toString().split(',').collect  { v -> v.toDouble() }
    def from_grid  = params.from_list.toString().split(',').collect  { v -> v.toDouble() }

    // `as int`: a command-line parameter arrives as a String, and 1.."2" builds
    // a range over character codes.
    def reps = (1..(params.replicates as int)).toList()

    Channel.fromList([n_grid, shape_grid, rate_grid, from_grid, reps].combinations())
    | SIMULATE
    | collect
    | COLLECT
}
