#!/bin/bash
module load snakemake samtools || exit 1

sbcmd="sbatch --cpus-per-task={threads} --mem={resources.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={log.out} --error={log.err} {cluster.extra}"

snakemake -pr --keep-going --local-cores 4 \
        --jobs 50 --cluster-config cluster.yaml --cluster "$sbcmd" \
            --latency-wait 120 all --configfile config.yaml --snakefile Snakefile
