# snakemake_GATKGenotyping
1. Prepare samplesheet. Samplesheet is a list of gvcf files, one file each row.

2. Prepare config files. Modify the samplesheet and workdir in the config.yaml file.

3. Run snakemake:
    A. Run on single node:
        snakemake -pr --keep-going --local-cores 24 --cores 24 --configfile config.yaml --snakefile Snakefile
    B. Run on slurm clusters:
        bash ./submit.sh
