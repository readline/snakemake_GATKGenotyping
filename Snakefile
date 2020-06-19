# Version 1.0
# 200618
# Kai Yu
from os.path import join
import pandas as pd
import yaml

configfile: 'config.yaml'
print('='*80)
print('Analysis configs:')
for i in config:
    print('%s:\t%s'%(i,config[i]))
print('='*80)


with open(config['samplesheet']) as infile:
    gvcflist = [i.strip() for i in infile.readlines()]

workdir: config['workdir']

    
rule all:
    input:
        "Genotype/Merge.flt.vqsr.vcf.anno/Merge.Anno.matrix.gz",
    threads:  1
    resources:
        mem  = 1*1024
    shell:
        "rm -rf VQSR"

        
rule GenomicsDBImport:
    input:
        bed="/data/yuk5/pipeline/wxs_phase1/ref/IDT_xGen_Exome_Research_Panel/hg38/intervals_withflank150_collapse/hg38.xGen.{itv}.bed",
        gvcf=lambda wildcards: ["{}".format(gvcf) for gvcf in gvcflist ]
    output:
        itvvcf=temp("Genotype/VQSR/Merge.itv_{itv}.vcf.gz"),
        itvmf=temp("Genotype/VQSR/Merge.itv_{itv}.mf.vcf.gz"),
        itvflt="Genotype/VQSR/Merge.itv_{itv}.flt.vcf.gz"
    log:
        out="logs/A1.GenomicsDBImport.{itv}.log",
        err="logs/A1.GenomicsDBImport.{itv}.err"
    message: "Running GenomicsDBImport of interval {wildcards.itv}."
    threads:  4
    resources:
        mem  = 32*1024
    run:
        inputs = " ".join("-V {}".format(f) for f in input.gvcf)
        shell(
        """
        module load {config[modules][gatk]}
        
        gatk --java-options "-Xmx4g -Xms4g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenomicsDBImport \
            {inputs} \
            -L {input.bed} \
            --genomicsdb-workspace-path /lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv}
            
        gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GenotypeGVCFs \
            -R {config[references][fasta]} \
            -O {output.itvvcf} \
            -D {config[references][gatk_dbsnp]} \
            -G StandardAnnotation \
            --only-output-calls-starting-in-intervals \
            --use-new-qual-calculator \
            -V gendb:///lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv} \
            -L {input.bed} 
            
        rm -rf /lscratch/$SLURM_JOB_ID/gdb.itv_{wildcards.itv}
        
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantFiltration \
            --filter-expression "ExcessHet>54.69" \
            --filter-name ExcessHet \
            -V {output.itvvcf} \
            -O {output.itvmf}
            
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            SelectVariants \
            --exclude-filtered \
            -V {output.itvmf} \
            -O {output.itvflt}
        """
        )
        
        
rule Genotyping:
    input:
        itvs=lambda wildcards: \
             ["Genotype/VQSR/Merge.itv_{}.flt.vcf.gz".format(itv) for itv in range(1, config['references']['interval']+1)]
    output:
        mvcf1="Genotype/VQSR/Merge.flt.vcf.gz",
        sidrecal="Genotype/VQSR/Merge.flt.sid.recal",
        snvrecal="Genotype/VQSR/Merge.flt.snv.recal",
        sidtranch="Genotype/VQSR/Merge.flt.sid.tranches",
        snvtranch="Genotype/VQSR/Merge.flt.snv.tranches",
        indrecal="Genotype/VQSR/Merge.flt.tmp.indel.recalibrated.vcf.gz",
        vqsr="Genotype/Merge.flt.vqsr.vcf.gz"
    log:
        out="logs/A2.Genotyping.log",
        err="logs/A2.Genotyping.err"
    threads:  4
    resources:
        mem  = 50*1024
    message: "Running Genotyping."
    run:
        inputs = " ".join("--INPUT {}".format(i) for i in input.itvs)
        shell(
        """
        module load {config[modules][gatk]} {config[modules][samtools]}
        gatk --java-options "-Xmx3g -Xms3g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            GatherVcfs \
            {inputs} \
            --OUTPUT {output.mvcf1}
            
        tabix -p vcf {output.mvcf1}
        
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantRecalibrator \
            -R {config[references][fasta]} \
            -V {output.mvcf1} \
            -O {output.sidrecal} \
            --tranches-file {output.sidtranch} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 \
            -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
            -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \
            -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            --max-gaussians 4 \
            -resource mills,known=false,training=true,truth=true,prior=12:{config[references][gatk_1000g]} \
            -resource axiomPoly,known=false,training=true,truth=false,prior=10:{config[references][gatk_axiom]} \
            -resource dbsnp,known=true,training=false,truth=false,prior=2:{config[references][gatk_dbsnp]}

        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            VariantRecalibrator \
            -V {output.mvcf1} \
            -O {output.snvrecal} \
            --tranches-file {output.snvtranch} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 \
            -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource hapmap,known=false,training=true,truth=true,prior=15:{config[references][gatk_hapmap]} \
            -resource omni,known=false,training=true,truth=true,prior=12:{config[references][gatk_omni]} \
            -resource 1000G,known=false,training=true,truth=false,prior=10:{config[references][gatk_1000hc]} \
            -resource dbsnp,known=true,training=false,truth=false,prior=7:{config[references][gatk_dbsnp]}
            
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            ApplyVQSR \
            -V {output.mvcf1} \
            -O {output.indrecal} \
            --recal-file {output.sidrecal} \
            --tranches-file {output.sidtranch} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode INDEL
        gatk --java-options "-Xmx24g -Xms24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
            ApplyVQSR \
            -V {output.indrecal} \
            -O {output.vqsr} \
            --recal-file {output.snvrecal} \
            --tranches-file {output.snvtranch} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode SNP
        """
        )
        
        
rule Annotation:
    input:
        "Genotype/Merge.flt.vqsr.vcf.gz"
    output:
        result="Genotype/Merge.flt.vqsr.vcf.anno/Merge.Anno.matrix.gz"
    log:
        out="logs/A3.Annotation.log",
        err="logs/A3.Annotation.err"
    threads:  48
    resources:
        mem  = 100*1024
    message: "Running Annotation."
    shell:
        """
        {config[bins][vcfanno]} \
            {input} \
            `dirname {output.result}` \
            {threads} n
        """
