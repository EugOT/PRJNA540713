"""individual"""
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from pandas import read_table
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.5.1")

##### load config and sample sheets #####
configfile: "config.yaml"

samples = read_table(config["samples"]).set_index(["run"], drop=False)

proj_dir_path=abspath("../")
refdata_dir_path=abspath("/mnt/pub/GENCODE/M27/")
whitelists_dir_path=abspath("/mnt/pub/whitelists/")

##### target rules #####

shell.executable("/bin/bash")
shell.prefix("source /home/etretiakov/.bashrc; export LC_ALL=en_US.utf-8 && export LANG=en_US.utf-8 &&")

rule all:
    input:
        expand(["kb_{sample}/counts_filtered/adata.h5ad"],
                sample=samples["run"])

# localrules: all, kallisto_g
ruleorder: kb_gc

##### load rules #####

rule kb_gc:
    input:
        r1=join("fastq_gz","{sample}_1.fastq.gz"),
        r2=join("fastq_gz","{sample}_2.fastq.gz")
    output:
        "kb_{sample}/counts_filtered/adata.h5ad"
    log:
        "logs/kallisto-bus_{sample}.log"
    params:
        idx=join(refdata_dir_path, "transcriptome.idx"),
        t2g=join(refdata_dir_path, "transcripts_to_genes.txt"),
        tech=lambda wildcards: str(samples["technology"][wildcards.sample]).upper(),
        outdir="kb_{sample}"
    conda: "kallisto-bus.yaml"
    threads: 10
    shell:
        ("kb count \
            -i {params.idx} \
            -g {params.t2g} \
            -o {params.outdir} \
            -x {params.tech} \
            -t {threads} \
            --h5ad \
            --filter bustools \
            {input.r1} {input.r2}")
