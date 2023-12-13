"""Rossi et al., 2019"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("Run", drop=False)


##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand(["piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx", "piscem_spliceu/{run}/{run}.h5ad"],
                run=samples["Run"]),
        expand("scrublet/{run}/{run}_initial_annotation_nc.h5ad",
                run=samples["Run"]),
        "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
        "output/tables/01-eda-whole_dataset-nc/parameters.json",
        "output/tables/01-eda-whole_dataset-nc/rossi2019-LHA-dropseq_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
        "output/tables/01-eda-whole_dataset-nc/rossi2019-LHA-dropseq_all_mrk-logreg_sct-combined-whole_dataset-nc.csv"

##### load rules #####

rule map_spliceu:
    input:
        r1="fastq/{run}_1.fastq.gz",
        r2="fastq/{run}_2.fastq.gz"
    output:
        map="piscem_spliceu/{run}/piscem_map/map.rad"
    params:
        prefix="piscem_spliceu/{run}/piscem_map",
        act="piscem_spliceu/{run}/.afhome"
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 32
    resources:
        mem_mb=32000
    benchmark:
        "benchmarks/piscem_spliceu_map/{run}.tsv"
    shell:
        ("export ALEVIN_FRY_HOME={params.act} \
        && simpleaf set-paths \
        && piscem map-sc \
        --index /data/GRCm39/index/piscem_idx \
        --threads {threads} \
        -o {params.prefix} \
        -1 {input.r1} \
        -2 {input.r2} \
        --geometry '1{{b[12]u[8]x:}}2{{r:}}' ")

rule quant_spliceu:
    input:
        map="piscem_spliceu/{run}/piscem_map/map.rad"
    output:
        result="piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx"
    params:
        prefix="piscem_spliceu/{run}/",
        map="piscem_spliceu/{run}/piscem_map",
        act="piscem_spliceu/{run}/.afhome",
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run]
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 32
    resources:
        mem_mb=32000
    benchmark:
        "benchmarks/piscem_spliceu/{run}.tsv"
    shell:
        ("export ALEVIN_FRY_HOME={params.act} \
        && simpleaf set-paths \
        && simpleaf quant \
        -c '1{{b[12]u[8]x:}}2{{r:}}' \
        -o {params.prefix} \
        -t {threads} \
        --map-dir {params.map} \
        -r cr-like \
        --knee \
        -m /data/GRCm39/index/t2g_3col.tsv ")

rule get_h5ad:
    input:
        gene="piscem_spliceu/{run}/af_quant/alevin/quants_mat.mtx"
    output:
        knee="output/figures/{run}_raw/knee-plot.pdf",
        h5ad="piscem_spliceu/{run}/{run}.h5ad"
    params:
        sample_run_name="{run}",
        expected_num_cells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        path="piscem_spliceu/{run}/af_quant"
    container:
        "docker://combinelab/usefulaf:0.9.3"
    threads: 4
    resources:
        mem_mb=16000
    script:
        "../code/pyroe.py"


rule doublets_call:
    input:
        filt_h5ad="piscem_spliceu/{run}/{run}.h5ad"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls_nc.tsv",
        h5ad="scrublet/{run}/{run}_initial_annotation_nc.h5ad"
    params:
        expected_dblt=0.1,
        sample_run_name="{run}",
        plots="output/figures/{run}_raw/doublets_call_nc"
    container:
        "docker://etretiakov/scrna-seq:jammy-2022.12.09-v0.0.1"
    threads: 8
    script:
        "../code/scrublet_cb.py"


rule exploratory_data_analysis:
    input:
        rmd="analysis/01-eda-whole_dataset-nc.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.001.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-nc.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-nc.pdf",
        "output/tables/01-eda-whole_dataset-nc/rossi2019-LHA-dropseq_all_mrk-MAST_sct-combined-whole_dataset-nc.csv",
        "output/tables/01-eda-whole_dataset-nc/rossi2019-LHA-dropseq_all_mrk-logreg_sct-combined-whole_dataset-nc.csv",
        "output/tables/01-eda-whole_dataset-nc/parameters.json"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 32
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")
