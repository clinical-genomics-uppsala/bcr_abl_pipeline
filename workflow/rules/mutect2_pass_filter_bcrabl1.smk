__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2022, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutect2_pass_filter_bcrabl1:
    input:
        vcf="snv_indels/gatk_mutect2/{sample}_{type}.merged.softfiltered.vcf.gz",
    output:
        vcf=temp("snv_indels/gatk_mutect2/{sample}_{type}.bcrabl1.merged.vcf.gz"),
    params:
        extra=config.get("mutect2_pass_filter_bcrabl1", {}).get("extra", ""),
    log:
        "bcr_abl_pipeline/mutect2_pass_filter_bcrabl1/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "bcr_abl_pipeline/mutect2_pass_filter_bcrabl1/{sample}_{type}.output.benchmark.tsv",
            config.get("mutect2_pass_filter_bcrabl1", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mutect2_pass_filter_bcrabl1", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mutect2_pass_filter_bcrabl1", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutect2_pass_filter_bcrabl1", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutect2_pass_filter_bcrabl1", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mutect2_pass_filter_bcrabl1", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutect2_pass_filter_bcrabl1", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mutect2_pass_filter_bcrabl1", {}).get("container", config["default_container"])
    conda:
        "../envs/mutect2_pass_filter_bcrabl1.yaml"
    message:
        "{rule}: Hardfilter all but PASS, multiallelic or strand_bias variants from {input.vcf} to {output.vcf}"
    script:
        "../scripts/mutect2_pass_filter_bcrabl1.py"
