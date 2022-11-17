__author__ = "Arielle R Munters"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"


rule cp_bam:
    input:
        "alignment/star/{sample}_T.bam",
    output:
        "Results/{sample}_T/{sample}_T.bam",
    params:
        extra=config.get("cp_bam", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.bam.log",
    benchmark:
        repeat("cp_results/{sample}_T.bam.benchmark.tsv", config.get("cp_bam", {}).get("benchmark_repeats", 1))
    threads: config.get("cp_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_bam", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"


rule cp_bai:
    input:
        "alignment/star/{sample}_T.bam.bai",
    output:
        "Results/{sample}_T/{sample}_T.bam.bai",
    params:
        extra=config.get("cp_bai", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.bam.bai.log",
    benchmark:
        repeat( "cp_results/{sample}_T.bam.bai.benchmark.tsv", config.get("cp_bai", {}).get("benchmark_repeats", 1))
    threads: config.get("cp_bai", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_bai", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_bai", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_bai", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_bai", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_bai", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_bai", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"


rule cp_vcf:
    input:
        "snv_indels/pisces/{sample}_T.normalized.sorted.vcf.gz",
    output:
        "Results/{sample}_T/{sample}_T.normalized.sorted.vcf.gz",
    params:
        extra=config.get("cp_vcf", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.normalized.sorted.vcf.gz.log",
    benchmark:
        repeat(
            "cp_results/{sample}_T.normalized.sorted.vcf.gz.benchmark.tsv", config.get("cp_vcf", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cp_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_vcf", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"


rule cp_tbi:
    input:
        "snv_indels/pisces/{sample}_T.normalized.sorted.vcf.gz.tbi",
    output:
        "Results/{sample}_T/{sample}_T.normalized.sorted.vcf.gz.tbi",
    params:
        extra=config.get("cp_tbi", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.normalized.sorted.vcf.gz.tbi.log",
    benchmark:
        repeat(
            "cp_results/{sample}_T.normalized.sorted.vcf.gz.tbi.benchmark.tsv",
            config.get("cp_tbi", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cp_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tbi", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"


rule cp_arrbia_fusions:
    input:
        "fusions/arriba/{sample}_T.fusions.tsv",
    output:
        "Results/{sample}_T/{sample}_T.arriba.fusions.tsv",
    params:
        extra=config.get("cp_arrbia_fusions", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.arriba.fusions.tsv.log",
    benchmark:
        repeat(
            "cp_results/{sample}_T.arriba.fusions.tsv.benchmark.tsv",
            config.get("cp_arrbia_fusions", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cp_arrbia_fusions", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_arrbia_fusions", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_arrbia_fusions", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_arrbia_fusions", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_arrbia_fusions", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_arrbia_fusions", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_arrbia_fusions", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"


rule cp_arrbia_pdf:
    input:
        "fusions/arriba_draw_fusion/{sample}_T.pdf",
    output:
        "Results/{sample}_T/{sample}_T.arriba.pdf",
    params:
        extra=config.get("cp_arrbia_pdf", {}).get("extra", ""),
    log:
        "cp_results/{sample}_T.arriba.pdf.log",
    benchmark:
        repeat(
            "cp_results/{sample}_T.arriba.pdf.benchmark.tsv",
            config.get("cp_arrbia_pdf", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("cp_arrbia_pdf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_arrbia_pdf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_arrbia_pdf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_arrbia_pdf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_arrbia_pdf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_arrbia_pdf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_arrbia_pdf", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    shell:
        "(cp {input} {output}) &> {log}"
