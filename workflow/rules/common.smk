__author__ = "Arielle R. Munters"
__copyright__ = "Copyright 2022, Arielle R. Munters"
__email__ = "arielle.munters@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

min_version("6.8.0")

### Set and validate config file


configfile: "config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


def generate_read_group_star(wildcards):
    return "--outSAMattrRGline  ID:{} SM:{} PL:Illumina PU:{} LB:{}".format(
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
        "{}_{}".format(wildcards.sample, wildcards.type),
    )


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def compile_output_list(wildcards):
    output_files = ["qc/multiqc/multiqc_RNA.html"]
    output_files.append(
        [
            "alignment/star/%s_%s.bam%s" % (sample, type, suffix)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
            for suffix in ["", ".bai"]
        ]
    )
    output_files.append(
        [
            "fusions/arriba/%s_%s.fusions.tsv" % (sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "fusions/arriba_draw_fusion/%s_%s.pdf" % (sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "fusions/star_fusion/%s_%s/star-fusion.fusion_predictions.tsv" % (sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    output_files.append(
        [
            "snv_indels/pisces/%s_%s.normalized.sorted.vcf.gz%s" % (sample, type, suffix)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
            for suffix in ["", ".tbi"]
        ]
    )
    output_files.append(
        [
            "Results/%s_%s_summary.xlsx" % (sample, type)
            for sample in get_samples(samples)
            for type in get_unit_types(units, sample)
        ]
    )
    return output_files
