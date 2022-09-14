import sys
from pysam import VariantFile

vcf_in = VariantFile(snakemake.input.vcf)
vcf_out = VariantFile(snakemake.output.vcf, 'w', header=vcf_in.header)

for record in vcf_in.fetch():
    if record.filter.keys() == ["PASS"] or record.filter.keys() == ["multiallelic"] or record.filter.keys() == ["strand_bias"] or record.filter.keys() == ["haplotype"] or ('strand_bias' in list(record.filter) and 'haplotype' in list(record.filter)):
        vcf_out.write(record)
