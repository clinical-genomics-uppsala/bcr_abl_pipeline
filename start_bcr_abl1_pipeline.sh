#!/bin/bash

#borde gora en specifik for abl pipeline!
source /projects/wp4/nobackup/workspace/arielle_test/twist/hydra/Workarea/220222_Test/venv/bin/activate
module load singularity/3.4.2 slurm-drmaa/1.1.1
echo "Module loaded"
seqrun=$1
snakemake_profile=/projects/wp2/nobackup/CGU_2022_12_BCR-ABL1/Bin/bcr_abl_pipeline/snakemake_profile_config.yaml

outbox_dir=
start_dir=$(pwd)
bin_dir=/projects/wp2/nobackup/CGU_2022_12_BCR-ABL1/Bin/bcr_abl_pipeline

hydra-genetics create-input-files -d ${start_dir}/fastq/ -p MiSeq \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -t R -s "(R[0-9]{2}-[0-9]{5})" -b NNNNNNNN && \
echo "Hydra genetics create input files done" && \
cp ${bin_dir}/config/*yaml ./ && \
echo "Cp configfiles to ${start_dir} done" && \
# Cp sample.tsv units.tsv resources.yaml config.yaml to scratch
mkdir -p /scratch/wp2/abl/${seqrun}/ && \
rsync -ruvp *tsv /scratch/wp2/abl/${seqrun}/ && \
rsync -ruvp *yaml /scratch/wp2/abl/${seqrun}/ && \

cd /scratch/wp2/${seqrun}/ && \
echo "Cp files to scratch and move to scratch done" && \
snakemake --profile ${snakemake_profile} --configfile config.yaml && \

echo "Snakemake done" && \
mkdir -p ${outbox_dir}/${seqrun} && \
rsync -ruvp -r Results/ ${outbox_dir}/${seqrun} && \
touch ${outbox_dir}/${seqrun}/Done.txt && \
echo "Cp to outbox done" && \
rsync -ruvp -r Results ${start_dir}/ && \
echo "Done!"
