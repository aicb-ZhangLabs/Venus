#!/bin/bash
#SBATCH --job-name=testing_bulkPairedHIV
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-07/testing_bulkPairedHIV.log

STAR_dir=/srv/disk00/cheyul1/Venus/STAR
Venus_dir=/srv/disk00/cheyul1/Venus/repo

indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/YaChi_GSE112576
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-07/testing_bulkPairedHIV

# Testing Detection Module
python3 ${Venus_dir}/module-detection.py \
    --read ${data_dir}/SRR6944349.1_1.fastq.gz ${data_dir}/SRR6944349.1_2.fastq.gz \
    --virusGenome ${out_dir}/new_virus.genomeDir \
    --humanGenome ${out_dir}/human.genomeDir \
    --out ${out_dir} \
    --readFilesCommand zcat \
    --thread 32


## Testing Indexing Module
#python3 ${Venus_dir}/module-index.py \
#    --humanGenome ${out_dir}/human.genomeDir \
#    --humanFASTA ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.fna \
#    --humanGTF ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.gtf \
#    --virusGenome ${out_dir}/virus.genomeDir \
#    --virusFASTA ${indices_dir}/new_virus.genomeDir/mega-virus.fa \
#    --out ${out_dir} \
#    --thread 32
