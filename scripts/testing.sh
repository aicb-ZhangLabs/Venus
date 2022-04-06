#!/bin/bash
#SBATCH --job-name=testing_bulkSingleMEGAV
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-06/testing_bulkSingleMEGAV.log

STAR_dir=/srv/disk00/cheyul1/Venus/STAR
Venus_dir=/srv/disk00/cheyul1/Venus/repo

indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/YaChi_GSE112576
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-06/testing_bulkSingleMEGAV

python3 ${Venus_dir}/module-detection.py \
    --read ${data_dir}/SRR6944349.1_1.fastq.gz \
    --virusGenome ${indices_dir}/new_virus.genomeDir \
    --humanGenome ${indices_dir}/human.genomeDir \
    --out ${out_dir}/ \
    --readFilesCommand zcat \
    --thread 32
