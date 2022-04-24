#!/bin/bash
#SBATCH --job-name=trimgalore_qc-sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-24/trimgalore_qc-sc.log

STAR_dir=/srv/disk00/cheyul1/Venus/STAR
Venus_dir=/srv/disk00/cheyul1/Venus/repo

indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/test_files
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-24/trimgalore_qc-sc

#mkdir -p ${out_dir}/single-end/
#trim_galore \
#--trim-n \
#--quality 5 \
#--phred33 \
#--length 20 \
#--output_dir ${out_dir}/single-end/ \
#--gzip \
#--cores 7 \
#${data_dir}/bulk_1.fastq.gz

#mkdir -p ${out_dir}/paired-end/
#trim_galore \
#--trim-n \
#--quality 5 \
#--phred33 \
#--length 20 \
#--output_dir ${out_dir}/paired-end/ \
#--gzip \
#--cores 7 \
#--paired \
#${data_dir}/bulk_1.fastq.gz ${data_dir}/bulk_2.fastq.gz

mkdir -p ${out_dir}/single-cell/
trim_galore \
--trim-n \
--quality 5 \
--phred33 \
--length 20 \
--output_dir ${out_dir}/single-cell/ \
--gzip \
--cores 8 \
--paired \
${data_dir}/singlecell_1cDNA.fastq.gz ${data_dir}/singlecell_2CB+UMI.fastq.gz
