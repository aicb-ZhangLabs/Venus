#!/bin/bash
#SBATCH --job-name=sameFQ
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/srv/disk00/cheyul1/logs/10-21-21/sameFQ.log
#SBATCH --partition=zhanglab.p

Venus_dir=/srv/disk00/cheyul1/Venus
out_dir=/srv/disk00/cheyul1/outputs/10-21-21

#for i in {2..70}; do
#    for file in ${out_dir}/fastqs/GSE126230-run$i/*.fastq.gz; do
#        mv "$file" "$(echo "$file" | sed s/_/_R/)"
#    done
#done
#
#for i in {0..39}; do
#    for file in ${out_dir}/fastqs/GSE126230-run$i/*.fastq.gz; do
#        mv "$file" "$(echo "$file" | sed s/_/_R/)"
#    done
#done

for i in {35..35}; do
    cd ${out_dir}
    mkdir GSE126230-run$i
    cd GSE126230-run$i
    
    bash ${Venus_dir}/chimFinder/cheyu-wrapper.sh \
    /srv/disk00/cheyul1/outputs/09-24-21/fastqs/GSE126230-run$i \
    GSE126230-run$i \
    ${Venus_dir}/bowtie2/indices/HIV/HIV \
    ${Venus_dir}/bowtie2/indices/human/human \
    ${Venus_dir}/STAR/indices/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    64 \
    ${Venus_dir}/hisat2/human.ss \
    ${Venus_dir}/hisat2/HIV_splicesites.txt \
    ${Venus_dir}/bowtie2/indices/HXB2/HXB2 \
    ${Venus_dir}/bowtie2/indices/HXB2/HXB2.fa \
    trim_galore \
    ${Venus_dir}/chimFinder/chimFinder.py \
    ${Venus_dir}/hisat2/indices/human/human \
    ${Venus_dir}/hisat2/indices/HXB2/HXB2
done



