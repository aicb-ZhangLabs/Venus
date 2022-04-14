#!/bin/bash
#SBATCH --job-name=testing_singleMEGAV
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-12/testing_singleMEGAV.log

STAR_dir=/srv/disk00/cheyul1/Venus/STAR
Venus_dir=/srv/disk00/cheyul1/Venus/repo

indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/Sngl_HIV-only
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-12/testing_singleMEGAV

# Testing Detection Module
#python3 ${Venus_dir}/module-detection.py \
#    --read ${data_dir}/SRR6944349.1_1.fastq.gz ${data_dir}/SRR6944349.1_2.fastq.gz \
#    --virusThreshold 5 \
#    --virusGenome ${indices_dir}/virus.genomeDir \
#    --humanGenome ${indices_dir}/human.genomeDir \
#    --out ${out_dir} \
#    --readFilesCommand zcat \
#    --thread 32
    
python3 ${Venus_dir}/module-detection.py \
    --read ${data_dir}/SRR12165309.1_3.fastq.gz ${data_dir}/SRR12165309.1_2.fastq.gz \
    --virusThreshold 5 \
    --virusGenome ${indices_dir}/virus.genomeDir \
    --humanGenome ${indices_dir}/human.genomeDir \
    --out ${out_dir} \
    --readFilesCommand zcat \
    --thread 32 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist ${indices_dir}/3M-february-2018.txt


## Testing Indexing Module
#python3 ${Venus_dir}/module-index.py \
#    --humanGenome ${out_dir}/human.genomeDir \
#    --humanFASTA ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.fna \
#    --humanGTF ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.gtf \
#    --virusGenome ${out_dir}/virus.genomeDir \
#    --virusFASTA ${indices_dir}/new_virus.genomeDir/mega-virus.fa \
#    --out ${out_dir} \
#    --thread 32

## Printing Test
#python3 module-index.py \
#    --out out_dir \
#    --humanFASTA human.fa \
#    --humanGTF human.gff \
#    --virusFASTA virus.fa \
#    --virusGTF virus.gff \
#    --humanGenome humanGenome_dir \
#    --virusGenome virusGenome_dir
#python3 module-detection.py \
#    --read read1.fq read2.fq \
#    --virusGenome virusGenome \
#    --humanGenome humanGenome \
#    --out ../data
#python3 module-detection.py \
#    --read read1.fq read2.fq \
#    --virusThreshold 5 \
#    --virusGenome virusGenome \
#    --humanGenome humanGenome \
#    --out ../data \
#    --singleCellBarcode 1 16 \
#    --singleUniqueMolIdent 17 10 \
#    --singleWhitelist whitelist.txt
