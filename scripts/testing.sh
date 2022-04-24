#!/bin/bash
#SBATCH --job-name=testing_indexHybrid
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-22/testing_indexHybrid.log

STAR_dir=/srv/disk00/cheyul1/Venus/STAR
Venus_dir=/srv/disk00/cheyul1/Venus/repo

indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/test_files
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-22/testing_indexHybrid

# Testing Detection Module
#python3 ${Venus_dir}/module-detection.py \
#    --read ${data_dir}/SRR6944349.1_1.fastq.gz ${data_dir}/SRR6944349.1_2.fastq.gz \
#    --virusThreshold 5 \
#    --virusChrRef repo/reference_files/virus_chr-ref.tsv \
#    --virusGenome ${indices_dir}/virus.genomeDir \
#    --humanGenome ${indices_dir}/human.genomeDir \
#    --out ${out_dir} \
#    --readFilesCommand zcat \
#    --thread 32
    
#python3 ${Venus_dir}/module-detection.py \
#    --read ${data_dir}/SRR12165309.1_3.fastq.gz ${data_dir}/SRR12165309.1_2.fastq.gz \
#    --virusThreshold 5 \
#    --virusChrRef repo/reference_files/virus_chr-ref.tsv \
#    --virusGenome ${indices_dir}/virus.genomeDir \
#    --humanGenome ${indices_dir}/human.genomeDir \
#    --out ${out_dir} \
#    --readFilesCommand zcat \
#    --thread 32 \
#    --singleCellBarcode 1 16 \
#    --singleUniqueMolIdent 17 10 \
#    --singleWhitelist ${indices_dir}/3M-february-2018.txt

# Testing Integration Module
#python3 ${Venus_dir}/module-integration.py \
#    --read ${data_dir}/singlecell_1cDNA.fastq.gz \
#    --virusGenome ${indices_dir}/HIV.genomeDir \
#    --hybridGenome ${indices_dir}/hg38HIV.genomeDir \
#    --guideFASTA ${indices_dir}/integrSeq.genomeDir/integrSeq.fna \
#    --readFilesCommand zcat \
#    --out ${out_dir} \
#    --virusChr NC_001802.1 \
#    --thread 32 \
#    --geneBed /srv/disk00/cheyul1/genes.bed


## Testing Indexing Module
python3 ${Venus_dir}/module-index.py \
    --hGenome ${out_dir}/hybrid.genomeDir \
    --humanFASTA ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF ${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusGenome ${out_dir}/HIV.genomeDir \
    --virusFASTA ${indices_dir}/HIV.genomeDir/NC_001802.fna \
    --virusGTF ${indices_dir}/HIV.genomeDir/NC_001802.gtf \
    --module integration \
    --out ${out_dir} \
    --thread 32

### Printing Test
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
#    --virusChrRef virus_chr-ref.tsv \
#    --humanGenome humanGenome \
#    --out ../data
#python3 module-detection.py \
#    --read read1.fq read2.fq \
#    --virusThreshold 5 \
#    --virusGenome virusGenome \
#    --virusChrRef virus_chr-ref.tsv \
#    --humanGenome humanGenome \
#    --out ../data \
#    --singleCellBarcode 1 16 \
#    --singleUniqueMolIdent 17 10 \
#    --singleWhitelist whitelist.txt
#python3 module-integration.py \
#    --read read1.fq \
#    --virusGenome virusGenome \
#    --hybridGenome hybridGenome \
#    --guideFASTA guide.fa \
#    --out ../data \
#    --virusChr NC_001802.1
