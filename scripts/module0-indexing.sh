#!/bin/bash
#SBATCH --job-name=integrSeq-index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/02-02-22/integrSeq-index.log


out_dir=/srv/disk00/cheyul1/outputs/02-02-22/integrSeq-index
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices


#### Creating a fna index ###
#samtools faidx /srv/disk00/cheyul1/Venus/STAR/indices/hg38HIV.genomeDir/hg38HIV.fna

STAR \
--runThreadN 16 \
--outFileNamePrefix ${out_dir}/ \
--runMode genomeGenerate \
--genomeDir ${indices_dir}/integrSeq.genomeDir/ \
--genomeFastaFiles ${indices_dir}/integrSeq.genomeDir/integrSeq.fna \
--genomeSAindexNbases 2
#--sjdbGTFfile ${indices_dir}/CoV.genomeDir/CoV.gtf \
##--sjdbGTFfeatureExon mature_protein_region_of_CDS

#bowtie2-build \
#/srv/disk00/cheyul1/Venus/bowtie2/indices/HXB2/HXB2.fa \
#/srv/disk00/cheyul1/Venus/bowtie2/indices/HXB2/HXB2 \
#--threads 32

#hisat2-build \
${indices_dir}/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.fna \
/srv/disk00/cheyul1/Venus/hisat2/indices/human/human \
-p 32

hisat2-build -p 32 \
--exon /srv/disk00/cheyul1/Venus/hisat2/human.exon \
--ss /srv/disk00/cheyul1/Venus/hisat2/human.ss \
/srv/disk00/cheyul1/Venus/STAR/indices/human.genomeDir/GCF_000001405.39_GRCh38.p13_genomic.fna \
/srv/disk00/cheyul1/Venus/hisat2/indices/human/human


