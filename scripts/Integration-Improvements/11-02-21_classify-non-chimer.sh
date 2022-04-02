#!/bin/bash
#SBATCH --job-name=sanityCK-rmFASTQ
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/11-02-21/sanityCK-rmFASTQ.log

run_prefix=GSE112576
out_dir=/srv/disk00/cheyul1/outputs/11-02-21/sanityCK-rmFASTQ
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

run_dir=${out_dir}/${run_prefix}

STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/unspliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-unspliced.fastq \
--chimSegmentMin 10 \
--chimOutType SeparateSAMold \
--alignIntronMin 20 --alignIntronMax 500000 --winBinNbits 7 \
--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 \
--outSAMmultNmax 1

STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/spliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-spliced.fastq \
--chimSegmentMin 10 \
--chimOutType SeparateSAMold \
--alignIntronMin 20 --alignIntronMax 500000 --winBinNbits 7 \
--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 \
--outSAMmultNmax 1
