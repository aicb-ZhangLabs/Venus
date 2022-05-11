#!/bin/bash
#SBATCH --job-name=testing1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-05-10/testing1.log

repo_dir=/srv/disk00/cheyul1/Venus/outputs/22-05-10/Venus
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-05-10/testing1


## Testing Indexing Module
python3 ${repo_dir}/src/module-index/module-index.py \
    --humanFASTA ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA ${repo_dir}/reference_files/mega-virus.fasta \
    --module detection \
    --thread 32 \
    --out ${out_dir}/indices \
    --hGenome ${out_dir}/indices/human.genomeDir \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir

python3 ${repo_dir}/src/module-index/module-index.py \
    --humanFASTA ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA ${repo_dir}/reference_files/NC_001802.fna \
    --virusGTF ${repo_dir}/reference_files/NC_001802.gtf \
    --module detection \
    --thread 32 \
    --out ${out_dir}/indices \
    --hGenome ${out_dir}/indices/human2.genomeDir \
    --virusGenome ${out_dir}/indices/HIV2.genomeDir
    
python3 ${repo_dir}/src/module-index/module-index.py \
    --humanFASTA ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF ${repo_dir}/reference_files/GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA ${repo_dir}/reference_files/NC_001802.fna \
    --virusGTF ${repo_dir}/reference_files/NC_001802.gtf \
    --module integration \
    --thread 32 \
    --out ${out_dir}/indices \
    --hGenome ${out_dir}/indices/hybrid.genomeDir \
    --virusGenome ${out_dir}/indices/HIV.genomeDir



# Testing Detection Module
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --out ${out_dir}/detection/bulk_single-end
    
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz ${repo_dir}/test_data/bulk_2.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --out ${out_dir}/detection/bulk_paired-end
    
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/singlecell_1cDNA.fastq.gz ${repo_dir}/test_data/singlecell_2CB+UMI.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist ${repo_dir}/reference_files/3M-february-2018.txt \
    --out ${out_dir}/detection/single-cell



# Testing Integration Module
python3 ${repo_dir}/src/module-integration/module-integration.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz \
    --virusGenome ${out_dir}/indices/HIV.genomeDir \
    --hybridGenome ${out_dir}/indices/hybrid.genomeDir \
    --guideFASTA ${repo_dir}/reference_files/integrSeq.fna \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed ${repo_dir}/reference_files/genes.bed \
    --out ${out_dir}/integration/single-end

python3 ${repo_dir}/src/module-integration/module-integration.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz ${repo_dir}/test_data/bulk_2.fastq.gz \
    --virusGenome ${out_dir}/indices/HIV.genomeDir \
    --hybridGenome ${out_dir}/indices/hybrid.genomeDir \
    --guideFASTA ${repo_dir}/reference_files/integrSeq.fna \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed ${repo_dir}/reference_files/genes.bed \
    --out ${out_dir}/integration/paired-end
