#!/bin/bash
#SBATCH --job-name=miss-fusion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/11-09-21/miss-fusion.log

run_prefix=GSE112576
out_dir=/srv/disk00/cheyul1/outputs/11-09-21/miss-fusion
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

run_dir=${out_dir}/${run_prefix}

### Chimeric mapping
STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/unspliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-unspliced.fastq \
--outSAMtype BAM Unsorted --outSAMmultNmax 1 \
--chimOutType Junctions --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
--alignInsertionFlush Right


### Chimeric mapping
STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/spliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-spliced.fastq \
--outSAMtype BAM Unsorted --outSAMmultNmax 1 \
--chimOutType Junctions --chimSegmentMin 8 --chimJunctionOverhangMin 8 \
--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
--alignInsertionFlush Right


### Chimeric mapping + Bowtie settings (except with splicing)
STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/unspliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-unspliced.fastq \
--outSAMtype BAM Unsorted --outSAMmultNmax 1 \
--chimOutType Junctions --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
--alignInsertionFlush Right \

--alignIntronMin 20 --winBinNbits 7 \
--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27


### Chimeric mapping + Bowtie settings (except with splicing)
STAR \
--runThreadN 32 \
--outFileNamePrefix ${run_dir}/spliced/ \
--genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
--readFilesIn /srv/disk00/cheyul1/missed-spliced.fastq \
--outSAMtype BAM Unsorted --outSAMmultNmax 1 \
--chimOutType Junctions --chimSegmentMin 8 --chimJunctionOverhangMin 8 \
--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
--alignInsertionFlush Right \

--alignIntronMin 20 --winBinNbits 7 \
--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27


########################################################################################################
######================================ Collects integration genes ================================######
########################################################################################################
grep "^SRR" ${out_dir}/integration.sam | sort -t $'\t' -k 3,3 -k 4,4 -V | cut -f 3,4 | uniq | \
 sed -e "s/NC_000001.11/chr1/g" | sed -e "s/NC_000002.12/chr2/g" | sed -e "s/NC_000003.12/chr3/g" | \
 sed -e "s/NC_000004.12/chr4/g" | sed -e "s/NC_000005.10/chr5/g" | sed -e "s/NC_000006.12/chr6/g" | \
 sed -e "s/NC_000007.14/chr7/g" | sed -e "s/NC_000008.11/chr8/g" | sed -e "s/NC_000009.12/chr9/g" | \
 sed -e "s/NC_000010.11/chr10/g" | sed -e "s/NC_000011.10/chr11/g" | sed -e "s/NC_000012.12/chr12/g" | \
 sed -e "s/NC_000013.11/chr13/g" | sed -e "s/NC_000014.9/chr14/g" | sed -e "s/NC_000015.10/chr15/g" | \
 sed -e "s/NC_000016.10/chr16/g" | sed -e "s/NC_000017.11/chr17/g" | sed -e "s/NC_000018.10/chr18/g" | \
 sed -e "s/NC_000019.10/chr19/g" | sed -e "s/NC_000020.11/chr20/g" | sed -e "s/NC_000021.9/chr21/g" | \
 sed -e "s/NC_000022.11/chr22/g" | sed -e "s/NC_000023.11/chrX/g" | sed -e "s/NC_000024.10/chrY/g" \
 > ${out_dir}/integration.sites
 
while read site
do
    cand=( $site )
    echo $site
#    echo ${cand[1]}
    while read gene
    do
        position=( $gene )
        if [[ ${cand[0]} == ${position[0]} ]]; then
            if [[ ${cand[1]} -gt ${position[1]} ]]; then
                if [[ ${cand[1]} -lt ${position[2]} ]]; then
                    echo ${position[0]}
                    echo ${position[1]}
                    echo ${position[2]}
                    echo ${position[3]}
                    echo ${position[3]} >> ${out_dir}/integration-all.genes
                fi
            fi
        fi
    done < /srv/disk00/cheyul1/genes.bed
done < ${out_dir}/integration.sites
sort ${out_dir}/integration-all.genes | uniq > ${out_dir}/integration.genes

# Visualize integration.sam in IGV
#grep "^SRR" integration.sam.old > integration.sam
#samtools view -bt /srv/disk00/cheyul1/Venus/STAR/indices/hg38HIV.genomeDir/hg38HIV.fna.fai -u ${out_dir}/integration.sam | \
samtools sort -o ${out_dir}/integration.sortedByCoord.bam
#samtools index -b ${out_dir}/integration.sortedByCoord.bam ${out_dir}/integration.sortedByCoord.bam.bai

grep -v "@" Aligned.out.sam | cut -f 3,4 | sort -t $'\t' -k 1,1 -k 2,2 -V | uniq            # Aligned.out.sam
cut -f 1,2,4,5 Chimeric.out.junction | sort -t $'\t' -k 1,1 -k 3,3 -k 2,2 -k 4,4 -V | uniq  # Chimeric.out.junction
