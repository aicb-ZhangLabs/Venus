#!/bin/bash
#SBATCH --job-name=fsn-bam
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/11-10-21/fsn-bam.log

run_prefix=GSE112576
out_dir=/srv/disk00/cheyul1/outputs/11-10-21/fsn-bam
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices


for i in {0..39}; do
    run_dir=${out_dir}/${run_prefix}-run${i}
    
    ###############################################################################################
    ######================================ Trims Low-Q Reads ================================######
    ###############################################################################################
    mkdir --parents ${run_dir}/R1/trim_fastq/
    trim_galore \
    --trim-n \
    --quality 5 \
    --phred33 \
    --length 20 \
    --output_dir ${run_dir}/R1/trim_fastq/ \
    --gzip \
    --cores 7 \
    /srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run${i}/*_R1.fastq.gz

    mkdir --parents ${run_dir}/R2/trim_fastq/
    trim_galore \
    --trim-n \
    --quality 5 \
    --phred33 \
    --length 20 \
    --output_dir ${run_dir}/R2/trim_fastq/ \
    --gzip \
    --cores 7 \
    /srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run${i}/*_R2.fastq.gz
    
    
    
    
    ###############################################################################################
    ######================================ Maps to HIV index ================================######
    ###############################################################################################
    
#    STAR \
#    --runThreadN 32 \
#    --readFilesCommand zcat \
#    --outFileNamePrefix ${run_dir}/R1/aln_fastq/ \
#    --genomeDir ${indices_dir}/HIV.genomeDir/ \
#    --readFilesIn ${run_dir}/R1/trim_fastq/*_trimmed.fq.gz \
#    --outSAMtype BAM Unsorted \
#    --outSAMmultNmax 1 \
#    --alignIntronMax 1 --winBinNbits 7 \
#    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
#    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
#    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27
#
#    STAR \
#    --runThreadN 32 \
#    --readFilesCommand zcat \
#    --outFileNamePrefix ${run_dir}/R2/aln_fastq/ \
#    --genomeDir ${indices_dir}/HIV.genomeDir/ \
#    --readFilesIn ${run_dir}/R2/trim_fastq/*_trimmed.fq.gz \
#    --outSAMtype BAM Unsorted \
#    --outSAMmultNmax 1 \
#    --alignIntronMax 1 --winBinNbits 7 \
#    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
#    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
#    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27
    
    
    
    
    #########################################################################################################################
    #######================================ Converts aligned BAM file into fastq files ================================######
    #########################################################################################################################
    
#    samtools fastq -@ 31 ${run_dir}/R1/aln_fastq/Aligned.out.bam > ${run_dir}/R1/aln_fastq/aligned_R1.fastq
#    samtools view -h -o ${run_dir}/R1/aln_fastq/Aligned.out.sam ${run_dir}/R1/aln_fastq/Aligned.out.bam
#
#    samtools fastq -@ 31 ${run_dir}/R2/aln_fastq/Aligned.out.bam > ${run_dir}/R2/aln_fastq/aligned_R2.fastq
#    samtools view -h -o ${run_dir}/R2/aln_fastq/Aligned.out.sam ${run_dir}/R2/aln_fastq/Aligned.out.bam
    
    
    
    
    ##################################################################################################
    ######================================ Maps to Hybrid index ================================######
    ##################################################################################################
    
    STAR \
    --runThreadN 32 \
    --outFileNamePrefix ${run_dir}/R1/ \
    --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
    --readFilesIn /srv/disk00/cheyul1/outputs/10-27-21/chimr2-star/${run_prefix}-run${i}/R1/aln_fastq/aligned_R1.fastq \
    --alignIntronMin 20 --winBinNbits 7 \
    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 \
    --outSAMtype BAM Unsorted --outSAMmultNmax 1 \
    --chimOutType WithinBAM SoftClip --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
    --chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
    --alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
    --alignInsertionFlush Right
    
#    # Chimeric.out.junction
#    echo -en '\n' >> ${out_dir}/integration.junction
#    echo "Run ${i} read ONE (1)" >> ${out_dir}/integration.junction
#    grep "NC_001802.1" ${run_dir}/R1/Chimeric.out.junction | cut -f 1,2 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
#    grep "NC_001802.1" ${run_dir}/R1/Chimeric.out.junction | cut -f 4,5 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
    
    # Aligned.out.bam
    samtools view -h -o ${run_dir}/R1/Aligned.out.sam ${run_dir}/R1/Aligned.out.bam
    grep -v "@" ${run_dir}/R1/Aligned.out.sam >> ${out_dir}/integration.sam


    STAR \
    --runThreadN 32 \
    --outFileNamePrefix ${run_dir}/R2/ \
    --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
    --readFilesIn /srv/disk00/cheyul1/outputs/10-27-21/chimr2-star/${run_prefix}-run${i}/R2/aln_fastq/aligned_R2.fastq \
    --alignIntronMin 20 --winBinNbits 7 \
    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 \
    --outSAMtype BAM Unsorted --outSAMmultNmax 1 \
    --chimOutType WithinBAM SoftClip --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
    --chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 --chimNonchimScoreDropMin 10 \
    --alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
    --alignInsertionFlush Right
    
#    # Chimeric.out.junction
#    echo -en '\n' >> ${out_dir}/integration.junction
#    echo "Run ${i} read TWO (2)" >> ${out_dir}/integration.junction
#    grep "NC_001802.1" ${run_dir}/R2/Chimeric.out.junction | cut -f 1,2 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
#    grep "NC_001802.1" ${run_dir}/R2/Chimeric.out.junction | cut -f 4,5 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
    
    # Aligned.out.bam
    samtools view -h -o ${run_dir}/R2/Aligned.out.sam ${run_dir}/R2/Aligned.out.bam
    grep -v "@" ${run_dir}/R2/Aligned.out.sam >> ${out_dir}/integration.sam
done




########################################################################################################
######================================ Collects integration genes ================================######
########################################################################################################
#grep "NC_" ${out_dir}/integration.junction | sort -t $'\t' -k 1,1 -k 2,2 -V | uniq | \
# sed -e "s/NC_000001.11/chr1/g" | sed -e "s/NC_000002.12/chr2/g" | sed -e "s/NC_000003.12/chr3/g" | \
# sed -e "s/NC_000004.12/chr4/g" | sed -e "s/NC_000005.10/chr5/g" | sed -e "s/NC_000006.12/chr6/g" | \
# sed -e "s/NC_000007.14/chr7/g" | sed -e "s/NC_000008.11/chr8/g" | sed -e "s/NC_000009.12/chr9/g" | \
# sed -e "s/NC_000010.11/chr10/g" | sed -e "s/NC_000011.10/chr11/g" | sed -e "s/NC_000012.12/chr12/g" | \
# sed -e "s/NC_000013.11/chr13/g" | sed -e "s/NC_000014.9/chr14/g" | sed -e "s/NC_000015.10/chr15/g" | \
# sed -e "s/NC_000016.10/chr16/g" | sed -e "s/NC_000017.11/chr17/g" | sed -e "s/NC_000018.10/chr18/g" | \
# sed -e "s/NC_000019.10/chr19/g" | sed -e "s/NC_000020.11/chr20/g" | sed -e "s/NC_000021.9/chr21/g" | \
# sed -e "s/NC_000022.11/chr22/g" | sed -e "s/NC_000023.11/chrX/g" | sed -e "s/NC_000024.10/chrY/g" \
# > ${out_dir}/integration.sites
#
#while read site
#do
#    cand=( $site )
#    echo $site
##    echo ${cand[1]}
#    while read gene
#    do
#        position=( $gene )
#        if [[ ${cand[0]} == ${position[0]} ]]; then
#            if [[ ${cand[1]} -gt ${position[1]} ]]; then
#                if [[ ${cand[1]} -lt ${position[2]} ]]; then
#                    echo ${position[0]}
#                    echo ${position[1]}
#                    echo ${position[2]}
#                    echo ${position[3]}
#                    echo ${position[3]} >> ${out_dir}/integration-all.genes
#                fi
#            fi
#        fi
#    done < /srv/disk00/cheyul1/genes.bed
#done < ${out_dir}/integration.sites
#sort ${out_dir}/integration-all.genes | uniq > ${out_dir}/integration.genes




###########################################################################################################
######================================ Prepares visualization in IGV ================================######
###########################################################################################################
samtools view -@ 31 -bt /srv/disk00/cheyul1/Venus/STAR/indices/hg38HIV.genomeDir/hg38HIV.fna.fai -u ${out_dir}/integration.sam | \
samtools sort -@ 32 -o ${out_dir}/integration.Sorted.bam
samtools index -@ 31 -b ${out_dir}/integration.Sorted.bam ${out_dir}/integration.Sorted.bam.bai
