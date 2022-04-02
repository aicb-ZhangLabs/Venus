#!/bin/bash
#SBATCH --job-name=2GSE126230
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/10-06-21/GSE126230.log

run_prefix=GSE126230
out_dir=/srv/disk00/cheyul1/outputs/10-06-21
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

fastq_links=(  )
metatable=${STAR_dir}/new_virus.species.txt

for i in {40..70}; do
    ######================================ Create directories ================================######
    cd ${out_dir}/
    mkdir "${run_prefix}_run$i"
    cd "${run_prefix}_run$i"
    mkdir "human_R1"
    mkdir "human_R2"
    mkdir "virus_R1"
    mkdir "virus_R2"
#    mkdir "hybrid"
#    mkdir "hybrid_sam"
    cd human_R1
    mkdir "fastqs"
    mkdir "output"
    cd ../human_R2
    mkdir "fastqs"
    mkdir "output"
    cd ../virus_R1
    mkdir "output"
    cd ../virus_R2
    mkdir "output"
#    cd ../hybrid
#    mkdir "output"
#    cd ../hybrid_sam
#    mkdir "output"

    run_dir="${out_dir}/${run_prefix}_run$i"
    human_R1_dir=${run_dir}/human_R1
    human_R2_dir=${run_dir}/human_R2
    virus_R1_dir=${run_dir}/virus_R1
    virus_R2_dir=${run_dir}/virus_R2
#    hybrid_dir=${run_dir}/hybrid
#    hybrid_sam_dir=${run_dir}/hybrid_sam

#    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/
#    wget -nv ${fastq_links[$i]}
#    parallel-fastq-dump -s ${human_dir}/fastqs/* --threads 32 --split-files --tmpdir ${human_dir}/fastqs/tmp/ --outdir ${human_dir}/fastqs/
#    for file in *; do
#        fastq-dump $file --split-files --outdir .;
#    done

#    for file in *; do
#        if [[ $file != *.fastq ]]; then
#            rm -rf $file
#        fi
#    done

    fastqs="/srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run$i"
    cd ${fastqs}/
    numfiles=$(ls | wc -l)
    declare -a filenames
    for file in *.fastq.gz; do
        filenames=(${filenames[@]} "$file")
    done


    ############################################################ VIRUS DETECTION ############################################################


    ######================================ Maps to HIV index ================================######
    fastqs="/srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run$i"
    cd ${fastqs}/

    STAR \
    --runThreadN 32 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${virus_R1_dir}/output/ \
    --genomeDir ${indices_dir}/HIV.genomeDir/ \
    --readFilesIn ${fastqs}/${filenames[0]} \
    --outSAMtype BAM Unsorted \
    --seedSearchStartLmax 20
    
    STAR \
    --runThreadN 32 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${virus_R2_dir}/output/ \
    --genomeDir ${indices_dir}/HIV.genomeDir/ \
    --readFilesIn ${fastqs}/${filenames[1]} \
    --outSAMtype BAM Unsorted \
    --seedSearchStartLmax 20


    #######================================ Converts aligned BAM file into fastq files ================================######
    cd ${virus_R1_dir}/output/
    samtools fastq \
    -@ 32 \
    ${virus_R1_dir}/output/Aligned.out.bam > ${human_R1_dir}/fastqs/aligned_R1.fastq
    
    cd ${virus_R2_dir}/output/
    samtools fastq \
    -@ 32 \
    ${virus_R2_dir}/output/Aligned.out.bam > ${human_R2_dir}/fastqs/aligned_R2.fastq
    
    
    #######================================ Map aligned fastq to Human index ================================######
    fastqs=${human_R1_dir}/fastqs
    cd ${fastqs}/
    STAR \
    --runThreadN 32 \
    --outFileNamePrefix ${human_R1_dir}/output/ \
    --genomeDir ${indices_dir}/human.genomeDir/ \
    --readFilesIn ${fastqs}/aligned_R1.fastq \
    --seedSearchStartLmax 20
    
    
    fastqs=${human_R2_dir}/fastqs
    cd ${fastqs}/
    STAR \
    --runThreadN 32 \
    --outFileNamePrefix ${human_R2_dir}/output/ \
    --genomeDir ${indices_dir}/human.genomeDir/ \
    --readFilesIn ${fastqs}/aligned_R2.fastq \
    --seedSearchStartLmax 20



#    ############################################################ VIRUS INTEGRATION ############################################################
#
#    ######================================ Maps to hybrid index for Chimeric Junctions ================================######
#    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/
#    if [[ $numfiles == 1 ]]; then
#        STAR \
#        --runThreadN 32 \
#        --outFileNamePrefix ${hybrid_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} \
#        --outReadsUnmapped Fastx \
#        --outSAMtype None \
#        --chimSegmentMin 10 \
#        --chimOutType Junctions
#    else
#        STAR \
#        --runThreadN 32 \
#        --outFileNamePrefix ${hybrid_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} ${fastqs}/${filenames[1]} \
#        --outReadsUnmapped Fastx \
#        --outSAMtype None \
#        --chimSegmentMin 10 \
#        --chimOutType Junctions
#    fi
#
#
#    ######================================ Gets integration sites in Chimeric Junctions files ================================######
#    cd ${hybrid_dir}/output/
#    grep "NC_001802.1" Chimeric.out.junction > Chimeric.filtered.junction
#    awk '$1 != $4' Chimeric.filtered.junction > Chimeric.discordant.tsv
#
#
#    ######================================ Maps to hybrid index for Chimeric SAM file ================================######
#    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/
#    if [[ $numfiles == 1 ]]; then
#        STAR \
#        --runThreadN 32 \
#        --outFileNamePrefix ${hybrid_sam_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} \
#        --outReadsUnmapped Fastx \
#        --outSAMtype None \
#        --chimSegmentMin 10 \
#        --chimOutType SeparateSAMold
#    else
#        STAR \
#        --runThreadN 32 \
#        --outFileNamePrefix ${hybrid_sam_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38HIV.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} ${fastqs}/${filenames[1]} \
#        --outReadsUnmapped Fastx \
#        --outSAMtype None \
#        --chimSegmentMin 10 \
#        --chimOutType SeparateSAMold
#    fi
#
#    for file in *.fastq; do
#        rm -f $file
#    done
#
#
#    ######================================ Gets Chimeric SAM files for visual in IGV ================================######
#    cd ${hybrid_sam_dir}/output/
#    grep "NC_001802.1" Chimeric.out.sam > Chimeric.filtered.sam     # Can only do this for bulk paired-end seq, because need to have mate read
#    filesize=$(stat -c%s "Chimeric.filtered.sam")
#    if [[ $filesize > 0 ]]; then
#        samtools view -bt /srv/disk00/cheyul1/Venus/STAR/indices/hg38HIV.genomeDir/hg38HIV.fna.fai -u Chimeric.filtered.sam | \
#        samtools sort -o Chimeric.filtered.sortedByCoord.bam
#        samtools index -b -@ 32 Chimeric.filtered.sortedByCoord.bam Chimeric.filtered.sortedByCoord.bam.bai
#    else
#        echo "No Chimeric Transcripts Found"
#    fi
    
    unset filenames
done


