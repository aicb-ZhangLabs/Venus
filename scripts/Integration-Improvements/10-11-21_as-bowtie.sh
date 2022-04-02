#!/bin/bash
#SBATCH --job-name=JunctionFilter
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/10-21-21/JunctionFilter.log

run_prefix=GSE126230
out_dir=/srv/disk00/cheyul1/outputs/10-21-21/JunctionFilter
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

fastq_links=(  )
metatable=${STAR_dir}/new_virus.species.txt

for i in {0..70}; do
    ######================================ Create directories ================================######
    cd ${out_dir}/
    mkdir "${run_prefix}_run$i"
    cd "${run_prefix}_run$i"
    mkdir "human_R1"
    mkdir "human_R2"
    mkdir "splice_human_R1"
    mkdir "splice_human_R2"
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
    splice_human_R1_dir=${run_dir}/splice_human_R1
    splice_human_R2_dir=${run_dir}/splice_human_R2
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
    --outSAMmultNmax 1 \
    --alignIntronMax 1 --winBinNbits 7 \
    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27
    
    STAR \
    --runThreadN 32 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${virus_R2_dir}/output/ \
    --genomeDir ${indices_dir}/HIV.genomeDir/ \
    --readFilesIn ${fastqs}/${filenames[1]} \
    --outSAMtype BAM Unsorted \
    --outSAMmultNmax 1 \
    --alignIntronMax 1 --winBinNbits 7 \
    --scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 \
    --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 \
    --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27


    #######================================ Converts aligned BAM file into fastq files ================================######
    fastqs="/srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run$i"
    
    cd ${virus_R1_dir}/output/
    samtools fastq -@ 32 ${virus_R1_dir}/output/Aligned.out.bam > ${human_R1_dir}/fastqs/aligned_R1.fastq
    samtools view -h -o ${virus_R1_dir}/output/Aligned.out.sam ${virus_R1_dir}/output/Aligned.out.bam
    grep "^SRR" ${virus_R1_dir}/output/Aligned.out.sam | cut -f 1 > ${virus_R1_dir}/output/R1.reads
    
    cd ${virus_R2_dir}/output/
    samtools fastq -@ 32 ${virus_R2_dir}/output/Aligned.out.bam > ${human_R2_dir}/fastqs/aligned_R2.fastq
    samtools view -h -o ${virus_R2_dir}/output/Aligned.out.sam ${virus_R2_dir}/output/Aligned.out.bam
    grep "^SRR" ${virus_R2_dir}/output/Aligned.out.sam | cut -f 1 > ${virus_R2_dir}/output/R2.reads

#    for reads in R2's Aligned.out.sam not in R1's Aligned.out.sam
#        Add to R1's fastq
    cd ${human_R1_dir}/fastqs/
    while read name
    do
        if ! grep -Fxq "$name" ${virus_R1_dir}/output/R1.reads ; then
            echo $name
            zgrep -A 3 -m 1 -w "$name" ${fastqs}/${filenames[0]} >> ${human_R1_dir}/fastqs/aligned_R1.fastq
        fi
    done < ${virus_R2_dir}/output/R2.reads
    
    cd ${human_R2_dir}/fastqs/
    while read name
    do
        if ! grep -Fxq "$name" ${virus_R2_dir}/output/R2.reads ; then
            echo $name
            zgrep -A 3 -m 1 -w "$name" ${fastqs}/${filenames[1]} >> ${human_R2_dir}/fastqs/aligned_R2.fastq
        fi
    done < ${virus_R1_dir}/output/R1.reads
    
    
    #######================================ Maps to hybrid index ================================######
    
    
    unset filenames
done



