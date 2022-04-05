#!/bin/bash
#SBATCH --job-name=generateTest
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/Venus/logs/22-04-05/generateTest.log

run_prefix=generateTest
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-04-05
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices
data_dir=/srv/disk00/cheyul1/Venus/datasets/YaChi_GSE112576
numfiles=1

fastq_links=( https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR6944369/SRR6944369.1 )
metatable=${STAR_dir}/new_virus.species.txt

for i in ${!fastq_links[@]}; do
    ######================================ Downloads SRA runs and create directories ================================######

    run_dir="${out_dir}/${run_prefix}_run$i"
    human_dir=${run_dir}/human
    virus_dir=${run_dir}/virus
    mega_virus_dir=${run_dir}/mega_virus
#    hybrid_dir=${run_dir}/hybrid
#    hybrid_sam_dir=${run_dir}/hybrid_sam


    ############################################################ VIRUS DETECTION ############################################################

    ######================================ Maps to Human index ================================######
    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/

    if [[ $numfiles == 1 ]]; then
        STAR \
        --runThreadN 32 \
        --readFilesCommand zcat \
        --outFileNamePrefix ${human_dir}/output/ \
        --genomeDir ${indices_dir}/human.genomeDir/ \
        --readFilesIn ${data_dir}/SRR6944349.1_1.fastq.gz \
        --outReadsUnmapped Fastx \
        --outSAMtype None
    else
        STAR \
        --runThreadN 32 \
        --readFilesCommand zcat \
        --outFileNamePrefix ${human_dir}/output/ \
        --genomeDir ${indices_dir}/human.genomeDir/ \
        --readFilesIn ${fastqs}/${filenames[0]} ${fastqs}/${filenames[1]} \
        --outReadsUnmapped Fastx \
        --outSAMtype None
    fi


    #######================================ Maps to a Target Virus index ================================######
    fastqs=${human_dir}/output
#    cd ${fastqs}/
    mv ${fastqs}/Unmapped.out.mate1 ${fastqs}/Unmapped.out.mate1.fastq

    if test -f "Unmapped.out.mate2"; then
        mv Unmapped.out.mate2 Unmapped.out.mate2.fastq

        STAR \
        --runThreadN 32 \
        --outFileNamePrefix ${virus_dir}/output/ \
        --genomeDir ${indices_dir}/HIV.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate
    else
        STAR \
        --runThreadN 32 \
        --outFileNamePrefix ${virus_dir}/output/ \
        --genomeDir ${indices_dir}/HIV.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype SAM
    fi


    #######================================ Maps to Mega-Virus index ================================#######
    fastqs=${human_dir}/output
#    cd ${fastqs}/

    if test -f "Unmapped.out.mate2.fastq"; then
        STAR \
        --runThreadN 32 \
        --outFileNamePrefix ${mega_virus_dir}/output/ \
        --genomeDir ${indices_dir}/new_virus.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype SAM
    else
        STAR \
        --runThreadN 32 \
        --outFileNamePrefix ${mega_virus_dir}/output/ \
        --genomeDir ${indices_dir}/new_virus.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype SAM
    fi

    for file in *.fastq; do
        rm -f $file
    done
    
    unset filenames
done


######================================ Creates mega-index virus species output ================================######
cd ${mega_virus_dir}/output/
#samtools view -h -o Aligned.sortedByCoord.sam Aligned.sortedByCoord.out.bam

grep -v @ Aligned.out.sam |cut -f 3|sort|uniq > uniq_aligned.txt
total=$(grep -v @ Aligned.sortedByCoord.sam|wc -l)
while read line
do
    counts=$(grep -v @ Aligned.sortedByCoord.sam | grep ${line} |wc -l)
    species=$(grep ${line} ${metatable}|cut -f 2)
    echo -e "$line\t$species\t$counts\t$total\t"$(bc -l <<< "(($counts / $total)*100)")"%" >> species.counts.txt
done<uniq_aligned.txt

sort -t $'\t' -k 3,3 -r -V species.counts.txt > ordered_species.counts.txt



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


