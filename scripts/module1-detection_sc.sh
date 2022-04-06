#!/bin/bash
#SBATCH --job-name=CoVlungs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/09-11-21/CoVlungs.log

run_prefix=CoVlungs
out_dir=/srv/disk00/cheyul1/outputs/09-11-21
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

fastq_links=( https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11181954/SRR11181954.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR11181955/SRR11181955.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR11181957/SRR11181957.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11181956/SRR11181956.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11181958/SRR11181958.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-21/SRR11181959/SRR11181959.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-13/SRR11537950/SRR11537950.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11537949/SRR11537949.2 \
https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-19/SRR11537951/SRR11537951.2 )
metatable=${STAR_dir}/new_virus.species.txt

for i in ${!fastq_links[@]}; do
    ######================================ Downloads SRA runs and create directories ================================######
    cd ${out_dir}/
    mkdir "${run_prefix}_run$i"
    cd "${run_prefix}_run$i"
    mkdir "human"
    mkdir "virus"
    mkdir "mega_virus"
#    mkdir "hybrid"
#    mkdir "hybrid_sam"
    cd human
    mkdir "fastqs"
    mkdir "output"
    cd ../virus
    mkdir "output"
    cd ../mega_virus
    mkdir "output"
#    cd ../hybrid
#    mkdir "output"
#    cd ../hybrid_sam
#    mkdir "output"

    run_dir="${out_dir}/${run_prefix}_run$i"
    human_dir=${run_dir}/human
    virus_dir=${run_dir}/virus
    mega_virus_dir=${run_dir}/mega_virus
#    hybrid_dir=${run_dir}/hybrid
#    hybrid_sam_dir=${run_dir}/hybrid_sam

    fastqs=${human_dir}/fastqs
    cd ${fastqs}/
    wget -nv ${fastq_links[$i]}
    parallel-fastq-dump -s ${human_dir}/fastqs/* --threads 16 --split-files --tmpdir ${human_dir}/fastqs/tmp/ --outdir ${human_dir}/fastqs/
#    for file in *; do
#        fastq-dump $file --split-files --outdir .
#    done

    for file in *; do
        if [[ $file != *.fastq ]]; then
            rm -rf $file
        fi
    done

    cd ${fastqs}/
    numfiles=$(ls | wc -l)
    declare -a filenames
    for file in *.fastq; do
        filenames=(${filenames[@]} "$file")
    done



    ############################################################ VIRUS DETECTION ############################################################
    
#    ###############################################################################################
#    ######================================ Trims Low-Q Reads ================================######
#    ###############################################################################################
#    mkdir --parents ${run_dir}/R1/trim_fastq/
#    trim_galore \
#    --trim-n \
#    --quality 5 \
#    --phred33 \
#    --length 20 \
#    --output_dir ${run_dir}/R1/trim_fastq/ \
#    --gzip \
#    --cores 7 \
#    /srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run${i}/*_R1.fastq.gz
#
#    mkdir --parents ${run_dir}/R2/trim_fastq/
#    trim_galore \
#    --trim-n \
#    --quality 5 \
#    --phred33 \
#    --length 20 \
#    --output_dir ${run_dir}/R2/trim_fastq/ \
#    --gzip \
#    --cores 7 \
#    /srv/disk00/cheyul1/outputs/09-24-21/fastqs/${run_prefix}-run${i}/*_R2.fastq.gz


    ######================================ Maps to Human index ================================######
    fastqs=${human_dir}/fastqs
    cd ${fastqs}/
    
    if [[ $numfiles == 1 ]]; then
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${human_dir}/output/ \
        --genomeDir ${indices_dir}/human.genomeDir/ \
        --readFilesIn ${fastqs}/${filenames[0]} \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --outSAMattributes NH HI nM AS GX GN CR CB CY UR UB UY sS sQ sM
    else
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${human_dir}/output/ \
        --genomeDir ${indices_dir}/human.genomeDir/ \
        --readFilesIn ${fastqs}/${filenames[1]} ${fastqs}/${filenames[0]} \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --outSAMattributes NH HI nM AS GX GN CR CB CY UR UB UY sS sQ sM
    fi
    
    
    #######================================ Maps to a Target Virus index ================================######
    fastqs=${human_dir}/output
    cd ${fastqs}/
    mv Unmapped.out.mate1 Unmapped.out.mate1.fastq
    
    if test -f "Unmapped.out.mate2"; then
        mv Unmapped.out.mate2 Unmapped.out.mate2.fastq
    
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${virus_dir}/output/ \
        --genomeDir ${indices_dir}/SARSCoV2.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --outSAMattributes NH HI nM AS GX GN CR CB CY UR UB UY sS sQ sM
    else
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${virus_dir}/output/ \
        --genomeDir ${indices_dir}/SARSCoV2.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq \
        --outFilterMultimapNmax 1 \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS GX GN CR CB CY UR UB UY sS sQ sM
    fi


    #######================================ Creates BAM index for visual in IGV ================================#######
    cd ${virus_dir}/output/
    samtools index -b -@ 16 Aligned.sortedByCoord.out.bam Aligned.sortedByCoord.out.bam.bai
    
    
    #######================================ Maps to Mega-Virus index ================================#######
    fastqs=${human_dir}/output
    cd ${fastqs}/

    if test -f "Unmapped.out.mate2.fastq"; then
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${mega_virus_dir}/output/ \
        --genomeDir ${indices_dir}/new_virus.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq \
        --outFilterMultimapNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --soloType CB_samTagOut \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloCBmatchWLtype 1MM \
        --outSAMattributes NH HI nM AS CB UR
    else
        STAR \
        --runThreadN 16 \
        --outFileNamePrefix ${mega_virus_dir}/output/ \
        --genomeDir ${indices_dir}/new_virus.genomeDir/ \
        --readFilesIn ${fastqs}/Unmapped.out.mate1.fastq \
        --outFilterMultimapNmax 1 \
        --soloType CB_samTagOut \
        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloCBmatchWLtype 1MM \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CB UR
    fi
    
    fastqs=${human_dir}/output
    cd ${fastqs}/
    for file in *.fastq; do
        rm -f $file
    done


    ######================================ Creates mega-index virus species output ================================######
    cd ${mega_virus_dir}/output/
    samtools view -h -o Aligned.sortedByCoord.sam Aligned.sortedByCoord.out.bam

    grep -v @ Aligned.sortedByCoord.sam |cut -f 3|sort|uniq > uniq_aligned.txt
    
    filesize=$(stat -c%s "uniq_aligned.txt")
    if [[ $filesize > 0 ]]; then
        total=$(grep -v @ Aligned.sortedByCoord.sam|wc -l)
        while read line
        do
            counts=$(grep -v @ Aligned.sortedByCoord.sam | grep ${line} |wc -l)
            species=$(grep ${line} ${metatable}|cut -f 2)
            echo -e "$line\t$species\t$counts\t$total\t"$(bc -l <<< "(($counts / $total)*100)")"%" >> species.counts.txt
        done<uniq_aligned.txt

        sort -t $'\t' -k 3,3 -r -V species.counts.txt > ordered_species.counts.txt
    else
        echo "No NCBI Viruses Found"
    fi
    unset filesize



#    ############################################################ VIRUS INTEGRATION ############################################################
#
#    ######================================ Maps to hybrid index for Chimeric Junctions ================================######
#    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/
#    if [[ $numfiles == 1 ]]; then
#        STAR \
#        --runThreadN 16 \
#        --outFileNamePrefix ${hybrid_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38SARSCoV2.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} \
#        --outSAMtype None \
#        --soloType CB_UMI_Simple \
#        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
#        --soloCBstart 1 \
#        --soloCBlen 16 \
#        --soloUMIstart 17 \
#        --soloUMIlen 10 \
#        --soloBarcodeReadLength 0 \
#        --chimSegmentMin 10 \
#        --chimOutType Junctions
#    else
#        STAR \
#        --runThreadN 16 \
#        --outFileNamePrefix ${hybrid_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38SARSCoV2.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[1]} ${fastqs}/${filenames[0]} \
#        --outSAMtype None \
#        --soloType CB_UMI_Simple \
#        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
#        --soloCBstart 1 \
#        --soloCBlen 16 \
#        --soloUMIstart 17 \
#        --soloUMIlen 10 \
#        --soloBarcodeReadLength 0 \
#        --chimSegmentMin 10 \
#        --chimOutType Junctions
#    fi
#
#
#    #######================================ Gets integration sites in Chimeric Junctions files ================================######
#    cd ${hybrid_dir}/output/
#    grep "NC_045512.2" Chimeric.out.junction > Chimeric.filtered.junction
#    awk '$1 != $4' Chimeric.filtered.junction > Chimeric.discordant.tsv
#
#
#    #######================================ Maps to hybrid index for Chimeric SAM file ================================#######
#    fastqs=${human_dir}/fastqs
#    cd ${fastqs}/
#    if [[ $numfiles == 1 ]]; then
#        STAR \
#        --runThreadN 16 \
#        --outFileNamePrefix ${hybrid_sam_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38SARSCoV2.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[0]} \
#        --outSAMtype None \
#        --soloType CB_UMI_Simple \
#        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
#        --soloCBstart 1 \
#        --soloCBlen 16 \
#        --soloUMIstart 17 \
#        --soloUMIlen 10 \
#        --soloBarcodeReadLength 0 \
#        --chimSegmentMin 10 \
#        --chimOutType SeparateSAMold
#    else
#        STAR \
#        --runThreadN 16 \
#        --outFileNamePrefix ${hybrid_sam_dir}/output/ \
#        --genomeDir ${indices_dir}/hg38SARSCoV2.genomeDir/ \
#        --readFilesIn ${fastqs}/${filenames[1]} ${fastqs}/${filenames[0]} \
#        --outSAMtype None \
#        --soloType CB_UMI_Simple \
#        --soloCBwhitelist ${indices_dir}/737K-august-2016.txt \
#        --soloCBstart 1 \
#        --soloCBlen 16 \
#        --soloUMIstart 17 \
#        --soloUMIlen 10 \
#        --soloBarcodeReadLength 0 \
#        --chimSegmentMin 10 \
#        --chimOutType SeparateSAMold
#    fi

    fastqs=${human_dir}/fastqs
    cd ${fastqs}/
    for file in *.fastq; do
        rm -f $file
    done


#    #######================================ Gets Chimeric SAM files for visual in IGV ================================#######
#    cd ${hybrid_sam_dir}/output/
#    grep "NC_045512.2" Chimeric.out.sam > Chimeric.filterered.sam     # Cannot filter only for virus only in single-cell seq, because mate read is barcode
#    filesize=$(stat -c%s "Chimeric.filterered.sam")
#    if [[ $filesize > 0 ]]; then
#        samtools view -bt /srv/disk00/cheyul1/Venus/STAR/indices/hg38SARSCoV2.genomeDir/hg38SARSCoV2.fna.fai -u Chimeric.filterered.sam | \
#        samtools sort -o Chimeric.filterered.sortedByCoord.bam
#        samtools index -b -@ 16 Chimeric.filterered.sortedByCoord.bam Chimeric.filterered.sortedByCoord.bam.bai
#    else
#        echo "No Chimeric Transcripts Found"
#    fi
#    unset filesize

    unset filenames
done
