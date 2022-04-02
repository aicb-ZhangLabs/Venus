#!/bin/bash
#SBATCH --job-name=integrSeq-QS-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=zhanglab.p
#SBATCH --output=/srv/disk00/cheyul1/logs/03-09-22/integrSeq-QS-all.log

run_prefix=GSE112576
out_dir=/srv/disk00/cheyul1/outputs/03-09-22/integrSeq-QS-all
STAR_dir=/srv/disk00/cheyul1/Venus/STAR
indices_dir=${STAR_dir}/indices

for i in {0..39}; do
    run_dir=${out_dir}/${run_prefix}-run${i}
    
    for j in {1..2}; do
        mkdir -p ${run_dir}/R${j}/input
        awk -F '\t' '$3 != "NC_001802.1"' /srv/disk00/cheyul1/outputs/11-10-21/fsn-bam/GSE112576-run${i}/R${j}/Aligned.out.sam \
        >> ${run_dir}/R${j}/input/hybrid_R${j}.sam
        samtools view -S -b ${run_dir}/R${j}/input/hybrid_R${j}.sam > ${run_dir}/R${j}/input/hybrid_R${j}.bam
        samtools fastq -@ 28 -F 1 ${run_dir}/R${j}/input/hybrid_R${j}.bam \
        > ${run_dir}/R${j}/input/hybrid_R${j}.fastq
        
        STAR \
        --runThreadN 28 \
        --genomeDir ${indices_dir}/integrSeq.genomeDir/ \
        --outFileNamePrefix ${run_dir}/R${j}/ \
        --readFilesIn ${run_dir}/R${j}/input/hybrid_R${j}.fastq \
        --outReadsUnmapped None --outSAMmultNmax 1 \
        --outFilterScoreMin 20 --outFilterScoreMinOverLread 0 --outFilterMatchNmin 25 --outFilterMatchNminOverLread 0 \
        --alignIntronMax 1 --scoreGenomicLengthLog2scale 0 --winBinNbits 7 \
        --outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000
        
        grep -v "@" ${run_dir}/R${j}/Aligned.out.sam | cut -f 1 > ${run_dir}/R${j}/class1.reads
        
        # Genes
        while read read_name
        do
            grep ${read_name} /srv/disk00/cheyul1/outputs/11-10-21/fsn-junction/GSE112576-run${i}/R${j}/Chimeric.out.junction | \
            grep "NC_001802.1" | cut -f 1,2 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
            
            grep ${read_name} /srv/disk00/cheyul1/outputs/11-10-21/fsn-junction/GSE112576-run${i}/R${j}/Chimeric.out.junction | \
            grep "NC_001802.1" | cut -f 4,5 | grep -v "NC_001802.1" | sort | uniq >> ${out_dir}/integration.junction
        done < ${run_dir}/R${j}/class1.reads
        
        # BAM file for Visual in IGV
        while read read_name
        do
            grep ${read_name} /srv/disk00/cheyul1/outputs/11-10-21/fsn-bam/GSE112576-run${i}/R${j}/Aligned.out.sam | \
            awk -F '\t' '$3 != "NC_001802.1"' >> ${out_dir}/integration.sam
#            echo ${read_name} "done!"
        done < ${run_dir}/R${j}/class1.reads
        
    done

done



########################################################################################################
######================================ Collects integration genes ================================######
########################################################################################################
grep "NC_" ${out_dir}/integration.junction | sort -t $'\t' -k 1,1 -k 2,2 -V | uniq | \
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



###########################################################################################################
######================================ Prepares visualization in IGV ================================######
###########################################################################################################
samtools view -@ 28 -bt /srv/disk00/cheyul1/Venus/STAR/indices/hg38HIV.genomeDir/hg38HIV.fna.fai -u ${out_dir}/integration.sam | \
samtools sort -@ 28 -o ${out_dir}/integration.Sorted.bam
samtools index -@ 28 -b ${out_dir}/integration.Sorted.bam ${out_dir}/integration.Sorted.bam.bai
