for name in $(grep "^SRR" outputs/10-21-21/bowtie24/R1/Aligned.out.sam | cut -f 1); do
    if grep -Fq "$name" outputs/10-21-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.hisat.sam; then
        echo "found!"
    else
        echo $name
    fi
done

#######################################################################################################################################
# HISAT reads missing in STAR spliced
for name in $(grep "^SRR" v | cut -f 1); do
    if ! grep -Fq "$name" outputs/10-19-21/WGsplice-a/GSE126230_run35/splice_human_R1/output/Aligned.out.sam; then
        echo $name
        grep "$name" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.hisat.sam >> missed_reads.txt
        echo -en "\n" >> missed_reads.txt
#    else
#        echo "found!"
    fi
done


# STAR spliced reads in extra to HISAT
for name in $(grep "^SRR" outputs/10-19-21/WGsplice-a/GSE126230_run35/splice_human_R1/output/Aligned.out.sam | cut -f 1); do
    if ! grep -Fq "$name" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.hisat.sam; then
        echo $name
        grep "$name" outputs/10-19-21/WGsplice-a/GSE126230_run35/splice_human_R1/output/Aligned.out.sam >> extra_reads.txt
        echo -en "\n" >> extra_reads.txt
#    else
#        echo "found!"
    fi
done


# Reads found in both HISAT and STAR spliced
for name in $(grep "^SRR" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.hisat.sam | cut -f 1); do
    if grep -Fq "$name" outputs/10-19-21/WGsplice-a/GSE126230_run35/splice_human_R1/output/Aligned.out.sam; then
        echo $name
        grep "$name" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.hisat.sam >> both_reads.txt
        echo -en "\n" >> both_reads.txt
#    else
#        echo "not found!"
    fi
done


#######################################################################################################################################
# bowtie reads missing in STAR <--- No Bowtie reads missing
for name in $(grep "^SRR" /srv/disk00/cheyul1/outputs/09-24-21/GSE112576-run20/GSE112576-run20_R1/tmp/alns/SRR6944369.1.hiv.bowtie.sam | cut -f 1); do
    if ! grep -Fq "$name" /srv/disk00/cheyul1/outputs/10-27-21/chimr2-star/GSE112576-run20/R1/aln_fastq/Aligned.out.sam; then
        echo $name
        grep "$name" /srv/disk00/cheyul1/outputs/09-24-21/GSE112576-run20/GSE112576-run20_R1/tmp/alns/SRR6944369.1.hiv.bowtie.sam >> dna_missed_reads.txt
        echo -en "\n" >> dna_missed_reads.txt
    else
        echo "found!"
    fi
done


# STAR reads in extra to bowtie
for name in $(grep "^SRR" outputs/10-19-21/WGsplice-a/GSE126230_run35/human_R1/output/Aligned.out.sam | cut -f 1); do
    if ! grep -Fq "$name" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.bowtie.sam; then
        echo $name
        grep "$name" outputs/10-19-21/WGsplice-a/GSE126230_run35/human_R1/output/Aligned.out.sam >> dna_extra_reads.txt
        echo -en "\n" >> dna_extra_reads.txt
    else
        echo "found!"
    fi
done


# Reads found in both bowtie and STAR
for name in $(grep "^SRR" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.bowtie.sam | cut -f 1); do
    if grep -Fq "$name" outputs/10-19-21/WGsplice-a/GSE126230_run35/human_R1/output/Aligned.out.sam; then
        echo $name
        grep "$name" outputs/10-19-21/GSE126230-run35/GSE126230-run35_R1/tmp/alns/SRR8545436.1.hum.bowtie.sam >> dna_both_reads.txt
        echo -en "\n" >> dna_both_reads.txt
    else
        echo "not found!"
    fi
done

