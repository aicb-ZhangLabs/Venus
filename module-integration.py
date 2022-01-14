# ----------------------------------------------------------------------------------------------------------------------
# Venus: Virus infection detection and integration site discovery method using single-cell RNA-seq
# Integration module
#
# (C) 2022 Che Yu Lee, Irvine, California
# Released under GNU Public License (GPL)
# email cheyul1@uci.edu
# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# Import Libraries
########################################################################################################################
import os


########################################################################################################################
# Base scripts
########################################################################################################################
def script_base():
    """Base script that should be included in all Venus mappings"""
    script = "STAR " \
             + "--runThreadN 64 "
    return script


def script_base_singlecell():
    """Base script that should be included in all Venus bulk mappings"""
    script = "--soloCBwhitelist ${indices_dir}/737K-august-2016.txt " \
             + "--soloCBstart 1 " \
             + "--soloCBlen 16 " \
             + "--soloUMIstart 17 " \
             + "--soloUMIlen 10 " \
             + "--soloBarcodeReadLength 0 "
    return script


########################################################################################################################
# Extra Single-cell options
########################################################################################################################
def script_singlecell_gtf():
    """Script for mapping single-cell data with annotation gtf file in index"""
    script = "--soloType CB_UMI_Simple " \
             + "--outSAMattributes NH HI nM AS GX GN CR CB CY UR UB UY sS sQ sM "
    return script


def script_singlecell_nogtf():
    """Script for mapping single-cell data without annotation gtf file in index"""
    script = "--soloType CB_samTagOut " \
             + "--soloCBmatchWLtype 1MM " \
             + "--outSAMattributes NH HI nM AS CB UR "
    return script


########################################################################################################################
# Choosing which genome index to map to
########################################################################################################################
def script_targetvirus():
    """Script for mapping to the target virus genome index"""
    script = "--genomeDir ${indices_dir}/HIV.genomeDir/ " \
             + "--readFilesIn ${run_dir}/R1/trim_fastq/*_trimmed.fq.gz " \
             + "--outFileNamePrefix ${run_dir}/R1/aln_fastq/ " \
             + "--outFilterMultimapNmax 1 " \
             + "--outSAMtype BAM Unsorted "
    return script


########################################################################################################################
# Choosing which genome index to map to
########################################################################################################################
def script_sensitivity_shortreads():
    """Script for increasing sensitivity in detecting short reads"""
    script = "--outSAMmultNmax 1 " \
             + "--alignIntronMax 1 --winBinNbits 7 " \
             + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 " \
             + "--scoreInsOpen 0 --scoreDelBase -1 --scoreInsBase -1 " \
             + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 " \
             + "--seedPerReadNmax 30000 --alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
             + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --outFilterScoreMin 27 "
    return script


########################################################################################################################
# Combine scripts to command
########################################################################################################################
def script_combine():
    cmd = script_base() + script_targetvirus()
    return cmd


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # os.system(script_combine())
    print(script_combine())
