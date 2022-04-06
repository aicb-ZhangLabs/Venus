# ----------------------------------------------------------------------------------------------------------------------
# Venus: Virus infection detection and integration site discovery method using single-cell RNA-seq
# Detection module
#
# (C) 2022 Che Yu Lee, Irvine, California
# Released under GNU Public License (GPL)
# email cheyul1@uci.edu
#
# ----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
# Import Libraries
########################################################################################################################
import argparse
import os
import pysam


########################################################################################################################
# Main
########################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " + \
                                                 "Virus dEtecting in humaN bUlk and Single cell rna sequencing")

    parser.add_argument("--read", type=str, required=True,
                        help="read of RNA-seq")

    parser.add_argument("--virusGenome", type=str, required=True,
                        help="directory of virus genome index")

    parser.add_argument("--humanGenome", type=str, required=True,
                        help="directory of human genome index")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory of human genome index")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="number of parallel threads")

    args = parser.parse_args()

    def map_human():
        if args.readFilesCommand is None:
            cmd = "STAR " \
                  + "--runThreadN " + args.thread + " " \
                  + "--outFileNamePrefix " + args.out + "/human " \
                  + "--genomeDir " + args.humanGenome + " " \
                  + "--readFilesIn " + args.read + " " \
                  + "--outReadsUnmapped Fastx " \
                  + "--outSAMtype None"
        else:
            cmd = "STAR " \
                  + "--runThreadN " + args.thread + " " \
                  + "--readFilesCommand " + args.readFilesCommand + " " \
                  + "--outFileNamePrefix " + args.out + "/human " \
                  + "--genomeDir " + args.humanGenome + " " \
                  + "--readFilesIn " + args.read + " " \
                  + "--outReadsUnmapped Fastx " \
                  + "--outSAMtype None"

        # # Testing
        # f = open(args.out + "/human/Unmapped.out.mate1", "a")
        # f.close()

        return cmd

    def map_virus():
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus " \
              + "--genomeDir " + args.virusGenome + " " \
              + "--readFilesIn " + args.out + "/human/Unmapped.out.mate1.fastq" + " " \
              + "--outFilterMultimapNmax 1 " \
              + "--outSAMtype SAM"

        return cmd

    # # Testing
    # print(map_human())
    # os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")
    # print(map_virus())

    os.system(map_human())
    os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")
    os.system(map_virus())
    # pysam.view("-h", "-o",
    #            args.out + "/virus/Aligned.sortedByCoord.sam",
    #            args.out + "/virus/Aligned.sortedByCoord.bam")
    return


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
# Choosing which genome index to map to
########################################################################################################################
def script_human():
    """Script for mapping to the human genome index"""
    script = "--genomeDir ${indices_dir}/human.genomeDir/ " \
             + "--readFilesIn ${fastqs}/${filenames[0]} ${fastqs}/${filenames[1]} " \
             + "--outFileNamePrefix ${human_dir}/output/ " \
             + "--outReadsUnmapped Fastx " \
             + "--outSAMtype None "
    return script


def script_targetvirus():
    """Script for mapping to the target virus genome index"""
    script = "--genomeDir ${indices_dir}/HIV.genomeDir/ " \
             + "--readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq " \
             + "--outFileNamePrefix ${virus_dir}/output/ " \
             + "--outFilterMultimapNmax 1 " \
             + "--outSAMtype BAM SortedByCoordinate "
    return script


def script_megavirus():
    """Script for mapping to the mega virus genome index"""
    script = "--genomeDir ${indices_dir}/new_virus.genomeDir/ " \
             + "--readFilesIn ${fastqs}/Unmapped.out.mate1.fastq ${fastqs}/Unmapped.out.mate2.fastq " \
             + "--outFileNamePrefix ${mega_virus_dir}/output/ " \
             + "--outFilterMultimapNmax 1 " \
             + "--outSAMtype BAM SortedByCoordinate "
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
# Combine scripts to command
########################################################################################################################
def script_combine():
    cmd = script_base() + script_megavirus() + script_base_singlecell() + script_singlecell_nogtf()
    return cmd


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # os.system(script_combine())
    # print(script_combine())
    main()

