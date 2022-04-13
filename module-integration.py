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
import argparse
import os
import pysam


########################################################################################################################
# Main
########################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " + \
                                                 "Virus dEtecting in humaN bUlk and Single cell rna sequencing")

    parser.add_argument("--read", type=str, required=True, nargs='+',
                        help="read of RNA-seq \n(single-cell) first read should be cDNA, second should be CB+UMI")

    parser.add_argument("--virusGenome", type=str, required=True,
                        help="directory path of virus genome index")

    parser.add_argument("--hybridGenome", type=str, required=True,
                        help="directory path of hybrid genome index")

    parser.add_argument("--guideGenome", type=str, required=False,
                        help="directory path of guide sequence index")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="uncompression command")

    args = parser.parse_args()

    def map_virus():
        """Maps to the virus genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " "

        if args.readFilesCommand is not None:
            cmd = cmd + "--readFilesCommand " + args.readFilesCommand + " "

        return cmd

    def map_hybrid():
        """Maps the viral mapped reads to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/human/ " \
              + "--genomeDir " + args.humanGenome + " "

        return cmd

    def map_guide():
        """Classifies the integration sites based on guide sequences"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " " \
              + "--outFilterMultimapNmax 1 "
        
        return cmd

    ##################################################################################################
    # Action Steps
    ##################################################################################################
    print(map_virus)
    print(map_hybrid())
    print(map_guide())
    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
