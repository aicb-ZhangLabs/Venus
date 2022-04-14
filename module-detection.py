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
import pandas as pd


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

    parser.add_argument("--humanGenome", type=str, required=True,
                        help="directory path of human genome index")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--virusThreshold", type=str, required=False, default=0,
                        help="viral load threshold to filter out negligible viruses")

    parser.add_argument("--virusChrRef", type=str, required=True,
                        help="tsv file to map NC_* id to virus species name, e.g. NC_001802.1 --> HIV-1")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="uncompression command")

    parser.add_argument("--singleCellBarcode", type=str, required=False, nargs=2,
                        help="(single-cell) barcode specifications, e.g. '1 16' = start at 1, length is 16")

    parser.add_argument("--singleUniqueMolIdent", type=str, required=False, nargs=2,
                        help="(single-cell) umi specifications, e.g. '17 10' = start at 17, length is 10")

    parser.add_argument("--singleWhitelist", type=str, required=False,
                        help="(single-cell) barcode whitelist")

    args = parser.parse_args()

    # Determine sequencing read type
    if len(args.read) == 1:
        read_type = "single_end"
    else:
        read_type = "paired_end"

    # Determine sequencing resolution
    if args.singleCellBarcode is not None:
        seq_resol = "single_cell"
    else:
        seq_resol = "bulk"

    def map_human():
        """Maps to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/human/ " \
              + "--genomeDir " + args.humanGenome + " "

        if args.readFilesCommand is not None:
            cmd = cmd + "--readFilesCommand " + args.readFilesCommand + " "

        if read_type == "single_end":
            cmd = cmd + "--readFilesIn " + args.read[0] + " "
        elif read_type == "paired_end":
            cmd = cmd + "--readFilesIn " + args.read[0] + " " + args.read[1] + " "

        if seq_resol == "single_cell":
            cmd = cmd \
                  + "--soloType CB_UMI_Simple " \
                  + "--soloCBwhitelist " + args.singleWhitelist + " " \
                  + "--soloBarcodeReadLength 0 " \
                  + "--soloCBstart " + args.singleCellBarcode[0] + " " \
                  + "--soloCBlen " + args.singleCellBarcode[1] + " " \
                  + "--soloUMIstart " + args.singleUniqueMolIdent[0] + " " \
                  + "--soloUMIlen " + args.singleUniqueMolIdent[1] + " "

        cmd = cmd \
              + "--outSAMtype None " \
              + "--outReadsUnmapped Fastx "

        return cmd

    def map_virus():
        """Maps the leftover human reads to the virus genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " " \
              + "--outFilterMultimapNmax 1 "

        if read_type == "single_end":
            cmd = cmd + "--readFilesIn " + args.out + "/human/Unmapped.out.mate1.fastq "
        elif read_type == "paired_end":
            cmd = cmd + "--readFilesIn " + args.out + "/human/Unmapped.out.mate1.fastq " \
                  + args.out + "/human/Unmapped.out.mate2.fastq "

        if seq_resol == "single_cell":
            cmd = cmd \
                  + "--soloType CB_samTagOut " \
                  + "--soloCBwhitelist " + args.singleWhitelist + " " \
                  + "--soloCBmatchWLtype 1MM " \
                  + "--soloBarcodeReadLength 0 " \
                  + "--soloCBstart " + args.singleCellBarcode[0] + " " \
                  + "--soloCBlen " + args.singleCellBarcode[1] + " " \
                  + "--soloUMIstart " + args.singleUniqueMolIdent[0] + " " \
                  + "--soloUMIlen " + args.singleUniqueMolIdent[1] + " " \
                  + "--outSAMtype BAM Unsorted " \
                  + "--outSAMattributes NH HI nM AS CR UR "
        elif seq_resol == "bulk":
            cmd = cmd \
                  + "--outSAMtype SAM "

        return cmd

    def output_infection():
        """Produces the detection output file for mega-virus mode"""
        with open(args.out + "/virus/Aligned.out.sam", "r") as file:
            count = 0
            virus_species = []
            for line in file:
                if not line.startswith("@"):
                    virus_species.append(line.split("\t")[2])
                    count += 1
            total = len(virus_species)
            virus_species = list(set(virus_species))

            # # Testing
            # print(len(virus_species))

        species_name = []
        species_count = []
        species_ratio = []
        species_barcodes = []
        reference = pd.read_csv(args.virusChrRef, sep='\t', names=["id", "name"])
        reference = reference.set_index("id")

        # # Testing
        # stopper = 0
        for species in virus_species:
            virus_count = 0
            virus_barcodes = []
            with open(args.out + "/virus/Aligned.out.sam", "r") as file:
                for line in file:
                    if (not line.startswith("@")) and species in line:
                        virus_count += 1
                        if seq_resol == "single_cell":
                            virus_barcodes.append(line.split("\t")[-2].split(sep=":")[-1])
            if virus_count >= int(args.virusThreshold):
                species_count.append(virus_count)
                species_name.append(reference.loc[species, "name"])
                species_ratio.append(virus_count / total * 100.0)
                species_barcodes.append(virus_barcodes)

            # # Testing
            # if stopper >= 5:
            #     exit
            # stopper += 1
            # print(species)
            # print(species_name)
            # print(species_count)
            # print(species_ratio)

        if seq_resol == "single_cell":
            output = pd.DataFrame({"Name": species_name,
                                   "Count": species_count,
                                   "Percentage": species_ratio,
                                   "Barcodes": species_barcodes})
        else:
            output = pd.DataFrame({"Name": species_name,
                                   "Count": species_count,
                                   "Percentage": species_ratio})
        output.sort_values(by="Percentage", ascending=False, inplace=True)
        output.to_csv(path_or_buf=args.out + "/detection_output.tsv", sep="\t", index=False)

        # # Testing
        # print(virus_species[:5])
        # print(species_name[:5])
        # print(species_count[:5])
        # print(len(species_count))
        # print(len(species_name))

        return output

    # Testing
    # print(map_human())
    # f = open(args.out + "/human/Unmapped.out.mate1", "a")
    # f.close()
    # os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")
    # print(map_virus())

    ##################################################################################################
    # Action Steps
    ##################################################################################################
    os.system(map_human())  # map to human

    # prep unmapped human reads
    if read_type == "single_end":
        os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")
    elif read_type == "paired_end":
        os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")
        os.rename(args.out + "/human/Unmapped.out.mate2", args.out + "/human/Unmapped.out.mate2.fastq")

    os.system(map_virus())  # map to virus

    # input conversion for single-cell
    if seq_resol == "single_cell":
        infile = pysam.AlignmentFile(args.out + "/virus/Aligned.out.bam", "rb")
        outfile = pysam.AlignmentFile(args.out + "/virus/Aligned.out.sam", "w", template=infile)
        for s in infile:
            outfile.write(s)
        # pysam.view("-h", "-o",
        #            args.out + "/virus/Aligned.out.sam",
        #            args.out + "/virus/Aligned.out.bam")

    output_infection()  # makes detection output file
    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
