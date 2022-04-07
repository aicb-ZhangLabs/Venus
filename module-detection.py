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
                        help="read of RNA-seq")

    parser.add_argument("--virusGenome", type=str, required=True,
                        help="directory path of virus genome index")

    parser.add_argument("--humanGenome", type=str, required=True,
                        help="directory path of human genome index")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--virusThreshold", type=str, required=False, default=0,
                        help="viral load threshold to filter out negligible viruses")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="uncompression command")

    args = parser.parse_args()

    # reads = args.read.split(" ")
    if len(args.read) == 1:
        seq = "single"
    elif len(args.read) == 2:
        seq = "paired"
    else:
        print("Unable to determine sequencing read type!!!")

    def map_human():
        """Maps to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/human/ " \
              + "--genomeDir " + args.humanGenome + " " \
              + "--outReadsUnmapped Fastx " \
              + "--outSAMtype None "

        if args.readFilesCommand is not None:
            cmd = cmd + "--readFilesCommand " + args.readFilesCommand + " "

        if seq == "single":
            cmd = cmd + "--readFilesIn " + args.read[0] + " "
        elif seq == "paired":
            cmd = cmd + "--readFilesIn " + args.read[0] + " " + args.read[1] + " "

        return cmd

    def map_virus():
        """Maps the leftover human reads to the virus genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " " \
              + "--outFilterMultimapNmax 1 " \
              + "--outSAMtype SAM "

        if seq == "single":
            cmd = cmd + "--readFilesIn " + args.out + "/human/Unmapped.out.mate1.fastq "
        elif seq == "paired":
            cmd = cmd + "--readFilesIn " + args.out + "/human/Unmapped.out.mate1.fastq " \
                                         + args.out + "/human/Unmapped.out.mate2.fastq "

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
        reference = pd.read_csv('/srv/disk00/cheyul1/Venus/repo/new_virus.species.txt', sep='\t', names=["id", "name"])
        reference = reference.set_index("id")

        # # Testing
        # stopper = 0
        for species in virus_species:
            count = 0
            with open(args.out + "/virus/Aligned.out.sam", "r") as file:
                for line in file:
                    if (not line.startswith("@")) and species in line:
                        count += 1
            species_count.append(count)
            species_name.append(reference.loc[species, "name"])
            species_ratio.append(count / total * 100.0)

            # # Testing
            # if stopper >= 5:
            #     exit
            # stopper += 1
            # print(species)
            # print(species_name)
            # print(species_count)
            # print(species_ratio)

        output = pd.DataFrame({"Name": species_name,
                               "Count": species_count,
                               "Percentage": species_ratio})
        output.sort_values(by="Percentage", ascending=False, inplace=True)
        output.to_csv(path_or_buf=args.out + "/detection_output.csv", index=False)

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

    os.system(map_human())  # map to human
    if seq == "single":
        os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")  # prep input
    elif seq == "paired":
        os.rename(args.out + "/human/Unmapped.out.mate1", args.out + "/human/Unmapped.out.mate1.fastq")  # prep input
        os.rename(args.out + "/human/Unmapped.out.mate2", args.out + "/human/Unmapped.out.mate2.fastq")
    os.system(map_virus())  # map to virus
    # pysam.view("-h", "-o", args.out + "/virus/Aligned.out.sam", args.out + "/virus/Aligned.out.bam")  # prep input
    output_infection()  # makes detection output file
    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

