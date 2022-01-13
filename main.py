import argparse
import os
import sys
import csv

### Currently, the indices are un-changeable; but may be user-supplied in future ###
# This is the current path to the default human genome index directory
human_indexDir = "/srv/disk00/cheyul1/excessSTAR-2.7.8a/indices/human.genomeDir/"
# This is the current path to the default mega virus index directory
virus_indexDir = "/srv/disk00/cheyul1/excessSTAR-2.7.8a/indices/virus.genomeDir/"
# This is the current path to the virus species metatable
virus_species_metatable = "/srv/disk00/cheyul1/excessSTAR-2.7.8a/virus.species.txt"

def get_count(line):
    """
    Obtains the read counts from the initially unsorted venus.out.tsv line.
    """
    line_fields = line.strip().split("\t")
    return int(line_fields[1])


### Creates the Argument Parser object ###
parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " + \
                                             "Virus dEtecting in humaN bUlk and Single cell rna sequencing")


### Specifies the argument for bulk RNA-seq ###
parser.add_argument("--read1", type=str, required=True,
                    help="read1 of RNA-seq (barcode if single-cell)")
parser.add_argument("--read2", type=str, required=False,
                    help="read2 of RNA-seq (cDNA if single-cell)")
# parser.add_argument("--indexDir", type=str, required=False, default=human_indexDir, 
#                     help="user-specified genome index directory")
parser.add_argument("--outDir", type=str, required=False, default=os.getcwd(), 
                    help="directory to store output")
parser.add_argument("--singleCBstart", type=int, required=False, 
                    help="cell barcode's start position")
parser.add_argument("--singleCBlen", type=int, required=False,
                    help="cell barcode's length")
parser.add_argument("--singleUMIstart", type=int, required=False, 
                    help="UMI's start position")
parser.add_argument("--singleUMIlen", type=int, required=False, 
                    help="UMI's length")
parser.add_argument("--singleWhitelist", type=str, required=False, 
                    help="single-cell barcode whitelist")
args = parser.parse_args()


### Creates a python output file ###


### Organizes the output directories for human then virus mappings ###
try:
    os.chdir(args.outDir)
    sys.stdout = open("venus.log", "x")
    print("Current working directory: {0}".format(os.getcwd()))
except:
    print("An error has occured while changing into '{0}' directory".format(args.outDir))
    sys.exit()

try:
    human_outDir = os.getcwd() + "/human/" 
    os.mkdir(human_outDir)
    virus_outDir = os.getcwd() + "/virus/" 
    os.mkdir(virus_outDir)
except:
    print("An error has occured while making directories in '{0}' directory".format(args.outDir))
    sys.exit()



### Firstly maps the RNA-seq reads to the human genome ###

if args.singleCBstart:  # Single-Cell RNA-seq
    cmd="STAR " + \
        "--runThreadN 16 " + \
        "--outFileNamePrefix " + human_outDir + " " \
        "--genomeDir " + human_indexDir + " " \
        "--readFilesIn " + args.read2 + " " + args.read1 + " " \
        "--outReadsUnmapped Fastx " + \
        "--outSAMtype None" + " " \
        "--soloType CB_UMI_Simple" + " " \
        "--soloCBwhitelist " + str(args.singleWhitelist) + " " \
        "--soloCBstart " + str(args.singleCBstart) + " " \
        "--soloCBlen " + str(args.singleCBlen) + " " \
        "--soloUMIstart " + str(args.singleUMIstart) + " " \
        "--soloUMIlen " + str(args.singleUMIlen) + " " \
        "--soloBarcodeReadLength 0"

elif args.read2:  # Bulk Paired-end RNA-seq
        cmd="STAR " + \
            "--runThreadN 16 " + \
            "--outFileNamePrefix " + human_outDir + " " \
            "--genomeDir " + human_indexDir + " " \
            "--readFilesIn " + args.read1 + " " + args.read2 + " " \
            "--outReadsUnmapped Fastx " + \
            "--outSAMtype None"
else:   # Bulk Single-end RNA-seq
    cmd="STAR " + \
        "--runThreadN 16 " + \
        "--outFileNamePrefix " + human_outDir + " " \
        "--genomeDir " + human_indexDir + " " \
        "--readFilesIn " + args.read1 + " " \
        "--outReadsUnmapped Fastx " + \
        "--outSAMtype None"

os.system(cmd)  # Command run
print("Running " + cmd)




### Secondly maps the leftover reads to the viral genome ###

### Appropriately renames the read1 & read2 for the virus mapping ###
args.read1=human_outDir + "Unmapped.out.mate1.fastq"
os.rename(human_outDir + "Unmapped.out.mate1", args.read1)
if args.read2:
    args.read2=human_outDir + "Unmapped.out.mate2.fastq"
    os.rename(human_outDir + "Unmapped.out.mate2", args.read2)

if args.singleCBstart:  # Single-Cell RNA-seq
    cmd="STAR " + \
        "--runThreadN 16 " + \
        "--outFileNamePrefix " + virus_outDir + " " \
        "--genomeDir " + virus_indexDir + " " \
        "--readFilesIn " + args.read1 + " " + args.read2 + " " \
        "--outFilterMultimapNmax 1" + " " \
        "--soloType CB_UMI_Simple" + " " \
        "--soloCBwhitelist " + str(args.singleWhitelist) + " " \
        "--soloCBstart " + str(args.singleCBstart) + " " \
        "--soloCBlen " + str(args.singleCBlen) + " " \
        "--soloUMIstart " + str(args.singleUMIstart) + " " \
        "--soloUMIlen " + str(args.singleUMIlen) + " " \
        "--soloBarcodeReadLength 0"

elif args.read2:  # Bulk Paired-end RNA-seq
        cmd="STAR " + \
            "--runThreadN 16 " + \
            "--outFileNamePrefix " + virus_outDir + " " \
            "--genomeDir " + virus_indexDir + " " \
            "--readFilesIn " + args.read1 + " " + args.read2 + " " \
            "--outFilterMultimapNmax 1"

else:   # Bulk Single-end RNA seq
    cmd="STAR " + \
        "--runThreadN 16 " + \
        "--outFileNamePrefix " + virus_outDir + " " \
        "--genomeDir " + virus_indexDir + " " \
        "--readFilesIn " + args.read1 + " " \
        "--outFilterMultimapNmax 1"
os.system(cmd)  # Command run
print("Running " + cmd)



# 1. Creates the species_reads txt file.
#    This will be used to count the total reads
#    and to create the unique species txt file. 
os.chdir(args.outDir)
with open("species_reads.txt", "x") as species_reads:
    # with open("Aligned.out.sam") as aligned:                                    # For local testing
    with open(os.getcwd() + "/virus/" + "Aligned.out.sam", "r") as aligned:
        reads_total = 0
        for line in aligned.readlines():
            if "@" not in line:
                species_id = line.strip().split("\t")[2]
                species_reads.write(species_id)
                species_reads.write("\n")
                reads_total += 1


# 2. Creates the species txt file.
#    This will be used to loop over the unique species
#    thus counting the reads and writing the Venus output
with open("species_reads.txt") as species_reads: 
    with open("species.txt", "x") as species:
        species.writelines(set(species_reads.readlines()))


# 3. Creates the venus output tsv file.
#    Per row/specie, there will be:
#    specie name | reads count | reads total | reads % of total
venus = open("venus.out.tsv", "x")
venus.close()
with open("venus.out.tsv", "a") as venus:
    with open("species.txt") as species:
        for specie in species.readlines():
            specie = specie.replace("\n", "")
            with open(virus_species_metatable) as metatable:
                metareader = csv.reader(metatable, delimiter="\t")
                for row in metareader:
                    if specie == row[0]:
                        venus.write(row[1] + "\t")                                  # writes the specie name
            with open("species_reads.txt") as species_reads:
                reads_count = 0
                for line in species_reads.readlines():
                    if specie == line.replace("\n", ""):
                        reads_count += 1
                venus.write(str(reads_count) + "\t")                                # writes the reads count
            venus.write(str(reads_total) + "\t")                                    # writes the reads total
            venus.write(str(reads_count / (reads_total * 1.0)  * 100) + "%")        # writes the reads % of total
            venus.write("\n")



# 4. Sorts the venus output tsv file.
#    It first reads the lines to sort,
#    then it overwrites the venus output.
with open("venus.out.tsv") as venus:
    venus_lines = venus.readlines()
    venus_lines.sort(key=get_count, reverse=True)

with open("venus.out.tsv", "w") as venus:
    for line in venus_lines:
        venus.write(line)

os.remove(os.getcwd() + "/species_reads.txt")
os.remove(os.getcwd() + "/species.txt")

sys.stdout.close()


### local tests ###
# print("{} is read1 and {} is read2; {} is the genomeDir".format(args.read1, args.read2, args.genomeDir))
