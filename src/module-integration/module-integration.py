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
import pandas as pd
import pathlib


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

    parser.add_argument("--guideFASTA", type=str, required=True,
                        help="Fasta file for the guide sequence(s)")

    parser.add_argument("--virusChr", type=str, required=True,
                        help="Chromosome name for virus (e.g. NC_001802.1 for HIV)")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    parser.add_argument("--geneBed", type=str, required=True,
                        help="bed file to convert genomic coordinates to genes")

    parser.add_argument("--readFilesCommand", type=str, required=False,
                        help="uncompression command")

    parser.add_argument("--sensitivity", type=str, required=False, default="low",
                        help="sensitivity vs accuracy trade off, default option is low sensitivity, high accuracy")

    args = parser.parse_args()

    def quality_control():
        """Trims bad quality sequences"""
        pathlib.Path(args.out + "/quality_control/").mkdir(parents=True, exist_ok=True)

        cmd = "trim_galore " \
              + "--trim-n " \
              + "--quality 5 " \
              + "--phred33 " \
              + "--length 20 " \
              + "--output_dir " + args.out + "/quality_control/ " \
              + "--gzip "

        if int(args.thread) >= 7:
            cmd = cmd + "--cores 7 "
        else:
            cmd = cmd + "--cores " + args.thread + " "

        cmd = cmd + args.read

        return cmd

    def map_virus():
        """Maps to the virus genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/virus/ " \
              + "--genomeDir " + args.virusGenome + " "

        if args.readFilesCommand is not None:
            cmd = cmd + "--readFilesCommand " + args.readFilesCommand + " "

        cmd = cmd + "--readFilesIn " + args.out + "/quality_control/*_trimmed.fq.gz "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMax 1 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
                    + "--scoreDelBase -1 --scoreInsBase -1 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 "

        if args.sensitivity == "low":
            cmd = cmd \
                  + "--outFilterScoreMin 27 "
        else:
            cmd = cmd \
                  + "--outFilterScoreMin 25 "

        cmd = cmd \
              + "--outSAMtype BAM Unsorted "

        return cmd

    def map_hybrid():
        """Maps the viral mapped reads to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/hybrid/ " \
              + "--genomeDir " + args.hybridGenome + " "

        cmd = cmd + "--readFilesIn " + args.out + "/virus/Aligned.out.fastq "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMin 20 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
                    + "--scoreDelBase -1 --scoreInsBase -1 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 "

        if args.sensitivity == "low":
            cmd = cmd \
                  + "--outFilterScoreMin 27 "
        else:
            cmd = cmd \
                  + "--outFilterScoreMin 25 "

        cmd = cmd \
              + "--chimOutType WithinBAM SoftClip --chimSegmentMin 12 --chimJunctionOverhangMin 8 " \
              + "--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 " \
                    + "--chimNonchimScoreDropMin 10 " \
              + "--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 " \
                    + "--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 " \
              + "--alignInsertionFlush Right "

        cmd = cmd \
              + "--outSAMtype BAM Unsorted "

        return cmd

    def map_hybrid_junction():
        """Maps the viral mapped reads to the human genome"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/hybrid/junction/ " \
              + "--genomeDir " + args.hybridGenome + " "

        cmd = cmd + "--readFilesIn " + args.out + "/virus/Aligned.out.fastq "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMin 20 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 --scoreDelOpen 0 --scoreInsOpen 0 " \
                    + "--scoreDelBase -1 --scoreInsBase -1 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 "

        if args.sensitivity == "low":
            cmd = cmd \
                  + "--outFilterScoreMin 27 "
        else:
            cmd = cmd \
                  + "--outFilterScoreMin 25 "

        cmd = cmd \
              + "--chimOutType Junctions --chimSegmentMin 12 --chimJunctionOverhangMin 8 " \
              + "--chimMultimapNmax 10 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 0 " \
                    + "--chimNonchimScoreDropMin 10 " \
              + "--alignSJDBoverhangMin 10 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 " \
                    + "--alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 " \
              + "--alignInsertionFlush Right "

        cmd = cmd \
              + "--outSAMtype BAM Unsorted "

        return cmd

    def index_guide():
        """Generate the small genome index for guide sequence"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/guide/index/ " \
              + "--runMode genomeGenerate " \
              + "--genomeDir " + args.out + "/guide/index/ " \
              + "--genomeFastaFiles " + args.guideFASTA + " " \
              + "--genomeSAindexNbases 2 "

        return cmd

    def map_guide():
        """Classifies the integration sites based on guide sequences"""
        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/guide/ " \
              + "--genomeDir " + args.out + "/guide/index/ "

        cmd = cmd + "--readFilesIn " + args.out + "/hybrid/hybrid.out.fastq "

        cmd = cmd \
              + "--outSAMmultNmax 1 " \
              + "--alignIntronMax 1 --winBinNbits 7 " \
              + "--scoreGenomicLengthLog2scale 0 " \
              + "--outFilterMultimapNmax -1 --winAnchorMultimapNmax 3000 --seedPerReadNmax 30000 " \
                    + "--alignWindowsPerReadNmax 30000 --seedPerWindowNmax 1000 " \
              + "--outFilterMatchNmin 25 --outFilterMatchNminOverLread 0 " \
                    + "--outFilterScoreMin 20 --outFilterScoreMinOverLread 0 "

        cmd = cmd \
              + "--outReadsUnmapped None "

        return cmd

    def collect_gene():
        """Converts genomic coordinates into genes for final output"""
        with open(args.out + "/guide/Aligned.out.sam", "r") as file:
            quality_reads = []
            for line in file:
                if not line.startswith("@"):
                    quality_reads.append(line.split("\t")[0])
                    # print(line.split("\t")[0])
            quality_reads = list(set(quality_reads))

        integration_sites = pd.read_csv(args.out + "/hybrid/junction/Chimeric.out.junction", sep="\t")
        integrA = integration_sites[(integration_sites["read_name"].isin(quality_reads)) &
                                    (integration_sites["chr_donorA"] != args.virusChr)][["chr_donorA",
                                                                                         "brkpt_donorA",
                                                                                         "read_name"]]

        integrA.rename(columns={'chr_donorA': 'chr', 'brkpt_donorA': 'coord'}, inplace=True)
        integrB = integration_sites[(integration_sites["read_name"].isin(quality_reads)) &
                                    (integration_sites["chr_acceptorB"] != args.virusChr)][["chr_acceptorB",
                                                                                            "brkpt_acceptorB",
                                                                                            "read_name"]]
        integrB.rename(columns={'chr_acceptorB': 'chr', 'brkpt_acceptorB': 'coord'}, inplace=True)
        integration = pd.concat([integrA, integrB])

        integration.replace(["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
                             "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
                             "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
                             "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
                             "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
                             "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10"],
                            ["chr1", "chr2", "chr3", "chr4",
                             "chr5", "chr6", "chr7", "chr8",
                             "chr9", "chr10", "chr11", "chr12",
                             "chr13", "chr14", "chr15", "chr16",
                             "chr17", "chr18", "chr19", "chr20",
                             "chr21", "chr22", "chrX", "chrY"], inplace=True)
        integration["coord"] = integration["coord"].astype(int)

        gene_bed = pd.read_csv(args.geneBed, sep="\t", names=["chr", "start", "end", "gene"])
        gene_bed["start"] = gene_bed["start"].astype(int)
        gene_bed["end"] = gene_bed["end"].astype(int)

        genes = []
        for index, row in integration.iterrows():
            candidate = gene_bed[(row["chr"] == gene_bed["chr"]) &
                                 (row["coord"] >= gene_bed["start"]) &
                                 (row["coord"] <= gene_bed["end"])]
            if candidate.empty:
                genes.append("")
            elif candidate.shape[0] == 1:
                genes.append(candidate["gene"].to_string(index=False))
            else:
                gene_concat = ""
                for index2, gene in candidate["gene"].items():
                    gene_concat = gene_concat + gene + " "
                genes.append(gene_concat[:-1])
        integration["gene"] = genes
        class1 = integration[integration["gene"] != ""]
        class2 = integration[integration["gene"] == ""]
        # class1.sort_values(by=['chr', 'coord'], ascending=True, inplace=True)
        class1.reset_index(inplace=True, drop=True)
        class1.to_csv(path_or_buf=args.out + "/class1_integration.tsv", sep="\t", index=False)
        # class2.sort_values(by=['chr', 'coord'], ascending=True, inplace=True)
        class2.reset_index(inplace=True, drop=True)
        class2.to_csv(path_or_buf=args.out + "/class2_integration.tsv", sep="\t", index=False)

        return

    def action():
        # quality control
        os.system(quality_control())

        # virus (target) mapping
        os.system(map_virus())

        # prep input for hybrid mapping
        os.system("samtools fastq -@ " + args.thread + " -F 1 " + args.out + "/virus/Aligned.out.bam "
                  + "> " + args.out + "/virus/Aligned.out.fastq")

        # hybrid mapping
        os.system(map_hybrid())
        os.system("samtools sort -@ " + args.thread + " -o " + args.out + "/visual.bam "
                  + args.out + "/hybrid/Aligned.out.bam ")
        os.system("samtools index -@ " + args.thread + " -b " + args.out + "/visual.bam "
                  + args.out + "/visual.bam.bai ")
        os.system(map_hybrid_junction())

        # prep input for guide mapping
        os.system("samtools view -h -o " + args.out + "/hybrid/Aligned.out.sam "
                  + args.out + "/hybrid/Aligned.out.bam")
        outfile = open(args.out + "/hybrid/hybrid.out.sam", "w")
        with open(args.out + "/hybrid/Aligned.out.sam", "r") as file:
            for line in file:
                if line.startswith("@"):
                    outfile.write(line)
                elif line.split("\t")[2] != args.virusChr:
                    outfile.write(line)
        outfile.close()
        os.system("samtools view -S -b " + args.out + "/hybrid/hybrid.out.sam > "
                  + args.out + "/hybrid/hybrid.out.bam")
        os.system("samtools fastq -@ " + args.thread + " -F 1 " + args.out + "/hybrid/hybrid.out.bam > "
                  + args.out + "/hybrid/hybrid.out.fastq")

        # guide mapping
        os.system(index_guide())
        os.system(map_guide())

        # collect genes from coordinates
        collect_gene()

    ##################################################################################################
    # Action Steps
    ##################################################################################################
    if len(args.read) != 1:
        args.out_old = args.out
        args.read_old = args.read
        args.out = args.out_old + "/R1"
        args.read = args.read_old[0]
        action()
        args.out = args.out_old + "/R2"
        args.read = args.read_old[1]
        action()
    else:
        args.read = args.read[0]
        action()

    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
