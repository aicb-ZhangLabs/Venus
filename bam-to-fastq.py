import pysam
import os


def main():
    # pysam.sort("-o", "../data/conv_fastq-sorted.bam", "../data/conv_fastq.bam")
    # pysam.index("../data/conv_fastq-sorted.bam")
    # samfile_in = pysam.AlignmentFile("../data/conv_fastq-sorted.bam", "rb")
    #
    # with open("../data/conv_fastq.fastq", "w") as samfile_out:
    #     for read in samfile_in.fetch():
    #         if read.is_read1:
    #             samfile_out.write(str(read.qname) + "_1" + "\n" + str(read.query_sequence) +
    #                               "\n" + "+" + "\n" + str(read.qual) + "\n")
    #         elif read.is_read2:
    #             samfile_out.write(str(read.qname) + "_2" + "\n" + str(read.query_sequence) +
    #                               "\n" + "+" + "\n" + str(read.qual) + "\n")
    os.system("samtools fastq -@ 1 -F 1 ../data/conv_fastq.bam "
              + "> ../data/conv_fastq.fastq")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
