import argparse
import os
import pathlib


def main():
    parser = argparse.ArgumentParser(description="VENUS, a subtractive analysis software: " + \
                                                 "Virus dEtecting in humaN bUlk and Single cell rna sequencing")

    parser.add_argument("--humanGenome", type=str, required=True,
                        help="directory path of human genome index")

    parser.add_argument("--humanFASTA", type=str, required=True,
                        help="fasta file path for the human genome index")

    parser.add_argument("--humanGTF", type=str, required=True,
                        help="gtf file path for the human genome index")

    parser.add_argument("--virusGenome", type=str, required=True,
                        help="directory path of virus genome index")

    parser.add_argument("--virusFASTA", type=str, required=True,
                        help="fasta file path for the virus genome index")

    parser.add_argument("--virusGTF", type=str, required=False,
                        help="gtf file path for the virus genome index")

    parser.add_argument("--out", type=str, required=False, default=os.getcwd(),
                        help="directory path of output dir")

    parser.add_argument("--thread", type=str, required=False, default="1",
                        help="number of parallel threads")

    args = parser.parse_args()

    def index_human():
        pathlib.Path(args.humanGenome).mkdir(parents=True, exist_ok=True)

        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/ " \
              + "--runMode genomeGenerate " \
              + "--genomeDir " + args.humanGenome + " " \
              + "--genomeFastaFiles " + args.humanFASTA + " " \
              + "--genomeSAindexNbases 14 "

        if "gff" in args.humanGTF:
            cmd = cmd \
                  + "--sjdbGTFfile " + args.humanGTF + " " \
                  + "--sjdbGTFtagExonParentTranscript Parent "
        else:
            cmd = cmd \
                  + "--sjdbGTFfile " + args.humanGTF + " "

        return cmd

    def index_virus():
        pathlib.Path(args.virusGenome).mkdir(parents=True, exist_ok=True)

        cmd = "STAR " \
              + "--runThreadN " + args.thread + " " \
              + "--outFileNamePrefix " + args.out + "/ " \
              + "--runMode genomeGenerate " \
              + "--genomeDir " + args.virusGenome + " " \
              + "--genomeFastaFiles " + args.virusFASTA + " " \
              + "--genomeSAindexNbases 2 "

        if args.virusGTF is not None:
            if "gff" in args.virusGTF:
                cmd = cmd \
                      + "--sjdbGTFfile " + args.virusGTF + " " \
                      + "--sjdbGTFtagExonParentTranscript Parent "
            else:
                cmd = cmd  \
                      + "--sjdbGTFfile " + args.virusGTF + " "

        return cmd

    # Testing
    print(index_human())
    print(index_virus())

    # os.system(index_human())
    # os.system(index_virus())
    return


if __name__ == '__main__':
    main()

