# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/fig%202.png)

## Setup
This is how to set up Venus.

One can download all the virus fasta files from NCBI as so:
```
mkdir tmp  &&  cd tmp
wget ftp://ftp.ncbi.nlm.nih.gov//genomes/Viruses/all.fna.tar.gz
tar -xvf all.fna.tar.gz
cat */*.fna | sed "s/>/>Virus:/" > ../Virus.fa
cd ..
rm -r tmp
```

## Creating Index
This step creates the necessary human and viral indices directories in Venus. Of note, since some viruses may lack a gtf annotation file or biological splicing mechanisms, use of a vurs gtf file is optional.

To create indices:
```   
python3 module-index.py \
    --humanGenome path/to/human.genomeDir \
    --humanFASTA human.fna \
    --humanGTF human.gtf \
    --virusGenome path/to/virus.genomeDir \
    --virusFASTA mega-virus.fa \
    --out path/to/output/dir \
    --thread 32
```

## Virus Detection Module
This module detects viral load and will output a list of infecting viral species or infected cell barcodes, depending on the input and the viral index used. (*Note: For path/to/output/dir parameter, please do not include an end '/'.*) 

For bulk single-end sequencing:
```
python3 module-detection.py \
    --read SRR6944349.1_1.fastq.gz \
    --virusGenome path/to/virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32
```

## Integration Site Discovery Module
This is how to use Venus's integration site discovery module.
