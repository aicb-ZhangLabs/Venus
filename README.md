# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/fig%202.png)

## Setup
This is how to set up Venus.

Set up the conda envrionment. Note the spec-file.txt is for Linux platform.
```
conda create --name venus --file venus_spec-file.txt
```

One can make the mega-virus.fasta file from NCBI as so:
```
mkdir tmp  &&  cd tmp
wget ftp://ftp.ncbi.nlm.nih.gov//genomes/Viruses/all.fna.tar.gz
tar -xvf all.fna.tar.gz
cat */*.fna | sed "s/>/>Virus:/" > ../mega-virus.fasta
cd ..
rm -r tmp
```

## Creating Index
This step creates the necessary human and viral indices directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. (*Note: For path/to/output/dir parameter, please do not include an end '/'.*) 

To create indices (without a virus gtf):
```   
python3 module-index.py \
    --humanGenome path/to/human.genomeDir \
    --humanFASTA human.fasta \
    --humanGTF human.gtf \
    --virusGenome path/to/virus.genomeDir \
    --virusFASTA mega-virus.fasta \
    --out path/to/output/dir \
    --thread 32
```

## Virus Detection Module
This module detects viral load and will output a list of infecting viral species or infected cell barcodes, depending on the input and the viral index used. (*Note: For path/to/output/dir parameter, please do not include an end '/'.*) 

For bulk single-end sequencing:
```
python3 module-detection.py \
    --read read.fastq.gz \
    --virusThreshold 5 \
    --virusGenome path/to/virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32
```

For bulk paired-end sequencing (please separate paired reads by white space):
```
python3 module-detection.py \
    --read read_1.fastq.gz read_2.fastq.gz \
    --virusThreshold 5 \
    --virusGenome path/to/virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32
```

For single-cell sequencing. Please put cDNA read as first arg, CB+UMI read as second arg. 
Also, singleCellBarcode and singleUniqueMolIdent both specifiy a start position (int) and a length (int):
```
python3 module-detection.py \
    --read read_1-cDNA.fastq.gz read_2-CB+UMI.fastq.gz \
    --virusThreshold 5 \
    --virusGenome path/to/virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist whitelist.txt
```

## Integration Site Discovery Module
This is how to use Venus's integration site discovery module.
