# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/fig%202.png)

## Setup
This is how to set up Venus.

## Creating Index
This is how to create indices in Venus.

## Virus Detection Module
This module detects viral load and will output a list of infecting viral species or infected cell barcodes, depending on the input and the viral index used. (*Note: For out/file/name/prefix parameter, please do not include an end '/'.*) 

For bulk single-end sequencing:
```
python3 module-detection.py \
    --read SRR6944349.1_1.fastq.gz \
    --virusGenome virus.genomeDir \
    --humanGenome human.genomeDir \
    --out out/file/name/prefix \
    --readFilesCommand zcat \
    --thread 32
```

## Integration Site Discovery Module
This is how to use Venus's integration site discovery module.
