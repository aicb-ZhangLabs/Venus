# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Virus Detection Module
This module detects viral load and will output a list of infecting viral species or infected cell barcodes, depending on the input and the viral index used. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

### Bulk, single-end sequencing
For below tutorial test:
- virusChrRef file is in the repo's "reference_files/virus_chr-ref.tsv"
- GenomeDir's are directories created in the Creating Index section above
- bulk_1.fastq.gz (HIV) file is in the repo's "test_files/bulk_1.fastq.gz"

To map bulk single-end sequencing:
```
python3 module-detection.py \
    --read bulk_1.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef reference_files/virus_chr-ref.tsv \
    --virusGenome path/to/mega_virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32
```

### Bulk, paired-end sequencing
For below tutorial test:
- virusChrRef file is in the repo's "reference_files/virus_chr-ref.tsv"
- GenomeDir's are directories created in the Creating Index section above
- bulk_1.fastq.gz (HIV) file is in the repo's "test_files/bulk_1.fastq.gz"
- bulk_2.fastq.gz (HIV) file is in the repo's "test_files/bulk_2.fastq.gz"

To map bulk paired-end sequencing (please separate paired reads by white space):
```
python3 module-detection.py \
    --read bulk_1.fastq.gz bulk_2.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef reference_files/virus_chr-ref.tsv \
    --virusGenome path/to/mega_virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32
```

### Single-cell sequencing
For below tutorial test:
- virusChrRef file is in the repo's "reference_files/virus_chr-ref.tsv"
- GenomeDir's are directories created in the Creating Index section above
- singlecell_1cDNA.fastq.gz (HIV) file is in the repo's "test_files/singlecell_1cDNA.fastq.gz"
- singlecell_2CB+UMI.fastq.gz (HIV) file is in the repo's "test_files/singlecell_2CB+UMI.fastq.gz"
- singeWhitelist file download link [[singleWhitelist]](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz)

To map single-cell sequencing. Please put cDNA read as first arg, CB+UMI read as second arg. 
Also, `singleCellBarcode` and `singleUniqueMolIdent` both specifiy a start position (int) and a length (int):
```
python3 module-detection.py \
    --read singlecell_1cDNA.fastq.gz singlecell_2CB+UMI.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef reference_files/virus_chr-ref.tsv \
    --virusGenome path/to/mega_virus.genomeDir \
    --humanGenome path/to/human.genomeDir \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --thread 32 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist 3M-february-2018.txt
```
