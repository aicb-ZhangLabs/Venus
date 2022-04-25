# Detection Module
This module detects viral load and will output a list of infecting viral species (with infected cell barcodes if single-cell). (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Bulk sequencing
### Single-end reads
For below test:
- read bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- virusChrRef file is in the repo's "reference_files" [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- .genomeDir are directories created in the Creating Index module

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

### Paired-end reads
For below test:
- read bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- read bulk_2.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_2.fastq.gz)
- virusChrRef file is in the repo's "reference_files" [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- .genomeDir are directories created in the Creating Index module

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

## Single-cell sequencing
For single-cell, please ensure that the cDNA read is put first before the cell barcode + UMI read. Also, `singleCellBarcode` and `singleUniqueMolIdent` both specifiy a start position (int) and a length (int)

For below test:
- read singlecell_1cDNA.fastq.gz (HIV) file is in the repo's "test_data" [[singlecell-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/singlecell_1cDNA.fastq.gz)
- read singlecell_2CB+UMI.fastq.gz (HIV) file is in the repo's "test_data" [[singlecell-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/singlecell_2CB%2BUMI.fastq.gz)
- virusChrRef file is in the repo's "reference_files" [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- .genomeDir are directories created in the Creating Index module
- singeWhitelist file download link [[singleWhitelist]](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz)

To map single-cell sequencing. Please put cDNA read as first arg, CB+UMI read as second arg:
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
