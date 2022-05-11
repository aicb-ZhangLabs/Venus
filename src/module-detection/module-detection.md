# Detection Module
This module detects viral load and will output a list of infecting viral species (with infected cell barcodes if single-cell). (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Bulk sequencing
### Single-end reads
For below test:
- `read` bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- `virusThreshold` is the number of minimum viral transcripts to count viral species infection
- `virusChrRef` file is for mapping accession id to viral species name [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- `virusGenome` and `humanGenome` are directory paths created in the Creating Index module [[Creating Index]](../../src/module-index/module-index.md)
- `readFilesCommand` are commands necessary for gzipped reads
- `thread` allows for parallelization
- **(output)** `out` is the directory path for output files

To map bulk single-end sequencing:
```
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --out ${out_dir}/detection/bulk_single-end
```

### Paired-end reads
For below test:
- `read` (HIV) bulk_1.fastq.gz [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz) and bulk_2.fastq.gz [[bulk-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_2.fastq.gz) files are in the repo's "test_data" *(please separate by white space)*
- `virusThreshold` is the number of minimum viral transcripts to count viral species infection
- `virusChrRef` file is for mapping accession id to viral species name [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- `virusGenome` and `humanGenome` are directory paths created in the Creating Index module [[Creating Index]](../../src/module-index/module-index.md)
- `readFilesCommand` are commands necessary for gzipped reads
- `thread` allows for parallelization
- **(output)** `out` is the directory path for output files

To map bulk paired-end sequencing:
```
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/bulk_1.fastq.gz ${repo_dir}/test_data/bulk_2.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --out ${out_dir}/detection/bulk_paired-end
```

## Single-cell sequencing
For below test:
- `read` (HIV) singlecell_1cDNA.fastq.gz [[singlecell-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/singlecell_1cDNA.fastq.gz) and singlecell_2CB+UMI.fastq.gz [[singlecell-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/singlecell_2CB%2BUMI.fastq.gz) files are in the repo's "test_data" *(Please put cDNA read as 1st arg, CB+UMI read as 2nd arg)*
- `virusThreshold` is the number of minimum viral transcripts to count viral species infection
- `virusChrRef` file is for mapping accession id to viral species name [[virusChrRef]](../../reference_files/virus_chr-ref.tsv)
- `virusGenome` and `humanGenome` are directory paths created in the Creating Index module [[Creating Index]](../../src/module-index/module-index.md)
- `readFilesCommand` are commands necessary for gzipped reads
- `thread` allows for parallelization
- `singleCellBarcode` and `singleUniqueMolIdent` both specifiy a start position (int) and a length (int) for CB and UMI, respecitvely
- `singeWhitelist` is the barcode whitelist file [[singleWhitelist]](https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz)
- **(output)** `out` is the directory path for output files

To map single-cell sequencing:
```
python3 ${repo_dir}/src/module-detection/module-detection.py \
    --read ${repo_dir}/test_data/singlecell_1cDNA.fastq.gz ${repo_dir}/test_data/singlecell_2CB+UMI.fastq.gz \
    --virusThreshold 5 \
    --virusChrRef ${repo_dir}/reference_files/virus_chr-ref.tsv \
    --virusGenome ${out_dir}/indices/mega_virus.genomeDir \
    --humanGenome ${out_dir}/indices/human.genomeDir \
    --readFilesCommand zcat \
    --thread 32 \
    --singleCellBarcode 1 16 \
    --singleUniqueMolIdent 17 10 \
    --singleWhitelist ${repo_dir}/reference_files/3M-february-2018.txt \
    --out ${out_dir}/detection/single-cell
```
