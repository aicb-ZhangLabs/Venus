# Integration Module
This module detects viral integration sites with a hybrid index. The guideFASTA file is in place to classify fusion transcripts into 3 classes of integration sites. For single-cell sequencing, simply treat the cDNA read as a single-end sequencing experiment. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Single-end sequencing
For below test:
- `read` bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- `virusGenome` and `hybridGenome` are directory paths created in the Creating Index module [[Creating Index]](../../src/module-index/module-index.md)
- `guideFASTA` file is used to classify integration sites. This file should be edited for different viruses [[guideFASTA]](../../reference_files/integrSeq.fna)
- `readFilesCommand` are commands necessary for gzipped reads
- `virusChr` specifies the assembly name for the given viral species
- `thread` allows for parallelization
- `geneBed` is a file used to convert genomic coordinates to gene names [[geneBed]](../../reference_files/genes.bed)
- **(output)** `out` is the directory path for output files

To map single-end sequencing:
```
python3 module-integration.py \
    --read bulk_1.fastq.gz \
    --virusGenome path/to/HIV.genomeDir \
    --hybridGenome path/to/hybrid.genomeDir \
    --guideFASTA integrSeq.fna \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed genes.bed \
    --out path/to/output/dir
```

## Paired-end sequencing
For below test:
- `read` (HIV) bulk_1.fastq.gz [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz) and bulk_2.fastq.gz [[bulk-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_2.fastq.gz) files are in the repo's "test_data" *(separate paired reads by white space)*
- `virusGenome` and `hybridGenome` are directory paths created in the Creating Index module [[Creating Index]](../../src/module-index/module-index.md)
- `guideFASTA` file is used to classify integration sites. This file should be edited for different viruses [[guideFASTA]](../../reference_files/integrSeq.fna)
- `readFilesCommand` are commands necessary for gzipped reads
- `virusChr` specifies the assembly name for the given viral species
- `thread` allows for parallelization
- `geneBed` is a file used to convert genomic coordinates to gene names [[geneBed]](../../reference_files/genes.bed)
- **(output)** `out` is the directory path for output files

To map paired-end sequencing:
```
python3 module-integration.py \
    --read bulk_1.fastq.gz bulk_2.fastq.gz \
    --virusGenome path/to/HIV.genomeDir \
    --hybridGenome path/to/hybrid.genomeDir \
    --guideFASTA integrSeq.fna \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed genes.bed \
    --out path/to/output/dir
```
