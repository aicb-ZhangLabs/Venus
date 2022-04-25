# Integration Module
This module detects viral integration sites with a hybrid index. The guideFASTA file is in place to classify fusion transcripts into 3 classes of integration sites. For single-cell sequencing, simply treat the cDNA read as a single-end sequencing experiment. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Single-end sequencing
For below test:
- read bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- .genomeDir are directories created in the Creating Index module
- guideFASTA file is in the repo's "reference_files" [[guideFASTA]](https://github.com/aicb-ZhangLabs/Venus/blob/721b376f603c7917e88981e20b7098dd80a4aedf/reference_files/integrSeq.fna)
- virusChr parameter specifies the assembly name given for a specific viral species
- geneBed file, used to convert genomic coordinates to gene names, is in the repo's "reference_files" [[geneBed]](https://github.com/aicb-ZhangLabs/Venus/raw/main/reference_files/genes.bed)

To map single-end sequencing:
```
python3 module-detection.py \
    --read bulk_1.fastq.gz \
    --virusGenome path/to/HIV.genomeDir \
    --hybridGenome path/to/hybrid.genomeDir \
    --guideFASTA integrSeq.fna \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed genes.bed
```

## Paired-end sequencing
For below test:
- read bulk_1.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-1]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_1.fastq.gz)
- read bulk_2.fastq.gz (HIV) file is in the repo's "test_data" [[bulk-2]](https://github.com/aicb-ZhangLabs/Venus/raw/main/test_data/bulk_2.fastq.gz)
- .genomeDir are directories created in the Creating Index module
- guideFASTA file is in the repo's "reference_files" [[guideFASTA]](https://github.com/aicb-ZhangLabs/Venus/blob/721b376f603c7917e88981e20b7098dd80a4aedf/reference_files/integrSeq.fna)
- virusChr parameter specifies the assembly name given for a specific viral species
- geneBed file, used to convert genomic coordinates to gene names, is in the repo's "reference_files" [[geneBed]](https://github.com/aicb-ZhangLabs/Venus/raw/main/reference_files/genes.bed)

To map paired-end sequencing (please separate paired reads by white space):
```
python3 module-detection.py \
    --read bulk_1.fastq.gz bulk_2.fastq.gz \
    --virusGenome path/to/HIV.genomeDir \
    --hybridGenome path/to/hybrid.genomeDir \
    --guideFASTA integrSeq.fna \
    --out path/to/output/dir \
    --readFilesCommand zcat \
    --virusChr NC_001802.1 \
    --thread 32 \
    --geneBed genes.bed
```
