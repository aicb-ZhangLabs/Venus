# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/overview.png)

## Parameters
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/parameters.png)

## Setup
Clone the repository. 
```
git clone https://github.com/aicb-ZhangLabs/Venus.git
```

Please set `repo_dir` to your **own downloaded location** and `out_dir` to your **own output location** for easy *copy and pasting* when running Venus tutorial.
```
repo_dir=/srv/disk00/cheyul1/Venus/outputs/22-05-10/Venus
out_dir=/srv/disk00/cheyul1/Venus/outputs/22-05-10/testing
```

Set up the conda envrionment. Note the spec-file.txt is for Linux platform.
```
conda create --name venus --file ${repo_dir}/reference_files/venus_spec-file.txt
conda activate venus
```

One can download the mega-virus.fasta file directly from here [[link]](https://drive.google.com/file/d/1R3cGDuRj8IMSVFn0O6QJFTRkT5CLpo4f/view?usp=sharing). Please then move the file into the `reference_files` directory.
```
mv mega-virus.fasta ${repo_dir}/reference_files/
```

One should also download the human (vGRCh38) fasta, gtf files from NCBI in addition to the 10X single-cell CB whitelist we will use:
```
cd ${repo_dir}/reference_files
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
gunzip 3M-february-2018.txt.gz
```

## Running with test data
There are two main modules in Venus: Detection and Integration. However, one would need to create mapping indices first prior to running these two modules. The detailed steps to do so is stated below:
- Creating Index [[link]](src/module-index/module-index.md)
- Detection Module [[link]](src/module-detection/module-detection.md)
- Integration Module [[link]](src/module-integration/module-integration.md)

The expected outputs are given in the 'test_output' directory of the repository.
