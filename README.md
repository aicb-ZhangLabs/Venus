# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/overview.png)

## Setup
Clone the repository and set its path. 

*Please set `repo_dir` to your **own downloaded location** for easy copy and pasting when running the tutorial*
```
git clone https://github.com/aicb-ZhangLabs/Venus.git
repo_dir=/srv/disk00/cheyul1/Venus/outputs/22-05-10/Venus
```

Set up the conda envrionment. Note the spec-file.txt is for Linux platform.
```
conda create --name venus --file ${repo_dir}/reference_files/venus_spec-file.txt
conda activate venus
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

## Running with test data
There are two main modules in Venus: Detection and Integration. However, one would need to create mapping indices first prior to running these two modules. The detailed steps to do so is stated below:
- Creating Index [[link]](src/module-index/module-index.md)
- Detection Module [[link]](src/module-detection/module-detection.md)
- Integration Module [[link]](src/module-integration/module-integration.md)

The expected outputs are given in the 'test_output' directory of the repository.
