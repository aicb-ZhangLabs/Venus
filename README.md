# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Overview
![alt text](https://github.com/aicb-ZhangLabs/Venus/blob/main/overview.png)

## Setup
Clone the repository.
```
git clone https://github.com/aicb-ZhangLabs/Venus.git
```

Set up the conda envrionment. Note the spec-file.txt is for Linux platform.
```
conda create --name venus --file reference_files/venus_spec-file.txt
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
There are two main modules in Venus: Detection and Integration. However, one would need to create mapping indices first prior to running these two modules. The detailed steps to do so is detailed below:
- Creating Index [[link]](https://github.com/aicb-ZhangLabs/Venus/blob/87936f8f52769f0dd839925d94ddc229b052d5cf/src/module-index/module-index.md)
- Detection Module [[link]](https://github.com/aicb-ZhangLabs/Venus/blob/87936f8f52769f0dd839925d94ddc229b052d5cf/src/module-detection/module-detection.md)
- Integration Module [[link]](https://github.com/aicb-ZhangLabs/Venus/blob/87936f8f52769f0dd839925d94ddc229b052d5cf/src/module-integration/module-integration.md)
