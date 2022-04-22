# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

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

## Creating Index
This step creates the necessary human and viral indices directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

### Mega-virus index mode
For below tutorial test:
- humanFASTA latest version file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
- humanGTF latest version file download link *(file is in gff3 format, which Venus has already accomodated)* [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz)
- virusFASTA ~ the directions to create a "mega-virus.fasta" is in the Setup section above *(note: there will not be an associated gtf file)*

To create indices for mega-virus mode (without a virus gtf):
```   
python3 module-index.py \
    --humanGenome path/to/human.genomeDir \
    --humanFASTA GRCh38_latest_genomic.fna \
    --humanGTF GRCh38_latest_genomic.gff \
    --virusGenome path/to/mega_virus.genomeDir \
    --virusFASTA mega-virus.fasta \
    --out path/to/output/dir \
    --thread 32
```

### Single virus index mode
For below tutorial test:
- virusFASTA for HIV download link [[virusFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.fna.gz)
- virusGTF for HIV download link [[virusGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.gtf.gz)

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --humanGenome path/to/human.genomeDir \
    --humanFASTA GRCh38_latest_genomic.fna \
    --humanGTF GRCh38_latest_genomic.gff \
    --virusGenome path/to/HIV.genomeDir \
    --virusFASTA GCF_000864765.1_ViralProj15476_genomic.fna \
    --virusGTF GCF_000864765.1_ViralProj15476_genomic.gtf \
    --out path/to/output/dir \
    --thread 32
```

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

## Integration Site Discovery Module
This is how to use Venus's integration site discovery module.
