# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Creating Index
This step creates the necessary human and viral indices directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

### (Detection) Mega-virus index mode
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

### (Detection) Single virus index mode
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
