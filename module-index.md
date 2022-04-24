# Venus
**Venus:** an efficient Virus infection detection and fusion Site discovery method using single-cell and bulk RNA-seq data

## Creating Index
This step creates the necessary mapping indices directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. One would have to create separate indices for the two modules in Venus, Detection and Integration. Furthermore, there are two options of indices for the Detection module, Mega-virus and Single virus. User should choose the Single virus mode if the infected viral species is already known. If not, the Mega-virus mode is recommended. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

### Detection: Mega-virus index mode
For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
- humanGTF file download link *(file is in gff3 format, which Venus has already accomodated)* [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz)
- virusFASTA ~ the directions to create a "mega-virus.fasta" is in [[README]](https://github.com/aicb-ZhangLabs/Venus/blob/8b435737b902364b0a0b95813e6e5264b3d7753a/README.md) *(note: there will not be an associated gtf file)*

To create indices for mega-virus mode (without a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/human.genomeDir \
    --humanFASTA GRCh38_latest_genomic.fna \
    --humanGTF GRCh38_latest_genomic.gff \
    --virusGenome path/to/mega_virus.genomeDir \
    --virusFASTA mega-virus.fasta \
    --module detection \
    --out path/to/output/dir \
    --thread 32
```

### Detection: Single virus index mode
For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
- humanGTF file download link *(file is in gff3 format, which Venus has already accomodated)* [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz)
- virusFASTA for HIV download link [[virusFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.fna.gz)
- virusGTF for HIV download link [[virusGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.gtf.gz)

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/human.genomeDir \
    --humanFASTA GRCh38_latest_genomic.fna \
    --humanGTF GRCh38_latest_genomic.gff \
    --virusGenome path/to/HIV.genomeDir \
    --virusFASTA GCF_000864765.1_ViralProj15476_genomic.fna \
    --virusGTF GCF_000864765.1_ViralProj15476_genomic.gtf \
    --module detection \
    --out path/to/output/dir \
    --thread 32
```


### Integration
For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz)
- humanGTF file download link *(file is in gff3 format, which Venus has already accomodated)* [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz)
- virusFASTA for HIV download link [[virusFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.fna.gz)
- virusGTF for HIV download link [[virusGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_genomic.gtf.gz)

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/hybrid.genomeDir \
    --humanFASTA GRCh38_latest_genomic.fna \
    --humanGTF GRCh38_latest_genomic.gff \
    --virusGenome path/to/HIV.genomeDir \
    --virusFASTA GCF_000864765.1_ViralProj15476_genomic.fna \
    --virusGTF GCF_000864765.1_ViralProj15476_genomic.gtf \
    --module integration \
    --out path/to/output/dir \
    --thread 32
```
