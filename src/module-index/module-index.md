# Creating Index
This step creates the necessary mapping indices `*.genomeDir/` directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Detection
One would have to create separate indices for the two modules in Venus, Detection and Integration. Furthermore, there are two options of indices for the Detection module, Mega-virus and Single virus. User should choose the Single virus mode if the infected viral species is already known. If not, the Mega-virus mode is recommended.

### Mega-virus index mode
For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz)
- humanGTF file download link [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)
- virusFASTA ~ the directions to create a "mega-virus.fasta" is in the Setup section of [[README]](../../README.md) *(note: there will not be an associated gtf file)*

To create indices for mega-virus mode (without a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/human.genomeDir \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusGenome out_path/to/mega_virus.genomeDir \
    --virusFASTA mega-virus.fasta \
    --module detection \
    --out path/to/output/dir \
    --thread 32
```

### Single virus index mode
For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz)
- humanGTF file download link [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)
- virusFASTA for HIV download link [[virusFASTA]](../../reference_files/NC_001802.fna)
- virusGTF for HIV download link [[virusGTF]](../../reference_files/NC_001802.gtf)

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/human.genomeDir \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusGenome out_path/to/HIV.genomeDir \
    --virusFASTA NC_001802.fna \
    --virusGTF NC_001802.gtf \
    --module detection \
    --out path/to/output/dir \
    --thread 32
```


## Integration
Here, we would create a hybrid mapping index by combining human and virus reference files. Please ensure that both human and viral annotation files are in the same format (i.e. both gff or both gtf).

For below test:
- humanFASTA file download link [[humanFASTA]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz)
- humanGTF file download link [[humanGTF]](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz)
- virusFASTA for HIV download link [[virusFASTA]](../../reference_files/NC_001802.fna)
- virusGTF for HIV download link [[virusGTF]](../../reference_files/NC_001802.gtf)

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --hGenome path/to/hybrid.genomeDir \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusGenome out_path/to/HIV.genomeDir \
    --virusFASTA NC_001802.fna \
    --virusGTF NC_001802.gtf \
    --module integration \
    --out path/to/output/dir \
    --thread 32
```
