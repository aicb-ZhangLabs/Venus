# Creating Index
This step creates the necessary mapping indices `*.genomeDir/` directories in Venus. Of note, since some viruses may lack a gtf annotation file, use of a virus gtf file is optional. (*Note: For `out` parameter, please do not include an end '/' in path.*) 

## Detection
One would have to create separate indices for the two modules in Venus, Detection and Integration. Furthermore, there are two options of indices for the Detection module, Mega-virus and Single virus. User should choose the Single virus mode if the infected viral species is already known. If not, the Mega-virus mode is recommended.

### Mega-virus index mode
For below test:
- `humanFASTA` and `humanGTF` are the human reference sequence & annotation files, respectively
-  `virusFASTA` ~ the directions to create a "mega-virus.fasta" is in the Setup section of [[README]](../../README.md)
- `module` indicates which module are we generating the index directories for
- `thread` allows for parallelization
- **(output)** `out` is the directory path for extra output files
- **(output)** `hGenome` is the index directory for either the human or hybrid genome
- **(output)** `virusGenome` is the index directory for virus

To create indices for mega-virus mode (without a virus gtf):
```   
python3 module-index.py \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA mega-virus.fasta \
    --module detection \
    --thread 32 \
    --out path/to/output/dir \
    --hGenome path/to/human.genomeDir \
    --virusGenome out_path/to/mega_virus.genomeDir
```

### Single virus index mode
For below test:
- `humanFASTA` and `humanGTF` are the human reference sequence & annotation files, respectively
- `virusFASTA` and `virusGTF` are the virus reference sequence & annotation files, respectively, in this case HIV
- `module` indicates which module are we generating the index directories for
- `thread` allows for parallelization
- **(output)** `out` is the directory path for extra output files
- **(output)** `hGenome` is the index directory for either the human or hybrid genome
- **(output)** `virusGenome` is the index directory for virus

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA NC_001802.fna \
    --virusGTF NC_001802.gtf \
    --module detection \
    --thread 32 \
    --out path/to/output/dir \
    --hGenome path/to/human.genomeDir \
    --virusGenome out_path/to/HIV.genomeDir
```


## Integration
Here, we would create a hybrid mapping index by combining human and virus reference files. Please ensure that both human and viral annotation files are in the same format (i.e. both gff or both gtf).

For below test:
- `humanFASTA` and `humanGTF` are the human reference sequence & annotation files, respectively
- `virusFASTA` and `virusGTF` are the virus reference sequence & annotation files, respectively, in this case HIV
- `module` indicates which module are we generating the index directories for
- `thread` allows for parallelization
- **(output)** `out` is the directory path for extra output files
- **(output)** `hGenome` is the index directory for either the human or hybrid genome
- **(output)** `virusGenome` is the index directory for virus

To create indices for single-virus mode (with a virus gtf):
```   
python3 module-index.py \
    --humanFASTA GCF_000001405.39_GRCh38.p13_genomic.fna \
    --humanGTF GCF_000001405.39_GRCh38.p13_genomic.gtf \
    --virusFASTA NC_001802.fna \
    --virusGTF NC_001802.gtf \
    --module integration \
    --thread 32 \
    --out path/to/output/dir \
    --hGenome path/to/hybrid.genomeDir \
    --virusGenome out_path/to/HIV.genomeDir
```
