# chemobrain
The transcriptional alternations within the brain associated with chemo therapy 

| Group                     | Count   | Sex (F, M)|
| ------------------------- |:-------:|:-------:|
| SALINE                    | 12      | (6,6)   |
| LEUC                      | 8       | (4,4)   |
| MTX                       | 8       | (5,3)   |
| MTXLEUC                   | 14      | (5,9)   |

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing of bulk RNA-sequencing data
-   Generation of manuscript figures 
-   Generation of shiny app for exploration of the results, view app [here]( https://fryerlab.shinyapps.io/chemobrain_mice/)


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the conda environment:
```
conda env create -n chemobrain --file chemobrain.yml

# To activate this environment, use
#
#     $ conda activate chemobrain
#
# To deactivate an active environment, use
#
#     $ conda deactivate chemobrain
```

## Reference genome
Reference genome and annotation were downloaded prior to running snakemake. 
Ensembl refdata-gex-GRCm39-2024-A fasta

```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir refdata-gex-GRCm39-2024-A_STARv2.7.11_150sjdb  --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 150
```

## Snakemake for trimming reads and alignment to the reference genome
See the snakefile for specific details for each step. 
The config file was generated using the 01 and 02 scripts for obtaining sample information and creating the config file. 
```
Snakemake -s Snakefile 
```
Output includes STAR reads per gene, which can then be read into R for differential expression.

## Variance assessment and differential expression
```
R 04_differential_expression.Rmd
```
Output is differentially expressed genes for all pairwise comparisons.

## References
All packages used in this workflow are publicly available. If you use this workflow please cite the packages used. 
If you use the data in this workflow cite the following:
[et al. 2026]()

## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |