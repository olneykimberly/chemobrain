# chemobrain
The transcriptional alternations within the brain associated with chemo therapy 


#---
01_processing
processing raw fastq files
```
cd 01_processing/
sh 01_get_readgroup_info.sh
python 02_create_config.py
```

# snakemake 
```
sbatch 03_run_snakemake.sh
```