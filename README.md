# Mouse Mammary Gland Cell-Sorted RNA-seq

This repository contains the data-processing pipeline scripts and analysis code for mammary gland cell-sorted total RNA-seq in hybrid mice from reciprocal cross between CAST and B6 mouse strains. 

The main data processing pipeline is built using _snakemake_, and was run the University of Cambridge HPC using the following command: 

```bash
snakemake --use-conda --cluster-status workflow/snakemake_slurm_status.py --cluster "sbatch --parsable -A <ACCOUNT> --mail-type=FAIL -J {rulename} -o logs/slurm/{rulename}-job#%j.log -p icelake,icelake-himem,cclake,skylake,skylake-himem -c {threads} -t {params.runtime}" --jobs 10
```

The _snakemake_ scripts and config file can be found in the `workflow` folder.
The pipeline expects raw read data to be located in `data/raw/reads`, with sub-folders named according to the sequencing run. See the file `read_info.csv` for details about the input file locations. 
The read data can be downloaded from [add link to SRA repository upon publication].

Downstream data analysis and visualisation was done in _R_, and the code for this is available in the `analysis` folder (these are not part of the _snakemake_ workflow). 
