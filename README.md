# mtdna-dlp
scripts to generate figures/tables for mtDNA DLP+ paper

## Installation (in the command line)

Clone this repository (`git clone https://github.com/reznik-lab/mtdna-dlp.git`)

Install miniconda: https://docs.conda.io/en/latest/miniconda.html

Create a new conda environment and activate it __(make sure it is activated before installing the dependencies)__:
```
conda create -n [environment]
conda activate [environment]
```

Install the following:
- VEP (https://anaconda.org/bioconda/ensembl-vep, version 95.3 recommended)
- vt (https://anaconda.org/bioconda/vt)
- pybedtools (https://anaconda.org/bioconda/pybedtools)
- bcftools (https://anaconda.org/bioconda/bcftools)
- samtools (https://anaconda.org/bioconda/samtools)
- vcf2maf (https://anaconda.org/bioconda/vcf2maf)
- pysam (https://anaconda.org/bioconda/pysam)
- gatk4 (https://anaconda.org/bioconda/gatk4)
- biopython (https://anaconda.org/conda-forge/biopython)
- numpy (https://anaconda.org/anaconda/numpy)
- pandas (https://anaconda.org/anaconda/pandas)
- matplotlib (https://anaconda.org/conda-forge/matplotlib)

```
conda install -c bioconda ensembl-vep=95.3
conda install -c bioconda vt
conda install -c bioconda pybedtools
conda install -c bioconda bcftools
conda install -c bioconda samtools
conda install -c bioconda vcf2maf
conda install -c bioconda pysam
conda install -c bioconda gatk4
conda install -c conda-forge biopython
conda install numpy
conda install pandas
conda install -c conda-forge matplotlib
```

Additionally, a VEP offline cache needs to be installed __(NOTE: CACHE MUST BE SAME VERSION AS VEP VERSION)__. Please refer to https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html for instructions on how to install a VEP cache. Due to the size of the caches, it will likely take several hours to install.

## Running the single cell pipeline (in the command line)

Navigate to the directory with the `scMTpipeline.py` file and run the following (replace all brackets):
```
conda activate [environment]
python3 scMTpipeline.py -d [data_directory] -r [reference_fasta] -w [working_directory] -l [library_id] -re [results_directory] -vc [optional_vep_cache_directory] -q [optional_mapping_quality] -Q [optional_base_quality] -s [optional_strand] -p [optional_patternlist] -t [optional_threshold]
```

Parameter descriptions:

- __Data directory__: path to directory with input .bam files
- __Reference fasta__: path to fasta file (recommended to use `[working_directory]/reference/b37/b37_MT.fa` for GRCh37 or `[working_directory]/reference/GRCh38/genome_MT.fa` for GRCh38)
- __Working directory__: path to directory with `scMTpipeline.py` file in it
- __Library ID__: name of .bam file to use as input
- __Results directory__: path to directory where results will be stored
- __(OPTIONAL )VEP cache directory__: path to directory with VEP cache
- __(OPTIONAL) Mapping quality__: minimum mapping quality (default=20)
- __(OPTIONAL) Base qualtiy__: minimum base quality (default=20)
- __(OPTIONAL) Strand__: minimum number of reads mapping to forward and reverse strand to call mutation (default=2)
- __(OPTIONAL) Patternlist__: file containing a list of filenames to process at a time (to be used when there are many files to process)
- __(OPTIONAL) Threshold__: critical threshold for calling a cell wild-type (default=0.1)

For example, a call to run the single cell pipeline with the minimum paramaters could look like this:
```
python3 scMTpipeline.py -d /my_data/ -r /my_home/mtdna-dlp/python/reference/b37/b37_MT.fa -w /my_home/mtdna-dlp/python/ -l my_file -re /my_home/mtdna-dlp/results/
```

To run `scMTpipeline.py` on the provided example data:
```
python3 scMTpipeline.py -d /mtdna-dlp/python/example_data/ -r /mtdna-dlp/python/reference/b37/b37_MT.fa -w /mtdna-dlp/python/ -l SA1101-A96155C -re /mtdna-dlp/python/results
```