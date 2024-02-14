# mtdna-dlp

## Installation (in the command line)

Clone this repository:
```
git clone https://github.com/reznik-lab/mtdna-dlp.git
```

Install miniconda (at least Python v3.6): https://docs.conda.io/en/latest/miniconda.html

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

Additionally, a VEP offline cache needs to be installed __(NOTE: CACHE MUST BE SAME VERSION AS VEP VERSION)__. Please refer to http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html for instructions on how to install a VEP cache. Due to the size of the caches, it will likely take several hours to install.

You will also need to ensure that the reference genome fasta has a dictionary. Afer downloading the correct reference genome (make sure the mitochondrial chromosome name matches the one in your bam file, i.e. chrM or MT). 

To download the reference genome, create a directory [reference] to store it and run the following:

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar xzf refdata-gex-GRCh38-2020-A.tar.gz 
```
To make "genome.dict" file run the following in the same directory as your reference genome.fa file:

```
java -jar [your_directory]/mtdna-dlp/python/reference/picard.jar CreateSequenceDictionary R=genome.fa O=genome.dict
```


## Bulk pipeline

### Running the bulk pipeline (in the command line)

Subsetting for mitochondrial regions (where chrM can either be chrM or MT, depending on how it is annotated in the possorted_bam file):

```
conda activate [environment]
samtools view -b [library_id].bam chrM > [library_id]_chrM.bam
```
Index the input .bam file:
```
samtools index [library_id]_chrM.bam
```
*The index file needs to always be newer than the bam file, otherwise the pipeline will throw an error

Navigate to the directory with the `bulkpipeline.py` file and run the following (replace all brackets):
```
python3 bulkpipeline.py -d [data_directory] -w [working_directory] -l [library_id] -re [results_directory] -vc [optional_vep_cache_directory] -q [optional_mapping_quality] -Q [optional_base_quality] -s [optional_strand] -t [optional_threshold] -g [optional_genome] -r [optional_reference_fasta] -n [optional_normal] -nd [optional_normaldir] -m [optional_molecule] -c [optional_mincounts]
```

Parameter descriptions:

- __Tumor file__: Path of the tumor sample
- __Working directory__: path to directory with `bulkpipeline.py` file in it
- __Results directory__: path to directory where results will be stored
- __(OPTIONAL) VEP cache directory__: path to directory with VEP cache (default="$HOME/.vep")
- __(OPTIONAL) Mapping quality__: minimum mapping quality (default=20)
- __(OPTIONAL) Base quality__: minimum base quality (default=20)
- __(OPTIONAL) Strand__: minimum number of reads mapping to forward and reverse strand to call mutation (default=2)
- __(OPTIONAL) Threshold__: critical threshold for calling a cell wild-type (default=0.1)
- __(OPTIONAL) Genome__: genome version (supported genomes are GRCh37, GRCh38, GRCm38, and mm10, default="GRCh37")
- __(OPTIONAL) Reference fasta__: path to fasta file (by default will use a file from the reference folder that matches the genome, but best to use same reference used to align fastqs)
- __(OPTIONAL) Normal__: Path of the matched normal sample (default="")
- __(OPTIONAL) Molecule__: type of molecule (dna or rna, default="dna")
- __(OPTIONAL) Minimum counts__: minimum number of read counts for MTvariantpipeline (default=100)

For example, a call to run the bulk pipeline with the required paramaters could look like this:
```
python3 bulkpipeline.py -t /my_data/file.bam -w /my_home/mtdna-dlp/python/ -re /my_home/mtdna-dlp/results/
```

### Bulk pipeline outputs

The bulk pipeline should have the following files and subdirectories in the results directory:

- [sample].bam.maf
- [sample]_mutsig.tsv
- MuTect2_results
- MTvariant_results
- TEMPMAFfiles

## Single cell pipeline

### Splitting the bam file into one file per cell

Using the output from cellranger (i.e. /outs/filtered_feature_bc_matrix/barcodes.tsv.gz and /outs/possorted_genome_bam.bam), run the `split_bam.py` file. For example:
```
python3 split_bam.py possorted_genome_bam.bam output_directory --barcode_csv barcodes.tsv.gz
```

### Running the single cell pipeline (in the command line)
Note, check the limit of open files (ulimit -n), and adjust as needed. 

Index the split input .bam files:
```
conda activate [environment]
for each in [library_id]/*.bam; do samtools index "$each"; done
```

Navigate to the directory with the `scMTpipeline.py` file and run the following (replace all brackets):
```
python3 scMTpipeline.py -d [data_directory] -w [working_directory] -l [library_id] -re [results_directory] -vc [optional_vep_cache_directory] -q [optional_mapping_quality] -Q [optional_base_quality] -s [optional_strand] -t [optional_threshold] -g [optional_genome] -r [optional_reference_fasta] -m [optional_molecule] -c [optional_mincounts]
```

Parameter descriptions:

- __Data directory__: path to directory with split input .bam files
- __Working directory__: path to directory with `scMTpipeline.py` file in it
- __Library ID__: name of sample
- __Results directory__: path to directory where results will be stored
- __(OPTIONAL) VEP cache directory__: path to directory with VEP cache (default="$HOME/.vep")
- __(OPTIONAL) Mapping quality__: minimum mapping quality (default=20)
- __(OPTIONAL) Base quality__: minimum base quality (default=20)
- __(OPTIONAL) Strand__: minimum number of reads mapping to forward and reverse strand to call mutation (default=2)
- __(OPTIONAL) Threshold__: critical threshold for calling a cell wild-type (default=0.1)
- __(OPTIONAL) Genome__: genome version (supported genomes are GRCh37, GRCh38, GRCm38, and mm10, default="GRCh37)
- __(OPTIONAL) Reference fasta__: path to fasta file (by default will use a file from the reference folder that matches the genome, but best to use same reference used to align fastqs)
- __(OPTIONAL) Molecule__: type of molecule (dna or rna, default="dna")
- __(OPTIONAL) Minimum counts__: minimum number of read counts for MTvariantpipeline (default=100)

For example, a call to run the single cell pipeline with the required paramaters could look like this:
```
python3 scMTpipeline.py -d /my_data/ -w /my_home/mtdna-dlp/python/ -l my_sample -re /my_home/mtdna-dlp/results/
```

To run `scMTpipeline.py` on the provided example data:
```
python3 scMTpipeline.py -d /mtdna-dlp/python/example_data/ -w /mtdna-dlp/python/ -l SA1101-A96155C -re /mtdna-dlp/python/results/
```

### Single cell pipeline outputs

The single cell pipeline should have the following files and subdirectories in the results directory:

- filteredfiles
- merged
- mergedTEMPMAFfiles
- MuTect2_results
- MTvariant_results
- TEMPMAFfiles
- [sample]-merged.fillout
- [sample]_master.tsv
- [sample]_vaf.tsv
- [sample]_binary.tsv
- [sample]_variants.tsv
- [sample]_depth.tsv
- [sample]_mutprob.tsv
- [sample]_mutsig.tsv
- [sample]_heteroplasmy.pdf
- [sample]_haplogroups.txt (only for human DNA)
- filtered[sample]-merged.bam (only for human DNA)
