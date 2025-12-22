# SNPSnip

In this tutorial, we will filter a tiny SNP set from the 1001 Genomes project
(*A. thaliana*). SNPSnip is designed for the interactive filtering and
subsetting of massive raw SNP call datasets. It supports two modes of
operation, an "online" mode, which is somewhat simpler to use, however is only
usable on machines with a screen (e.g. your laptop), and therefore isn't
particularly useful for large datasets. To get around this, there is also an
"offline" mode, which requires slightly more manual work, however requires only
that you can copy files from/to a machine with a web browser, and therefore is
great for use on remote HPC servers, or even as batch jobs on a compute
cluster.

A note on filtering philosophy: each downstream analysis will tolerate
a different balance of SNP count, missing data, and sample population
structure. In general, I recommend creating a main dataset that has been
filtered relatively lightly, and where only very clearly outlying or failed
samples are removed. Then, if some analysis that requires e.g.
a highly-complete SNP matrix, further filtering and likely imputation can be
employed before that analysis specifically.


## SNPSnip pipeline stages

For both online and offline modes, SNPsnip performs the following steps of
a filtering pipeline. 

1. **Sample processing**: SNPSnip extracts a random subset of SNPs passing
   basic filters (e.g. `--maf 0.01`), and calculates per-sample stats and
   a sample PCA. You will be then presented with a web UI where you can set
   your sample filtering thresholds to exclude poor samples and (optionally)
   create subsets of samples. 
3. **Variant stats within subsets**: For each sample group (or for all samples,
   if there are no sample groups), SNPSnip will then calculate variant-level
   statistics from the random subset of SNPs from step 1. You can then set
   your thresholds per-group to exclude poor SNPs.
4. **Final Filtering**: The tool applies your sample and VCF filters to the
   full VCF file to generate filtered outputs for each group of samples.


You can specify absolute minimum thresholds on MAF, missingness, and variant
quality to subset the variants that are ever considered, which is useful if you have
a very large number of singleton or poor quality variants (use `--maf`,
`--max-missing`, `--min-qual` for this). Note that only variants passing these
thresholds are used for statistic calculation, so these parameters will
truncate the variant statistic distributions. Also note that these thresholds
are also applied to the output, so be careful not to be too harsh. 

Optionally, some predefined groups (e.g. populations, species, etc) can be
provided with the `--groups-file`, `--group-column` and `--sample-column`
arguments. These predefined groups can be used to set the default sample
grouping, which can then be refined based on the sample PCA or simply used
verbatim to define subsets of samples.


# Tutorial


## Obtain data

First we obtain a filtered subset of the 1001g data, which is only about 10MB.
Please feel free to substitute this for any standards compliant VCF, including
your own data. SNPsnip should work with any VCF or BCF files, however they
*must be indexed*! (just use bcftools index if they aren't already)

```
wget -O ath_filt-MAC5-MISS20.vcf.gz 'https://zenodo.org/records/16880445/files/ath_filt-MAC5-MISS20.vcf.gz?download=1'
wget -O ath_filt-MAC5-MISS20.vcf.gz.csi 'https://zenodo.org/records/16880445/files/ath_filt-MAC5-MISS20.vcf.gz.csi?download=1'
```

### Install

Here, I will install snpsnip in a virtual environment. For production purposes,
I would recommend using pipx to install a global binary named snpsnp, but have
the python innards in a separate virtual env to avoid version hell as SNPSnip
has quite a few dependencies.


```
python3 -m venv sspenv
source sspenv/bin/activate
pip install -e .
```


### Online mode


```
snpsnip --vcf ath_filt-MAC5-MISS20_ctg.vcf.gz --output-dir ath_online --region-size 100000000 --maf 0.001
```


### Offline mode

First, make a subset and calcuate PCA & Sample stats

```
snpsnip --vcf input.vcf.gz --offline
```

This generates a static HTML file you can download & play with to set your
thresholds. Then, you save a .json file to your PC and then copy that file
back to wherever you're running SNPsnip, then:

```
snpsnip --vcf input.vcf.gz --offline --next snpsnip_sample_filters.json
```

This makes the subsets, and calculates the SNP stats for each group of
samples you selected. This again generates a static HTML file you can use to
interactively make your SNP filtering threshold selections. Again save the
output, copy it back to where you're running SNPsnip, then:

```
snpsnip --vcf input.vcf.gz --offline --next snpsnip_variant_filters.json
```

This will generate the final files.


# All arguments

```
usage: snpsnip [-h] --vcf VCF [--output-dir OUTPUT_DIR]
               [--state-file STATE_FILE] [--temp-dir TEMP_DIR] [--host HOST]
               [--port PORT] [--offline] [--next NEXT] [--maf MAF]
               [--max-missing MAX_MISSING] [--min-qual MIN_QUAL]
               [--subset-freq SUBSET_FREQ] [--groups-file GROUPS_FILE]
               [--group-column GROUP_COLUMN] [--sample-column SAMPLE_COLUMN]
               [--processes PROCESSES] [--region-size REGION_SIZE]

SNPSnip - An interactive VCF filtering tool

options:
  -h, --help            show this help message and exit
  --vcf VCF             Input VCF/BCF file (must be indexed)
  --output-dir OUTPUT_DIR
                        Output directory for filtered VCFs
  --state-file STATE_FILE
                        State file for checkpointing
  --temp-dir TEMP_DIR   Directory for temporary files
  --host HOST           Host to bind the web server to
  --port PORT           Port for the web server
  --offline             Run in offline mode (generate static HTML)
  --next NEXT           JSON file with selections from previous step (for
                        offline mode)
  --maf MAF             Minimum minor allele frequency
  --max-missing MAX_MISSING
                        Maximum missingness rate
  --min-qual MIN_QUAL   Minimum variant quality
  --subset-freq SUBSET_FREQ
                        Fraction of SNPs to sample for interactive analysis
  --groups-file GROUPS_FILE
                        CSV or TSV file with sample and group columns for
                        predefined groups
  --group-column GROUP_COLUMN
                        Column in CSV or TSV file for predefined groups
  --sample-column SAMPLE_COLUMN
                        CSV or TSV file for sample
  --processes PROCESSES
                        Number of parallel processes to use
  --region-size REGION_SIZE
                        Size of each parallel region
```
