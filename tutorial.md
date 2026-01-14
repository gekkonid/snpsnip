# SNPSnip

In this tutorial, we will use SNPSnip's built-in example data generator to
create a synthetic VCF file and demonstrate the filtering workflow. SNPSnip is
designed for the interactive filtering and subsetting of massive raw SNP call
datasets. It supports two modes of operation, an "online" mode, which is
somewhat simpler to use, however is only usable on machines with a screen (e.g.
your laptop), and therefore isn't particularly useful for large datasets. To
get around this, there is also an "offline" mode, which requires slightly more
manual work, however requires only that you can copy files from/to a machine
with a web browser, and therefore is great for use on remote HPC servers, or
even as batch jobs on a compute cluster.

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


## Install

Here, I will install snpsnip in a local virtual environment. For production
purposes, I would recommend using pipx to install a global binary named snpsnp,
but have the python innards in a separate virtual env to avoid version hell as
SNPSnip has quite a few dependencies.

```bash
python3 -m venv sspenv
source sspenv/bin/activate
python3 -m pip install -e .
```

## Generate example data

SNPSnip includes a built-in data generator that creates synthetic VCF files for
testing and demonstration. This is the easiest way to get started and explore
SNPSnip's features. Of course, you can substitute this for any standards
compliant VCF, including your own data. SNPsnip should work with any VCF or BCF
files, however they *must be indexed* (just use bcftools index if they aren't
already), and ideally contain contig sizes in the header. If you supply a VCF
that we can't work with, the error message should tell you how to solve the
problem (e.g. by indexing the VCF, or supplying the fasta index of the
reference sequences). If you don't find the error message sufficiently helpful,
then your file is incompatible in some novel way, please contact me or raise an
issue on github.


First, let's generate a test dataset with 50 samples and 1M SNPs, including example filter
files:

```bash
snpsnip-generate-example -o example.vcf.gz -n 50 -m 1000000 --with-filters
```

This will create:
- `example.vcf.gz` and `example.vcf.gz.csi` (VCF and index)
- Two example filter files. (You will normally generate these files by loading
  a Web UI, making your selections, and saving the file from your browser. They
  are generated here for convenience and testing.)
    - `snpsnip_sample_filters.json`
    - `snpsnip_variant_filters.json`

If you want, you can customize the generated data using additional options. See
`snpsnip-generate-example --help`.

## Offline mode

Offline mode is ideal for remote servers or HPC environments where you can't run
a web browser. It generates static HTML files that you can download and view
locally.

If you generated example data with `--with-filters`, you can skip the
interactive HTML steps, and supply the JSON files made by the example data
generator in the steps below. For a more practical walk-through, open the HTML,
play with the filters and thresholds, and export your selections as you would
for your own data, saving and using the JSON files that the web UI generates.

### Step 1: Calculate sample statistics and PCA

```bash
snpsnip --vcf example.vcf.gz --output-dir example_offline --offline --maf 0.01
```

This generates a static HTML file (`example_offline/snpsnip_sample_filters.html`)
you can copy to your local machine (if you're running this remotely on a
server), and open in your browser. Use this web interface to set sample
filtering thresholds and define sample groups. Save your selections as a JSON
file.

### Step 2: Calculate variant statistics per group

Copy the resulting JSON file back to wherever you are running SNPSnip, then:

```bash
snpsnip --vcf example.vcf.gz --output-dir example_offline --offline --next snpsnip_sample_filters.json
```

This calculates SNP statistics for each sample group and generates another HTML
file (`example_offline/snpsnip_variant_filters.html`) for setting variant
filtering thresholds. Save your selections as a JSON file.

### Step 3: Generate final filtered VCFs

Copy the JSON file back to wherever you're running SNPSnip, then:

```bash
snpsnip --vcf example.vcf.gz --output-dir example_offline --offline --next snpsnip_variant_filters.json
```

This generates the final filtered VCF files for each sample group based on your
threshold selections.


## Online mode

Alternatively, SNPSnip has an online mode that launches an interactive web
interface on your local machine:

```
snpsnip --vcf example.vcf.gz --output-dir example_online --maf 0.01
```

This will do all of the steps above. The web interface allows you to
interactively explore sample PCA, set filtering thresholds, and define sample
groups all in one session, with VCF processing happening automatically when you
click buttons in your browser. The obvious limitation here is that your web
browser, your raw VCF, and a reasonably beefy CPU must be on the same machine,
which is not typically the case with large datasets.


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
