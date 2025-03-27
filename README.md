# SNPSnip - An interactive VCF Filtering Tool

SNPSnip is a command-line tool with an interactive web interface for filtering VCF files in multiple stages.

### Prerequisites

- Python 3.8 or higher
- bcftools must be installed and available in your PATH

### Install from source

```bash
pip install git+https://github.com/gekkonid/snpsnip.git
```

## Usage

### Online mode

```bash
snpsnip --output-dir filtered_results --vcf input.vcf.gz --maf 0.05 --max-missing 0.1 --min-qual 30
```

This will conduct three phases of analysis, with a web UI to select threholds etc between these steps.

1. **Initial Processing**: SNPSnip extracts a random subset of SNPs passing
   basic filters, and calculates per-sample stats and a sample PCA. You will be
   then presented with a web UI where you can set your sample filtering
   thresholds to exclude poor samples and (optionally) create subsets of
   samples.
3. **Variant Filtering**: For each sample group, SNPSnip will then calculate
   variant-level statistics from the random subset of SNPs (from step 1). You
   can then set your thresholds per-group to exclude poor SNPs.
4. **Final Filtering**: The tool applies your sample and VCF filters to the
   full VCF file to generate filtered outputs for each group of samples.

Optionally, some predefined groups (e.g. populations, species, etc) can be
provided with the `--groups-file`, `--group-column` and `--sample-column`
arguments. These predefined groups can then be refined based on the sample PCA
or simply used verbatim to define subsets of samples.

### Offline mode

You can also run this in an "offline" mode, useful for example on clusters

```bash
# First, make a subset and calcuate PCA & Sample stats
snpsnip --vcf input.vcf.gz --offline

# This generates a static HTML file you can download & play with to set your
# thresholds. This saves a .json file to your PC. Copy that file back to
# wherever you're running SNPsnip, then:

snpsnip --vcf input.vcf.gz --offline --next snpsnip_sample_filters.json

# This makes the subsets, and calculates the SNP stats for each group of
# samples you selected. This again generates a static HTML file you can use to
# interactively make your SNP filtering threshold selections. Again save the
# output, copy it back to where you're running SNPsnip, then:

snpsnip --vcf input.vcf.gz --offline --next snpsnip_variant_filters.json

# This will generate the final files.
```

## Workflow


## License

MPL2
