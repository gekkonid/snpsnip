# SNPSnip - Interactive VCF Filtering Tool

SNPSnip is a command-line tool with an interactive web interface for filtering VCF files in multiple stages.

### Prerequisites

- Python 3.7 or higher
- bcftools must be installed and available in your PATH

### Install from source

```bash
pip install git+https://github.com/gekkonid/snpsnip.git
```

## Usage

Basic usage:

```bash
snpsnip --output-dir filtered_results --vcf input.vcf.gz --maf 0.05 --max-missing 0.1 --min-qual 30
```

If you need, you can change the web server settings:

```bash
snpsnip --vcf input.vcf.gz --host 0.0.0.0 --port 8080
```

## Workflow

1. **Initial Processing**: SNPSnip extracts a random subset of SNPs passing basic filters.
2. **Sample Filtering**: The web UI allows you to filter samples based on quality metrics and PCA clustering.
3. **Variant Filtering**: For each sample group, set filtering thresholds for various metrics.
4. **Final Processing**: The tool applies your filters to the full VCF file to generate filtered outputs.

## License

MPL2
