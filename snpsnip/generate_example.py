#!/usr/bin/env python3
"""
Generate example VCF files for testing SNPSnip

Creates synthetic VCF files with:
- Configurable number of samples (N)
- Configurable number of SNPs (M)
- Uniform allele frequency distribution
- Exponentially distributed missing data
"""

import argparse
import logging
import numpy as np
import os
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_example_vcf(
    output_path,
    n_samples=20,
    n_snps=10000,
    missing_prop=0.05,
    n_chroms=2,
    compress=True,
    random_seed=None
):
    """
    Generate an example VCF file with synthetic data

    Parameters
    ----------
    output_path : str or Path
        Path to output VCF file (will be .vcf.gz if compress=True)
    n_samples : int, default=20
        Number of samples to generate
    n_snps : int, default=10000
        Total number of SNPs to generate
    missing_prop : float, default=0.05
        Proportion of missing genotypes (0.0 to 1.0).
        Missing data is exponentially distributed across samples.
    n_chroms : int, default=2
        Number of chromosomes to split SNPs across
    compress : bool, default=True
        Whether to compress and index the VCF with bcftools
    random_seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    str
        Path to the generated VCF file

    Notes
    -----
    - Allele frequencies are uniformly distributed between 1%-99%
    - Quality scores vary between 30 and 90
    - Read depth varies between 10 and 200
    - Missing data follows exponential distribution (some samples have more missing data)

    Examples
    --------
    >>> from snpsnip.generate_example import generate_example_vcf
    >>> vcf_path = generate_example_vcf(
    ...     "example.vcf.gz",
    ...     n_samples=50,
    ...     n_snps=50000,
    ...     missing_prop=0.1
    ... )
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    output_path = Path(output_path)

    # Check if bcftools is available if compression requested
    if compress:
        try:
            subprocess.run(["bcftools", "--version"],
                         capture_output=True, check=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError, subprocess.TimeoutExpired):
            logger.warning("bcftools not found - VCF will not be compressed")
            compress = False

    logger.info(f"Generating example VCF with {n_samples} samples and {n_snps} SNPs...")
    logger.info(f"Missing data proportion: {missing_prop:.2%} (exponentially distributed)")

    # Generate sample names
    # Use appropriate padding based on number of samples
    if n_samples < 100:
        samples = [f"sample_{i:02d}" for i in range(1, n_samples + 1)]
    elif n_samples < 10000:
        samples = [f"sample_{i:04d}" for i in range(1, n_samples + 1)]
    else:
        samples = [f"sample_{i:06d}" for i in range(1, n_samples + 1)]

    # Generate missing data probabilities for each sample (exponentially distributed)
    # This creates variation in missingness across samples
    if missing_prop > 0:
        # Generate exponentially distributed values
        # Scale parameter controls the mean of the distribution
        scale = 1.0
        missing_probs = np.random.exponential(scale, n_samples)
        # Normalize to get desired overall proportion
        missing_probs = missing_probs / missing_probs.sum() * missing_prop * n_samples
        # Cap at 1.0 (can't have more than 100% missing)
        missing_probs = np.minimum(missing_probs, 0.95)
    else:
        missing_probs = np.zeros(n_samples)

    logger.info(f"Missing data per sample: min={missing_probs.min():.2%}, "
                f"max={missing_probs.max():.2%}, mean={missing_probs.mean():.2%}")

    # Build VCF header
    vcf_lines = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    ]

    # Add contig headers
    chrom_length = 100_000_000  # 100 Mb per chromosome
    for i in range(1, n_chroms + 1):
        vcf_lines.append(f"##contig=<ID=chr{i},length={chrom_length}>")

    # Add INFO and FORMAT headers
    vcf_lines.extend([
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
    ])

    # Add column header
    vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))

    # Generate variants
    snps_per_chrom = n_snps // n_chroms
    remaining_snps = n_snps % n_chroms

    variant_idx = 0
    for chrom_idx in range(1, n_chroms + 1):
        # Distribute remaining SNPs across first chromosomes
        n_snps_this_chrom = snps_per_chrom + (1 if chrom_idx <= remaining_snps else 0)

        # Generate positions evenly spaced across chromosome
        positions = np.linspace(1000, chrom_length - 1000, n_snps_this_chrom, dtype=int)

        for pos in positions:
            variant_idx += 1
            chrom = f"chr{chrom_idx}"
            variant_id = f"{chrom}_{pos}"

            # Random quality (30-90)
            qual = int(np.random.uniform(30, 90))

            # Random allele frequency 
            af = np.random.uniform(0, 1)

            # Random total depth (50-200)
            total_dp = int(np.random.uniform(50, 200))

            # Random nucleotides for REF and ALT
            nucleotides = ['A', 'C', 'G', 'T']
            ref = np.random.choice(nucleotides)
            alt_choices = [n for n in nucleotides if n != ref]
            alt = np.random.choice(alt_choices)

            # Generate genotypes for all samples
            genotypes = []
            actual_alt_count = 0
            actual_total = 0

            for sample_idx in range(n_samples):
                # Determine if this genotype is missing based on sample's missing probability
                is_missing = np.random.random() < missing_probs[sample_idx]

                if is_missing:
                    gt = "./."
                    sample_dp = 0
                else:
                    # Generate genotype based on allele frequency
                    # Two alleles, each has probability af of being alternate
                    allele1 = 1 if np.random.random() < af else 0
                    allele2 = 1 if np.random.random() < af else 0
                    gt = f"{allele1}/{allele2}"

                    # Track actual AF in non-missing genotypes
                    actual_alt_count += allele1 + allele2
                    actual_total += 2

                    # Random read depth for this sample (10-50)
                    sample_dp = int(np.random.uniform(10, 50))

                genotypes.append(f"{gt}:{sample_dp}")

            # Calculate actual AF from generated genotypes
            if actual_total > 0:
                actual_af = actual_alt_count / actual_total
            else:
                actual_af = af  # All missing, use target AF

            # Build variant line
            info = f"DP={total_dp};AF={actual_af:.4f}"
            variant_line = f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\tPASS\t{info}\tGT:DP\t" + "\t".join(genotypes)
            vcf_lines.append(variant_line)

    logger.info(f"Generated {variant_idx} variants across {n_chroms} chromosomes")

    # Write uncompressed VCF
    if compress:
        uncompressed_path = output_path.with_suffix('')
        if uncompressed_path.suffix == '.vcf':
            pass  # Already .vcf
        else:
            uncompressed_path = output_path.with_suffix('.vcf')
    else:
        uncompressed_path = output_path
        if not str(uncompressed_path).endswith('.vcf'):
            uncompressed_path = Path(str(output_path) + '.vcf')

    logger.info(f"Writing VCF to {uncompressed_path}...")
    with open(uncompressed_path, 'w') as f:
        f.write('\n'.join(vcf_lines))
        f.write('\n')

    final_path = uncompressed_path

    # Compress and index if requested
    if compress:
        compressed_path = output_path if str(output_path).endswith('.gz') else Path(str(output_path) + '.gz')
        logger.info(f"Compressing with bcftools to {compressed_path}...")

        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", str(compressed_path), str(uncompressed_path)],
            check=True
        )

        logger.info("Indexing VCF...")
        subprocess.run(
            ["bcftools", "index", "-f", str(compressed_path)],
            check=True
        )

        # Remove uncompressed file
        uncompressed_path.unlink()
        final_path = compressed_path

    logger.info(f"Example VCF created: {final_path}")
    return str(final_path)


def main():
    """Command-line interface for generating example VCFs"""
    parser = argparse.ArgumentParser(
        description="Generate example VCF files for testing SNPSnip",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate default example (20 samples, 10k SNPs)
  python -m snpsnip.generate_example -o example.vcf.gz

  # Generate large example with more missing data
  python -m snpsnip.generate_example -o large.vcf.gz -n 100 -m 100000 --missing 0.15

  # Generate uncompressed VCF
  python -m snpsnip.generate_example -o example.vcf --no-compress

  # Reproducible example with seed
  python -m snpsnip.generate_example -o example.vcf.gz --seed 42
"""
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output VCF file path (use .vcf.gz for compressed)'
    )
    parser.add_argument(
        '-n', '--samples',
        type=int,
        default=20,
        help='Number of samples (default: 20)'
    )
    parser.add_argument(
        '-m', '--snps',
        type=int,
        default=10000,
        help='Number of SNPs (default: 10000)'
    )
    parser.add_argument(
        '--missing',
        type=float,
        default=0.05,
        help='Proportion of missing genotypes, exponentially distributed across samples (default: 0.05)'
    )
    parser.add_argument(
        '--chroms',
        type=int,
        default=2,
        help='Number of chromosomes (default: 2)'
    )
    parser.add_argument(
        '--no-compress',
        action='store_true',
        help='Do not compress and index the VCF'
    )
    parser.add_argument(
        '--seed',
        type=int,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s: %(message)s'
    )

    # Validate inputs
    if args.samples < 1:
        logger.error("Number of samples must be at least 1")
        return 1

    if args.snps < 1:
        logger.error("Number of SNPs must be at least 1")
        return 1

    if not 0 <= args.missing <= 1:
        logger.error("Missing proportion must be between 0 and 1")
        return 1

    try:
        vcf_path = generate_example_vcf(
            output_path=args.output,
            n_samples=args.samples,
            n_snps=args.snps,
            missing_prop=args.missing,
            n_chroms=args.chroms,
            compress=not args.no_compress,
            random_seed=args.seed
        )

        print(f"Successfully generated example VCF: {vcf_path}")
        return 0

    except Exception as e:
        logger.error(f"Error generating VCF: {e}")
        if args.verbose:
            raise
        return 1


if __name__ == '__main__':
    sys.exit(main())
