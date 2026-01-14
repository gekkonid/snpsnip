---
title: 'SNPSnip: Efficient interactive filtering and subsetting of genomic variant data'
tags:
  - Python
  - genomics
  - bioinformatics
  - variant calling
  - data filtering
authors:
  - name: [Author Name]
    orcid: [ORCID ID]
    affiliation: 1
affiliations:
 - name: [Institution Name]
   index: 1
date: [Date]
bibliography: paper.bib
---

# Summary

Genomic variation datasets nowadays contain thousands of samples genotyped at
millions of variants. Variant datasets require filtering and subsetting before
downstream analysis, however optimal filtering requires a degree of
interactivity incompatible with efficient batch processing. `SNPSnip` solves
this dilemma with the seamless combination of two aspects: fast, interactive
filtering on a subset of data, and computationally-efficient, parallelised
filtering of a large dataset once a user establishes optimal filter thresholds.
Using stand-alone, single-file, HTML interfaces, we separate efficient parallel
computation from interactivity, for example between a HPC server and a user's
laptop. SNPSnip is written in Python (with HTML/JS interactive interface), is
used via the command line and any modern web browser, and uses bcftools for
underlying bulk VCF/BCF processing.


# Statement of Need

Modern population and quantaitve genomics studies generate enormous variant
datasets, for example 1000 *Arabidopsis* samples at 6 million SNPs (@alonso
blanco). The variant callers used to generate such data do not produce
perfectly accurate calls, and different downstream analyses often have their
own optimal filtering and signal-to-noise tradeoffs. Compounding this,
industry-standard tools for variant filtering like BCFtools typically do not
scale linearly when given additional CPUs, as only some parts of the overall
computation are parallelised.

Typically, researchers will use a variant dataset in multiple downstream
analyses, perhaps with different degrees of stringency of sample or SNP
quality, or with different sample subsets. To determine the appropriate sample
groupings and sample/SNP quality thesholds, and select data-defined subsets, a
researcher must interactively plot quality measures and measures of sample
identity and relatedness (e.g. PCA). These statsitcs are typically calculated
from an aribtariy subset of a VCF, and statsitics are calcualted, plotted, and
deciions made in an interactive statisitcal environement e.g. R. Decisions must
then be ported back to filtering pipelines, which may run on remote servers or
bactch queue systems.

# Software Design

`SNPSnip` enables computationally-efficient interactive SNP filtering as follows:

1) Create a randomly-subset SNP dataset for interactive use, optionally applying baseline filters to samples and/or SNPs.
2) Generate summary plots (histograms of quality statstics and PCAs) that the user may use to filter samples, including optionally separating samples in to sub-groups by some metadata e.g. species.
3) Embed these plots in an interactive web page that allows selection of final sample sub-groups
4) Summarise SNP quality independenlty across each sample subset, and generate an interactive web page to set filters for each sample set.
5) Produce a final variant set containing all SNPs passing defined filters for each sample set. This computation is done region-parallel, scaling well up to very large numbers of CPUs (e.g. 256 or more) as each genome chunk is processed indpeendently.

Importantly, SNPSnip separates bulk processing from interactive analysis. When
user interaction is required, SNPSnip either opens an ephemeral webserver
(accssible on a remote machine via SSH port forwarding), or generates a static
HTML file that can be transferred e.g. using rsync. This file and resulting
selections are megabytes in size, even for datasets measured in terabytes. 

While industry-standard tools like bcftools offer some degree of parallelism,
this is typically limited to additional (de)compression worker threads. To
fully utilise the tens or hundreds of CPUs availble on modern HPC servers or
compute nodes, one must use region-parallism, where chunks of neighboring
variants are indpeendently processed on a single thread, with regions processed
in parallel across all available CPUs. The VCF/BCF format indexing makes such
region-parallelism highly efficient as each independent thread can effieiently
seek to its region. One then must concatenate the resulting chunks, however
this operation is relatively comptuationally efficient, and SNPSnip users are
provided guidance on an appropriate chunk size to best balance efficient
parallel filtering with the cost of merging too many files. BCFtools also
robustly handles variants on region borders, when merging regions back into a
single VCF/BCF.

# AI Usage Disclosure

The core code and documentation of SNPSnip was written by hand, however coding
of some javascript UI code and intergration tests, and some refactoring was
aided by LLMs, and all code was reviewed for bugs (and some bugs fixed) by
LLMs. We verfied all output of LLMs for correctness.



# References

