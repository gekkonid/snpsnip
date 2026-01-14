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
underlying bulk VCF/BCF processing. It is available from Github and PyPI under
the terms of the Mozilla Public License v2.


# Statement of Need

Modern population and quantitative genomics studies generate enormous variant
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
groupings and sample/SNP quality thresholds, and select data-defined subsets, a
researcher must interactively plot quality measures and measures of sample
identity and relatedness (e.g. PCA). These statistics are typically calculated
from an arbitrary subset of a VCF, and statistics are calculated, plotted, and
decisions made in an interactive statistical environment e.g. R. Decisions must
then be ported back to filtering pipelines, which may run on remote servers or
batch queue systems.

# State of the Field

Existing tools and approaches for variant filtering are typically either tools
for efficient batch processing of variant data, for example PLINK, VCFtools and
BCFtools, or generally less efficient tools for interactive exploration of
variant datasets, e.g. R packages VariantAnnotation and vcfR. Neither of these
classes of tool fully address the problem of user-friendly, interactive
filtering of very large datases, as they either simply filter based on
predefined thresholds with no interacrtive/exploratory phase (in the case of
tools like BCFtools and PLINK), or become inefficient with truy massive
datasets (in the case of the various scripting pacakges for interactive
explorations of VCFs). Some notable excptions are scripting packages designed
specifivally for very large VCF datasets, for example SNPRelate in R, or
sgkit/VCFZarr/bio2zarr in python. However, such tools achive their efficiency
by converting VCFs to a bespoke format, and typically provide the user with
tooling to perform certain analyses in their framework.


`SNPSnip` addresses these limitations by providing an integrated workflow that
combines interactive threshold selection with highly efficient and scalable
region-based filtering.


# Software Design

`SNPSnip` enables computationally-efficient interactive SNP filtering as follows:

1) Create a randomly-subset SNP dataset for interactive use, optionally applying baseline filters to samples and/or SNPs.
2) Generate summary plots (histograms of quality statstics and PCAs) that the user may use to filter samples, including optionally separating samples in to sub-groups by some metadata e.g. species.
3) Embed these plots in an interactive web page that allows selection of final sample sub-groups
4) Summarise SNP quality independently across each sample subset, and generate an interactive web page to set filters for each sample set.
5) Produce a final variant set containing all SNPs passing defined filters for each sample set. This computation is done region-parallel, scaling well up to very large numbers of CPUs (e.g. 256 or more) as each genome chunk is processed indpeendently.

Importantly, SNPSnip separates bulk processing from interactive analysis. When
user interaction is required, SNPSnip either opens an ephemeral webserver
(accessible on a remote machine via SSH port forwarding), or generates a static
HTML file that can be transferred e.g. using rsync. This file and resulting
selections are megabytes in size, even for datasets measured in terabytes. 

While industry-standard tools like bcftools offer some degree of parallelism,
this is typically limited to additional (de)compression worker threads. To
fully utilise the tens or hundreds of CPUs available on modern HPC servers or
compute nodes, one must use region-parallelism, where chunks of neighboring
variants are independently processed on a single thread, with regions processed
in parallel across all available CPUs. The VCF/BCF format indexing makes such
region-parallelism highly efficient as each independent thread can effieiently
seek to its region. One then must concatenate the resulting chunks, however
this operation is relatively computationally efficient, and SNPSnip users are
provided guidance on an appropriate chunk size to best balance efficient
parallel filtering with the cost of merging too many files. BCFtools also
robustly handles variants on region borders, when merging regions back into a
single VCF/BCF.

SNPSnip has an automated test suite that validates its operations, and is
packaged as an easily-installable python package. SNPSnip comes with an
interactive tutorial (that it itself also covered by automated tests),
including code to generate an example dataset so users have a small "sandpit"
dataset with which to learn how SNPSnip operates.


# Research Impact statement


The deveopment of SNPsnip was driven by our own work in plant popualtion
genomics, for which we need high-througput interactive filtering where data did
not need to move between remote workers and an HPC cluster. Despite being both
relatively new and so far quite poorly promoted, we are aware of several
research groups in quantiative, population, and evolutionary genetics using the
tool.



# AI Usage Disclosure

The core code and documentation of SNPSnip was written by hand, however coding
of some javascript UI code and integration tests, and some refactoring was
aided by LLMs, and all code was reviewed for bugs (and some bugs fixed) by
LLMs. We verified all output of LLMs for correctness.



# References
