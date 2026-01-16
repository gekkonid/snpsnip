---
title: 'SNPSnip: Efficient interactive filtering and subsetting of genomic variant data'
tags:
  - Python
  - genomics
  - bioinformatics
  - variant calling
  - data filtering
authors:
  - name: Kevin D. Murray
    orcid: 0000-0002-2466-1917
    affiliation: "1,2"
affiliations:
 - name: Gekkonid Scientific Pty. Ltd., Melbourne, AU
   index: 1
 - name: Max Planck Institute for Biology, Tübingen, DE
   index: 2
date: 2026-01-15
bibliography: paper.bib
---

# Summary

Genomic variant datasets now reguarly contain thousands of samples genotyped at
millions of variants. Researchers must filter and subset these datasets before
downstream analysis. Optimal filtering requires interactive exploratory
analyses incompatible with efficient batch computation. `SNPSnip` solves this
dilemma by combining fast, interactive filtering on a data subset with
parallelised processing of full datasets once users establish optimal
thresholds. Stand-alone HTML interfaces separate efficient parallel computation
from interactivity, enabling work across HPC servers and local machines.
`SNPSnip` is implemented in Python (with HTML/JS interface), operates via
command line and web browser, and uses `bcftools` for VCF/BCF processing. Users
can obtain `SNPSnip` from Github (<https://github.com/gekkonid/snpsnip>) or
PyPI (`pip install snpsnip`) under the Mozilla Public License v2.


# Statement of Need

Modern population and quantitative genomics studies generate enormous variant
datasets, such as 1000 *Arabidopsis* samples at 6 million SNPs
[@alonso-blanco16_1001genomes]. Variant callers do not produce perfectly
accurate calls, so reserchers must filter raw variant calls to avoid erroneous
inferences. Furthermore, researchers typically use a variant dataset in
multiple downstream analyses, and different downstream analyses have their own
optimal filtering and signal-to-noise tradeoffs. Determining appropriate sample
groupings and quality thresholds requires interactive plotting of quality
measures and sample relatedness [@marees18_tutorialconducting;
@ahrens21_regardingfword]. Researchers typically calculate these statistics,
then plot and analyse them in an interactive environment such as R. They must
then port decisions back to scalable filtering pipelines, which may run on
remote servers or batch queue systems.


# State of the Field

Existing variant filtering tools broadly fall into two categories: efficient
batch processing tools such as PLINK [@chang15_secondgenerationplink], VCFtools
[@danecek11_variantcall], and BCFtools [@danecek21_bcftools12years]; or less
efficient tools for interactive exploration, such as the R packages
VariantAnnotation [@obenchain14_variantannotationbioconductor] and vcfR
[@knaus17_vcfrpackage;@knaus16_vcfrpackage]. Neither category fully addresses
user-friendly, interactive filtering of very large datasets. Batch tools
require predefined thresholds with no exploratory phase, while scripting
packages become inefficient with massive datasets. Notable exceptions include
packages designed for large VCF datasets, such as SNPRelate in R
[@zheng12_highperformancecomputing] or sgkit/VCFZarr in Python
[@czech25_analysisreadyvcf]. However, these tools require users to convert VCFs
to bespoke formats and typically constrain analyses to those implemented within
their framework.

`SNPSnip` addresses these limitations by integrating interactive threshold
selection with highly efficient, scalable region-based filtering, that produces
standard-compliant VCF files suitable for any downstream tool.


# Software Design

`SNPSnip` enables efficient interactive SNP filtering through five steps:

1) Create a random subset of variants for interactive purposes, optionally applying baseline filters.
2) Summarise sample quality and relatedness (via PCA) on this subset, optionally separating samples into groups by metadata (e.g. species).
3) Embed plots in an interactive web page for sample sub-group selections, forming refined groups of only high-quality samples.
4) Summarise SNP quality independently across each sample group, generating an interactive web page for filter selection.
5) Produce final variant sets containing all SNPs passing defined filters for each sample group, using region-parallel computation that scales to hundreds of CPUs.

`SNPSnip` separates bulk processing from interactive analysis. When users need
to interact, `SNPSnip` either opens an ephemeral webserver (accessible remotely
via SSH port forwarding) or generates a static HTML file to transfer to a local
computer for interaction. These files and resulting selections remain megabytes
in size, even for terabyte-scale datasets.

While tools like `bcftools` offer some parallelism, for example additional
(de)compression threads, fully utilising modern HPC resources requires
region-parallelism, where single threads independently process chunks of
neighboring variants across all available CPUs. VCF/BCF indexing makes this
approach highly efficient, as each thread can seek directly to its assigned
region, and concatenating the results of these independently regions is
typically realtively fast. `SNPSnip` provides guidance on appropriate chunk
sizes to balance parallel filtering efficiency with merging costs. `bcftools`
robustly handles variants on region borders when merging regions.

`SNPSnip` includes an automated test suite, ships as an easily-installable
Python package, and provides an interactive tutorial (also covered by automated
tests) with example data for learning.


# Research Impact Statement

We developed `SNPSnip` for our work in plant population genomics, where we needed
high-throughput interactive filtering of VCFs without moving data between
HPC clusters and users' off-site laptops. Despite its recent release and
minimal promotion, several research groups in quantitative, population, and
evolutionary genetics now use the tool.


# AI Usage Disclosure

We wrote the core code and documentation of SNPSnip by hand. LLMs aided coding
of some JavaScript UI code and integration tests, and assisted with
refactoring. LLMs reviewed all code for bugs, and we fixed some bugs with LLM
assistance. We verified all LLM output for correctness.


# References
