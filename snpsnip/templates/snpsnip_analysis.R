#' ---
#' title: "SNPsnip R Analysis"
#' format: html
#' execute:
#'   warning: false
#' ---

#+ setup, include=FALSE
library(tidyverse)
library(hdf5r)   # install.packages("hdf5r")

h5 = H5File$new("snpsnip_data.h5", mode = "r")

samples         = h5[["samples"]]$read()
available_stats = h5[["available_stats"]]$read()
has_dp_matrix   = as.logical(h5[["has_dp_matrix"]]$read())

sample_stats = tibble(
    id           = h5[["sample_stats/id"]]$read(),
    missing_rate = h5[["sample_stats/missing_rate"]]$read(),
    mean_depth   = h5[["sample_stats/mean_depth"]]$read(),
    het_rate     = h5[["sample_stats/het_rate"]]$read()
)

# SNP info — always-present fields
snp_info = tibble(
    chrom     = h5[["snp_info/chrom"]]$read(),
    pos       = h5[["snp_info/pos"]]$read(),
    ref       = h5[["snp_info/ref"]]$read(),
    alt       = h5[["snp_info/alt"]]$read(),
    qual      = h5[["snp_info/qual"]]$read(),
    f_missing = h5[["snp_info/f_missing"]]$read()
)
# Optional SNP info fields (present only when the VCF contained them)
for (fld in c("dp", "af", "ac", "exhet")) {
    path = paste0("snp_info/", fld)
    if (h5$exists(path)) snp_info[[fld]] = h5[[path]]$read()
}

# PCA — HDF5 stores (n_samples, n_comp) in C order; hdf5r returns (n_comp, n_samples); t() fixes it
var_explained = h5[["pca/variance_explained"]]$read()
pca_data = as_tibble(
    t(h5[["pca/coordinates"]]$read()),
    .name_repair = ~ paste0("pc", seq_along(.))
) |> mutate(sample = h5[["pca/samples"]]$read(), .before = 1)

# Genotype matrix — same C/Fortran transpose: HDF5 (n_snps, n_samples) → hdf5r (n_samples, n_snps) → t()
# Encoding: 0=hom-ref, 1=het, 2=hom-alt, -1→NA=missing
gt_matrix = t(h5[["gt_matrix"]]$read())
gt_matrix[gt_matrix == -1L] = NA_integer_
colnames(gt_matrix) = samples

if (has_dp_matrix) {
    dp_matrix = t(h5[["dp_matrix"]]$read())
    dp_matrix[dp_matrix == -1L] = NA_integer_
    colnames(dp_matrix) = samples
} else {
    dp_matrix = NULL
}

h5$close_all()
rm(h5)

#' ## Section 1: Sample QC Thresholds
#'
#' **Edit the threshold values below**, then re-render to update the plots and
#' sample counts. Samples failing these thresholds will be excluded from all groups.

#+ sample-qc-thresholds
MAX_MISSING_RATE = 0.2    # exclude samples with missing rate  > this
MIN_MEAN_DEPTH   = 5.0    # exclude samples with mean depth    < this
MAX_HET_RATE     = 0.5    # exclude samples with het rate      > this

#+ sample-qc-histograms
p_missing_sample = ggplot(sample_stats, aes(x = missing_rate)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = MAX_MISSING_RATE, color = "red", linetype = "dashed",
               linewidth = 1) +
    labs(title = "Missing Rate per Sample",
         subtitle = sprintf("Red line = MAX_MISSING_RATE (%.2f)", MAX_MISSING_RATE),
         x = "Missing Rate", y = "Count") +
    theme_minimal()

p_depth_sample = ggplot(sample_stats, aes(x = mean_depth)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = MIN_MEAN_DEPTH, color = "blue", linetype = "dashed",
               linewidth = 1) +
    labs(title = "Mean Depth per Sample",
         subtitle = sprintf("Blue line = MIN_MEAN_DEPTH (%.1f)", MIN_MEAN_DEPTH),
         x = "Mean Depth", y = "Count") +
    theme_minimal()

p_het_sample = ggplot(sample_stats, aes(x = het_rate)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = MAX_HET_RATE, color = "red", linetype = "dashed",
               linewidth = 1) +
    labs(title = "Heterozygosity Rate per Sample",
         subtitle = sprintf("Red line = MAX_HET_RATE (%.2f)", MAX_HET_RATE),
         x = "Heterozygosity Rate", y = "Count") +
    theme_minimal()

print(p_missing_sample)
print(p_depth_sample)
print(p_het_sample)

passing_samples = sample_stats |>
    filter(
        missing_rate <= MAX_MISSING_RATE,
        mean_depth   >= MIN_MEAN_DEPTH,
        het_rate     <= MAX_HET_RATE
    ) |>
    pull(id)

message(sprintf("%d / %d samples pass QC thresholds",
                length(passing_samples), nrow(sample_stats)))

#' ## Section 2: Define Sample Groups
#'
#' **Edit the `groups` list below** to assign samples to groups.
#' Each entry is a named character vector of sample IDs. You can do this
#' however you wish, for example reading from a metadata file or similar.
#'
#' Only samples that pass the QC thresholds in Section 1 are retained.
#'
#' Example with two groups:
#' ```r
#' groups = list(
#'   "population_A" = c("sample1", "sample2", "sample3"),
#'   "population_B" = passing_samples   # all remaining passing samples
#' )
#' ```

#+ define-groups
groups = list(
    "all" = passing_samples   # all passing samples end up here
)

#+ apply-sample-thresholds
# Filter all group members to passing samples only
groups = lapply(groups, intersect, passing_samples)

group_summary = tibble(
    group     = names(groups),
    n_samples = sapply(groups, length)
)
knitr::kable(group_summary, caption = "Sample Group Sizes")

#' ## Section 3: PCA Plot

#+ pca-plot
group_assignments = imap_dfr(groups, \(samps, grp) tibble(id = samps, group = grp))
pca_plot_data = pca_data |>
    left_join(group_assignments, by = c("sample" = "id"))

p_pca = ggplot(pca_plot_data, aes(x = pc1, y = pc2, color = group)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_brewer(palette = "Set1", na.value = "grey50") +
    labs(
        title = "PCA of SNP Genotypes",
        x     = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y     = sprintf("PC2 (%.1f%%)", var_explained[2]),
        color = "Group"
    ) +
    theme_minimal()
print(p_pca)

#' *Tip: join your own metadata (population, batch, etc.) to `pca_data` before
#' plotting to add additional colour/shape aesthetics.*

#' ## Section 4: Per-Group SNP Statistics
#'
#' **Edit the threshold constants below**, then re-render to see threshold lines
#' on the distributions. Use `NULL` to leave a bound unset (no filter).
#'
#' Notes:
#'
#' - `missing`: fraction of samples missing a call (0–1); per group
#' - `af`: allele frequency computed per group (0–1)
#' - `depth`: mean per-sample depth; per group (only if depth data is available)
#' - `qual`: VCF QUAL score (site-level, same for all groups)
#' - `exhet`: excess heterozygosity in **-log10** scale (site-level, same for all groups)
#' - `ac`: allele count (site-level)
#'
#' Blue lines = minimum bound; red lines = maximum bound.

#+ snp-filter-thresholds
MIN_SNP_MISSING = NULL;  MAX_SNP_MISSING = 0.2
MIN_SNP_MAF     = 0.05;  MAX_SNP_MAF     = NULL
MIN_SNP_DEPTH   = NULL;  MAX_SNP_DEPTH   = NULL
MIN_SNP_QUAL    = NULL;  MAX_SNP_QUAL    = NULL
MIN_SNP_EXHET   = NULL;  MAX_SNP_EXHET   = NULL   # -log10 scale
MIN_SNP_AC      = NULL;  MAX_SNP_AC      = NULL

#+ snp-stat-helpers
# Add optional threshold vlines (blue = lower bound, red = upper bound; NULL = skip)
add_vlines = function(plot, lo = NULL, hi = NULL) {
    if (!is.null(lo)) plot = plot + geom_vline(xintercept = lo, color = "blue",
                                               linetype = "dashed", linewidth = 1)
    if (!is.null(hi)) plot = plot + geom_vline(xintercept = hi, color = "red",
                                               linetype = "dashed", linewidth = 1)
    plot
}

apply_bound = function(pass, vec, lo, hi) {
    if (!is.null(lo)) pass = pass & !is.na(vec) & (vec >= lo)
    if (!is.null(hi)) pass = pass & !is.na(vec) & (vec <= hi)
    pass
}

#+ per-group-snp-stats
# Process one group at a time to keep peak RAM to one matrix slice at a time.
# Only the three per-group statistics are accumulated (not full snp_info rows).
grp_stats_list = vector("list", length(groups))
names(grp_stats_list) = names(groups)

for (grp in names(groups)) {
    grp_samples = groups[[grp]]
    if (length(grp_samples) == 0L) {
        grp_stats_list[[grp]] = tibble()
        next
    }
    grp_cols = intersect(grp_samples, colnames(gt_matrix))

    # Avoid copying the full matrix when the group contains all samples
    gt_grp = if (length(grp_cols) == ncol(gt_matrix)) gt_matrix
             else gt_matrix[, grp_cols, drop = FALSE]

    af = rowMeans(gt_grp, na.rm = TRUE) / 2
    out = tibble(
        group   = grp,
        missing = rowMeans(is.na(gt_grp)),
        maf     = pmin(af, 1 - af)
    )
    rm(gt_grp)

    if (!is.null(dp_matrix)) {
        dp_grp = if (length(grp_cols) == ncol(dp_matrix)) dp_matrix
                 else dp_matrix[, grp_cols, drop = FALSE]
        out$depth = rowMeans(dp_grp, na.rm = TRUE)
        rm(dp_grp)
    }
    grp_stats_list[[grp]] = out
}
grp_stats = bind_rows(grp_stats_list)
rm(grp_stats_list)

# --- Faceted plots (per-group stats) ---

p_missing = ggplot(grp_stats, aes(x = missing)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    labs(title = "Missing Rate per SNP", x = "Missing Rate", y = "Count") +
    theme_minimal()
p_missing = add_vlines(p_missing, MIN_SNP_MISSING, MAX_SNP_MISSING)
print(p_missing)

p_maf = ggplot(grp_stats, aes(x = maf)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    labs(title = "Minor Allele Frequency per SNP", x = "MAF", y = "Count") +
    theme_minimal()
p_maf = add_vlines(p_maf, MIN_SNP_MAF, MAX_SNP_MAF)
print(p_maf)

if (!is.null(dp_matrix)) {
    p_depth = ggplot(grp_stats, aes(x = depth)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        facet_wrap(~ group, ncol = 1, scales = "free_y") +
        labs(title = "Mean Depth per SNP", x = "Mean Depth per Sample", y = "Count") +
        theme_minimal()
    p_depth = add_vlines(p_depth, MIN_SNP_DEPTH, MAX_SNP_DEPTH)
    print(p_depth)
}

# --- Site-level plots (identical across groups, shown once) ---

if ("qual" %in% names(snp_info)) {
    p_qual = ggplot(snp_info, aes(x = qual)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        labs(title = "QUAL (site-level)", x = "QUAL", y = "Count") +
        theme_minimal()
    p_qual = add_vlines(p_qual, MIN_SNP_QUAL, MAX_SNP_QUAL)
    print(p_qual)
}

if ("exhet" %in% available_stats && "exhet" %in% names(snp_info)) {
    exhet_df = tibble(exhet_log10 = -log10(pmax(snp_info$exhet, .Machine$double.eps))) |>
        filter(is.finite(exhet_log10))
    p_exhet = ggplot(exhet_df, aes(x = exhet_log10)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
        labs(title = "Excess Heterozygosity (site-level)",
             x = "-log10(ExcHet)", y = "Count") +
        theme_minimal()
    p_exhet = add_vlines(p_exhet, MIN_SNP_EXHET, MAX_SNP_EXHET)
    print(p_exhet)
}

# --- Pass count preview ---

for (grp in names(groups)) {
    rows = filter(grp_stats, group == grp)
    if (nrow(rows) == 0L) next
    pass = rep(TRUE, nrow(rows))
    pass = apply_bound(pass, rows$missing, MIN_SNP_MISSING, MAX_SNP_MISSING)
    pass = apply_bound(pass, rows$maf,     MIN_SNP_MAF,     MAX_SNP_MAF)
    if ("depth" %in% names(rows)) {
        pass = apply_bound(pass, rows$depth, MIN_SNP_DEPTH, MAX_SNP_DEPTH)
    }
    if ("qual" %in% names(snp_info)) {
        pass = apply_bound(pass, snp_info$qual, MIN_SNP_QUAL, MAX_SNP_QUAL)
    }
    if ("exhet" %in% names(snp_info)) {
        pass = apply_bound(pass, -log10(pmax(snp_info$exhet, .Machine$double.eps)),
                           MIN_SNP_EXHET, MAX_SNP_EXHET)
    }
    if ("ac" %in% names(snp_info)) {
        pass = apply_bound(pass, snp_info$ac, MIN_SNP_AC, MAX_SNP_AC)
    }
    message(sprintf("Group '%s': %d / %d SNPs pass filters (%.1f%%)",
                    grp, sum(pass, na.rm = TRUE), nrow(rows),
                    100 * mean(pass, na.rm = TRUE)))
}

#' ## Section 5: Write Selections JSON
#'
#' When you are satisfied with your groups (Section 2) and thresholds (Section 4),
#' run this section to write `snpsnip_selections.json`, then copy it back to
#' the cluster and run:
#' ```
#' snpsnip --vcf <your.vcf.gz> --output-dir <output-dir/> --next snpsnip_selections.json
#' ```

#+ write-json
snp_thresholds = lapply(names(groups), function(grp) {
    thr = list(
        missing = list(min = MIN_SNP_MISSING, max = MAX_SNP_MISSING),
        maf     = list(min = MIN_SNP_MAF,     max = MAX_SNP_MAF)
    )
    if (!is.null(dp_matrix)) thr$depth = list(min = MIN_SNP_DEPTH, max = MAX_SNP_DEPTH)
    if ("exhet" %in% available_stats) thr$exhet = list(min = MIN_SNP_EXHET, max = MAX_SNP_EXHET)
    if ("qual"  %in% names(snp_info)) thr$qual  = list(min = MIN_SNP_QUAL,  max = MAX_SNP_QUAL)
    if ("ac"    %in% available_stats) thr$ac     = list(min = MIN_SNP_AC,    max = MAX_SNP_AC)
    thr
})
names(snp_thresholds) = names(groups)

output_file = "snpsnip_selections.json"
write_json(
    list(groups = groups, thresholds = snp_thresholds),
    output_file,
    auto_unbox = TRUE,
    null = "null",
    pretty = TRUE
)
message(sprintf("Written: %s", output_file))
message(sprintf("Next step:\n  snpsnip --vcf <your.vcf.gz> --output-dir <output-dir/> --next %s",
                output_file))
