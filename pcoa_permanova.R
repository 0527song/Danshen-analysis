#!/usr/bin/env Rscript

# ------------------------------------------------------------
# pcoa_permanova.R
# - Compute Bray-Curtis distance from a feature table
# - PCoA ordination and plotting (color=group, shape=location)
# - PERMANOVA: location * tissue
# - Beta-dispersion tests (betadisper) for location and tissue
#
# Inputs:
#   --feature_table : TSV, rows=OTU/ASV, cols=samples, numeric counts/abundance
#   --metadata      : TSV with columns: sample, group (e.g., SXBS, SXRS, SXRT...)
# Outputs:
#   - pcoa_points.tsv
#   - pcoa_plot.pdf / pcoa_plot.png
#   - permanova_location_tissue.tsv
#   - permanova_location.tsv / permanova_tissue.tsv
#   - betadisper_location.tsv / betadisper_tissue.tsv
# ------------------------------------------------------------

# --------------------------
# 0) Simple arg parser
# --------------------------
parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--[^=]+=", a)) {
      key <- sub("^--([^=]+)=.*$", "\\1", a)
      val <- sub("^--[^=]+=(.*)$", "\\1", a)
      out[[key]] <- val
      i <- i + 1
    } else if (grepl("^--", a)) {
      key <- sub("^--", "", a)
      if (i == length(args)) stop("Missing value for argument: ", a)
      val <- args[i + 1]
      if (grepl("^--", val)) stop("Missing value for argument: ", a)
      out[[key]] <- val
      i <- i + 2
    } else {
      stop("Unrecognized argument: ", a)
    }
  }
  out
}

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "), "\n",
      "Install them before running. Example:\n",
      "  install.packages(c('ggplot2'))\n",
      "  install.packages(c('dplyr','readr'))\n",
      "  install.packages('vegan')\n",
      "  install.packages('ape')\n"
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

feature_table <- args[["feature_table"]]
metadata      <- args[["metadata"]]
outdir        <- args[["outdir"]]
distance_m    <- args[["distance"]]
ellipse_level <- args[["ellipse_level"]]
seed          <- args[["seed"]]

if (is.null(feature_table) || is.null(metadata) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/pcoa_permanova.R \\\n",
    "    --feature_table <feature_table.tsv> \\\n",
    "    --metadata <metadata.tsv> \\\n",
    "    --outdir <outdir> \\\n",
    "    [--distance bray] [--ellipse_level 0.68] [--seed 123]\n\n",
    "metadata.tsv must have columns: sample, group\n"
  )
  quit(status = 1)
}

if (!file.exists(feature_table)) stop("feature_table not found: ", feature_table)
if (!file.exists(metadata))      stop("metadata not found: ", metadata)

if (is.null(distance_m)) distance_m <- "bray"
if (is.null(ellipse_level)) ellipse_level <- "0.68"
if (is.null(seed)) seed <- "123"

ellipse_level <- as.numeric(ellipse_level)
seed <- as.integer(seed)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 1) Load packages
# --------------------------
need(c("ggplot2", "vegan", "ape", "dplyr"))

library(ggplot2)
library(vegan)
library(ape)
library(dplyr)

set.seed(seed)

# --------------------------
# 2) Read inputs
# --------------------------
cat("[1/6] Reading feature table and metadata...\n")

otu <- read.table(
  feature_table,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

meta <- read.table(
  metadata,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("sample", "group") %in% colnames(meta))) {
  stop("metadata must contain columns: sample, group")
}

# Keep intersecting samples and align order
samples_common <- intersect(colnames(otu), meta$sample)
if (length(samples_common) < 3) {
  stop("Too few overlapping samples between feature_table columns and metadata$sample.")
}

otu <- otu[, samples_common, drop = FALSE]
meta <- meta[match(samples_common, meta$sample), , drop = FALSE]

# Basic checks
if (any(is.na(meta$group))) stop("Some samples in metadata could not be matched/aligned; check sample names.")
if (any(otu < 0, na.rm = TRUE)) stop("feature_table contains negative values; cannot compute Bray distance.")
if (anyNA(otu)) stop("feature_table contains NA; please impute/remove before analysis.")

# --------------------------
# 3) Distance + PCoA
# --------------------------
cat("[2/6] Computing distance and running PCoA...\n")

# Samples x features matrix for vegan::vegdist
comm <- t(as.matrix(otu))

# If counts, it's generally safer to standardize to relative abundance for Bray
comm_rel <- vegan::decostand(comm, method = "total")

dist_obj <- vegan::vegdist(comm_rel, method = distance_m)

# PCoA (ape::pcoa)
pcoa_res <- ape::pcoa(dist_obj)
points <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
colnames(points) <- c("PC1", "PC2")
points$sample <- rownames(points)

# eigenvalues for variance explained (use positive eigenvalues sum)
eig <- pcoa_res$values$Eigenvalues
eig_pos_sum <- sum(eig[eig > 0])
pc1_pct <- 100 * eig[1] / eig_pos_sum
pc2_pct <- 100 * eig[2] / eig_pos_sum

# Attach metadata
points <- points %>%
  left_join(meta, by = c("sample" = "sample"))

# Derive location/tissue from group like "SXBS"
points$location <- substr(points$group, 1, 2)  # SX, YY, YZ
points$tissue   <- substr(points$group, 3, 4)  # BS, RS, RT

# Shape group = location (your original逻辑)
points$shape_group <- factor(points$location, levels = c("SX", "YY", "YZ"))

# Save points
points_out <- file.path(outdir, "pcoa_points.tsv")
write.table(points, points_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: ", points_out, "\n", sep = "")

# --------------------------
# 4) Plot (color=group, shape=location)
# --------------------------
cat("[3/6] Plotting PCoA...\n")

# Default color palette for 9 groups (edit as needed)
group_palette <- c(
  "SXBS" = "#1b9e77", "SXRS" = "#d95f02", "SXRT" = "#7570b3",
  "YYBS" = "#e7298a", "YYRS" = "#66a61e", "YYRT" = "#e6ab02",
  "YZBS" = "#a6761d", "YZRS" = "#666666", "YZRT" = "#ff7f00"
)

# If there are groups not in palette, fall back to ggplot default colors
use_manual_colors <- all(unique(points$group) %in% names(group_palette))

p <- ggplot(points, aes(x = PC1, y = PC2, color = group, shape = shape_group)) +
  geom_point(alpha = 0.8, size = 3) +
  stat_ellipse(level = ellipse_level, linetype = 2) +
  labs(
    x = paste0("PCoA 1 (", format(pc1_pct, digits = 4), "%)"),
    y = paste0("PCoA 2 (", format(pc2_pct, digits = 4), "%)"),
    title = paste0(toupper(distance_m), " PCoA")
  ) +
  theme_classic() +
  scale_shape_manual(values = c(SX = 16, YY = 17, YZ = 15))

if (use_manual_colors) {
  p <- p + scale_color_manual(values = group_palette)
}

pdf_file <- file.path(outdir, "pcoa_plot.pdf")
png_file <- file.path(outdir, "pcoa_plot.png")

ggsave(pdf_file, p, width = 7.5, height = 6, limitsize = FALSE)
ggsave(png_file, p, width = 7.5, height = 6, dpi = 300, limitsize = FALSE)

cat("  Saved:\n")
cat("   - ", pdf_file, "\n", sep = "")
cat("   - ", png_file, "\n", sep = "")

# --------------------------
# 5) PERMANOVA (adonis2)
# --------------------------
cat("[4/6] Running PERMANOVA (adonis2)...\n")

sub_design <- meta
sub_design$location <- factor(substr(sub_design$group, 1, 2), levels = c("SX", "YY", "YZ"))
sub_design$tissue   <- factor(substr(sub_design$group, 3, 4), levels = c("BS", "RS", "RT"))

# Ensure rownames align to distance labels
rownames(sub_design) <- sub_design$sample
sub_design <- sub_design[labels(dist_obj), , drop = FALSE]

# Sequential terms (same behavior as你的示例注释)
perm_full <- vegan::adonis2(dist_obj ~ location * tissue, data = sub_design, permutations = 999)

perm_loc  <- vegan::adonis2(dist_obj ~ location, data = sub_design, permutations = 999)
perm_tis  <- vegan::adonis2(dist_obj ~ tissue,   data = sub_design, permutations = 999)

write.table(as.data.frame(perm_full),
            file.path(outdir, "permanova_location_tissue.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(as.data.frame(perm_loc),
            file.path(outdir, "permanova_location.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(as.data.frame(perm_tis),
            file.path(outdir, "permanova_tissue.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("  Saved PERMANOVA tables under: ", outdir, "\n", sep = "")

# --------------------------
# 6) Beta dispersion (betadisper)
# --------------------------
cat("[5/6] Testing beta dispersion (betadisper)...\n")

disp_location <- vegan::betadisper(dist_obj, sub_design$location)
disp_tissue   <- vegan::betadisper(dist_obj, sub_design$tissue)

bd_loc <- anova(disp_location)
bd_tis <- anova(disp_tissue)

write.table(as.data.frame(bd_loc),
            file.path(outdir, "betadisper_location.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(as.data.frame(bd_tis),
            file.path(outdir, "betadisper_tissue.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("[6/6] Done.\n")
