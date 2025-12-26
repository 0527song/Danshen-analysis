#!/usr/bin/env Rscript

# ------------------------------------------------------------
# FigB: Top N genus abundance heatmap (DOR vs NOR)
# Input: CSV/TSV where rows=Genus/ID and columns=samples (e.g., DOR1..DOR4, NOR1..NOR4)
# Output: PDF + PNG heatmap
#
# Features:
# - Select top N taxa by mean abundance across all samples
# - Optional log10(x + pseudocount)
# - Row-wise Z-score (scale="row") and cap to [-cap_z, cap_z] to match legend like -2..2
# - Cluster rows, keep columns in fixed order (DOR then NOR)
# - Column annotation bar (Group: DOR/NOR)
# ------------------------------------------------------------

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--[^=]+=", a)) {
      k <- sub("^--([^=]+)=.*$", "\\1", a)
      v <- sub("^--[^=]+=(.*)$", "\\1", a)
      out[[k]] <- v; i <- i + 1
    } else if (grepl("^--", a)) {
      k <- sub("^--", "", a)
      if (i == length(args)) stop("Missing value for: ", a)
      v <- args[i + 1]
      if (grepl("^--", v)) stop("Missing value for: ", a)
      out[[k]] <- v; i <- i + 2
    } else stop("Unrecognized argument: ", a)
  }
  out
}

need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop(
      "Missing R packages: ", paste(miss, collapse = ", "), "\n",
      "Install, e.g.\n",
      "  install.packages(c('pheatmap','RColorBrewer'))\n"
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

input  <- args[["input"]]
outdir <- args[["outdir"]]
top_n  <- args[["top_n"]]
cap_z  <- args[["cap_z"]]
log10t <- args[["log10"]]
pseudocount <- args[["pseudocount"]]
width  <- args[["width"]]
height <- args[["height"]]

if (is.null(input) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/figB_top30_genus_heatmap.R \\\n",
    "    --input <streptomyces_abundance.csv> \\\n",
    "    --outdir <outdir> \\\n",
    "    [--top_n 30] [--cap_z 2] [--log10 TRUE] [--pseudocount 1e-6] [--width 6] [--height 7]\n",
    sep = ""
  )
  quit(status = 1)
}

if (!file.exists(input)) stop("Input not found: ", input)

if (is.null(top_n)) top_n <- "30"
if (is.null(cap_z)) cap_z <- "2"
if (is.null(log10t)) log10t <- "TRUE"
if (is.null(pseudocount)) pseudocount <- "1e-6"
if (is.null(width)) width <- "6"
if (is.null(height)) height <- "7"

top_n <- as.integer(top_n)
cap_z <- as.numeric(cap_z)
pseudocount <- as.numeric(pseudocount)
width <- as.numeric(width)
height <- as.numeric(height)

need(c("pheatmap", "RColorBrewer"))
library(pheatmap)
library(RColorBrewer)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 1) Read data (CSV/TSV)
# --------------------------
is_tsv <- grepl("\\.tsv$|\\.txt$", tolower(input))
df <- if (is_tsv) {
  read.table(input, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
} else {
  read.csv(input, header = TRUE, check.names = FALSE)
}

# If first column is ID, use it as rownames
if ("ID" %in% colnames(df)) {
  rownames(df) <- df$ID
  df$ID <- NULL
} else if (!is.null(df[[1]]) && !anyDuplicated(df[[1]]) && !is.numeric(df[[1]])) {
  # fallback: treat first column as ID if it looks like names
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
}

# Ensure numeric matrix
mat <- as.matrix(df)
storage.mode(mat) <- "numeric"

if (nrow(mat) < 2 || ncol(mat) < 2) stop("Input table is too small to plot a heatmap.")

# --------------------------
# 2) Column order & group annotation (DOR then NOR)
# --------------------------
cols <- colnames(mat)

dor_cols <- cols[grepl("^DOR", cols, ignore.case = FALSE)]
nor_cols <- cols[grepl("^NOR", cols, ignore.case = FALSE)]

if (length(dor_cols) == 0 || length(nor_cols) == 0) {
  stop("Column names must include DOR* and NOR* samples (e.g., DOR1.., NOR1..).")
}

# Keep order DOR1.. then NOR1..
# If suffix numbers exist, order by numeric part
order_by_suffix <- function(x) {
  suf <- suppressWarnings(as.numeric(gsub("^[A-Za-z]+", "", x)))
  if (all(!is.na(suf))) x[order(suf)] else x
}
dor_cols <- order_by_suffix(dor_cols)
nor_cols <- order_by_suffix(nor_cols)

mat <- mat[, c(dor_cols, nor_cols), drop = FALSE]

annotation_col <- data.frame(
  Group = factor(c(rep("DOR", length(dor_cols)), rep("NOR", length(nor_cols))),
                 levels = c("DOR", "NOR"))
)
rownames(annotation_col) <- colnames(mat)

ann_colors <- list(
  Group = c(DOR = "#f26b6b", NOR = "#39b7b5")  # 贴近你图里的红/青
)

# --------------------------
# 3) Select Top N genera by mean abundance
# --------------------------
mean_ab <- rowMeans(mat, na.rm = TRUE)
ord <- order(mean_ab, decreasing = TRUE)

top_n_use <- min(top_n, nrow(mat))
mat_top <- mat[ord[1:top_n_use], , drop = FALSE]

# --------------------------
# 4) Transform + row Z-score + cap
# --------------------------
if (toupper(log10t) == "TRUE") {
  mat_top <- log10(mat_top + pseudocount)
}

# Row-wise Z-score
z <- t(scale(t(mat_top)))
z[is.na(z)] <- 0  # handle zero-variance rows

# Cap to [-cap_z, cap_z] to resemble legend -2..2
z <- pmax(pmin(z, cap_z), -cap_z)

# --------------------------
# 5) Color palette (blue - yellow - red)
# --------------------------
# Similar to RdYlBu reversed (blue->yellow->red)
col_fun <- colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(101)

# Set breaks exactly to show -cap_z..cap_z
breaks <- seq(-cap_z, cap_z, length.out = length(col_fun) + 1)

# --------------------------
# 6) Plot & save
# --------------------------
pdf_file <- file.path(outdir, "FigB_top30_genus_heatmap.pdf")
png_file <- file.path(outdir, "FigB_top30_genus_heatmap.png")

# PDF
pdf(pdf_file, width = width, height = height)
pheatmap(
  z,
  color = col_fun,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,             # keep DOR then NOR
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  border_color = NA,
  main = "TOP 30 genus abundance heatmap"
)
dev.off()

# PNG
png(png_file, width = width * 300, height = height * 300, res = 300)
pheatmap(
  z,
  color = col_fun,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  border_color = NA,
  main = "TOP 30 genus abundance heatmap"
)
dev.off()

# Save the processed matrix for reproducibility
write.table(z, file.path(outdir, "FigB_heatmap_matrix_zscore.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("Done:\n", pdf_file, "\n", png_file, "\n", sep = "")
