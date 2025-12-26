#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Fig5A: Similarity matrix heatmap (ComplexHeatmap)
# - Input: similarity_matrix.csv (rows=IDs, cols=IDs, numeric similarity)
# - Type mapping: type_map.tsv (id, type)
# - Optional score: score.tsv (id, score) for right barplot
# - Output: Fig5A.pdf / Fig5A.png + ordered tables
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
      "  install.packages(c('circlize'))\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('ComplexHeatmap')\n"
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

matrix_file <- args[["matrix"]]
type_file   <- args[["type_map"]]
score_file  <- args[["score"]]      # optional
outdir      <- args[["outdir"]]
title       <- args[["title"]]

col_min <- args[["col_min"]]
col_mid <- args[["col_mid"]]
col_max <- args[["col_max"]]

width  <- args[["width"]]
height <- args[["height"]]

if (is.null(matrix_file) || is.null(type_file) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/fig5a_similarity_heatmap.R \\\n",
    "    --matrix <similarity_matrix.csv> \\\n",
    "    --type_map <type_map.tsv> \\\n",
    "    --outdir <outdir> \\\n",
    "    [--score <score.tsv>] \\\n",
    "    [--title \"BGC similarity\"] \\\n",
    "    [--col_min 0 --col_mid 0.5 --col_max 1] \\\n",
    "    [--width 9 --height 6]\n",
    sep = ""
  )
  quit(status = 1)
}

if (!file.exists(matrix_file)) stop("Matrix file not found: ", matrix_file)
if (!file.exists(type_file))   stop("Type map file not found: ", type_file)
if (!is.null(score_file) && !file.exists(score_file)) stop("Score file not found: ", score_file)

if (is.null(title)) title <- "Similarity"
if (is.null(width)) width <- "9"
if (is.null(height)) height <- "6"
width <- as.numeric(width)
height <- as.numeric(height)

need(c("ComplexHeatmap", "circlize"))
library(ComplexHeatmap)
library(circlize)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 1) Read similarity matrix
# --------------------------
sim <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
sim <- as.matrix(sim)
storage.mode(sim) <- "numeric"

if (nrow(sim) != ncol(sim)) {
  stop("Similarity matrix must be square: nrow != ncol")
}
if (!all(rownames(sim) == colnames(sim))) {
  # Not strictly required, but strongly recommended for similarity matrices
  warning("Row names and column names are not identical or not in the same order. Proceeding, but please verify.")
}

# --------------------------
# 2) Read type map (id, type)
# --------------------------
type_map <- read.table(type_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("id", "type") %in% colnames(type_map))) {
  stop("type_map must contain columns: id, type")
}

# Ensure all IDs exist
missing_ids <- setdiff(rownames(sim), type_map$id)
if (length(missing_ids) > 0) {
  stop("These matrix row IDs are missing in type_map.tsv (id column):\n", paste(missing_ids, collapse = ", "))
}

type_vector <- setNames(type_map$type[match(rownames(sim), type_map$id)], rownames(sim))

# --------------------------
# 3) Score vector (right barplot)
# --------------------------
if (!is.null(score_file)) {
  score_df <- read.table(score_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("id", "score") %in% colnames(score_df))) {
    stop("score.tsv must contain columns: id, score")
  }
  missing_s <- setdiff(rownames(sim), score_df$id)
  if (length(missing_s) > 0) {
    stop("These matrix row IDs are missing in score.tsv (id column):\n", paste(missing_s, collapse = ", "))
  }
  score_vector <- setNames(as.numeric(score_df$score[match(rownames(sim), score_df$id)]), rownames(sim))
} else {
  # Default: row mean similarity
  score_vector <- rowMeans(sim, na.rm = TRUE)
  names(score_vector) <- rownames(sim)
}

# --------------------------
# 4) Type colors (edit if needed)
# --------------------------
type_levels <- unique(type_vector)  # keep appearance order by default

# Default palette aligned to your categories (can be extended)
type_colors <- c(
  "Ectoine"        = "#0571b0",
  "Hybrid"         = "#8c6bb1",
  "LAP"            = "#41ab5d",
  "Lassopeptide"   = "#238b45",
  "NI-siderophore" = "#f16913",
  "NRPS"           = "#1d91c0",
  "T1PKS"          = "#de2d26",
  "Terpene"        = "#fd8d3c"
)

# Ensure all types have a color (fallback to grey)
missing_types <- setdiff(type_levels, names(type_colors))
if (length(missing_types) > 0) {
  warning("No color provided for types: ", paste(missing_types, collapse = ", "),
          ". They will be colored as grey.")
  for (tp in missing_types) type_colors[tp] <- "grey70"
}

# --------------------------
# 5) Color function for similarity
# --------------------------
# If user provides col_min/mid/max, use them; otherwise infer from data range.
rng <- range(sim, na.rm = TRUE)

if (is.null(col_min) || is.null(col_mid) || is.null(col_max)) {
  # Reasonable default: treat as [0,1] if within that; else use data range.
  if (rng[1] >= 0 && rng[2] <= 1) {
    col_min <- 0; col_mid <- 0.5; col_max <- 1
  } else {
    col_min <- rng[1]
    col_mid <- mean(rng)
    col_max <- rng[2]
  }
} else {
  col_min <- as.numeric(col_min)
  col_mid <- as.numeric(col_mid)
  col_max <- as.numeric(col_max)
}

col_fun <- circlize::colorRamp2(
  c(col_min, col_mid, col_max),
  c("white", "skyblue", "navy")
)

# --------------------------
# 6) Ordering: split by type, keep within-type order as matrix order
# --------------------------
row_split <- factor(type_vector, levels = unique(type_vector))

# Left annotation: Type
ha_left <- rowAnnotation(
  Type = type_vector,
  col = list(Type = type_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

# Right annotation: score bar + compound text
ha_right <- rowAnnotation(
  score = anno_barplot(score_vector, border = FALSE, gp = gpar(fill = "red"), width = unit(2, "cm")),
  compound = anno_text(rownames(sim), just = "left", gp = gpar(fontsize = 8)),
  annotation_name_side = "top"
)

# --------------------------
# 7) Plot
# --------------------------
ht <- Heatmap(
  sim,
  name = title,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  row_title_side = "left",
  left_annotation = ha_left,
  right_annotation = ha_right,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 45
)

pdf_file <- file.path(outdir, "Fig5A_similarity_heatmap.pdf")
png_file <- file.path(outdir, "Fig5A_similarity_heatmap.png")

pdf(pdf_file, width = width, height = height)
draw(ht, merge_legend = TRUE)
dev.off()

png(png_file, width = width * 300, height = height * 300, res = 300)
draw(ht, merge_legend = TRUE)
dev.off()

# Save ordered metadata for record
meta_out <- data.frame(
  id = rownames(sim),
  type = type_vector,
  score = score_vector,
  stringsAsFactors = FALSE
)
write.table(meta_out, file.path(outdir, "Fig5A_row_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done:\n", pdf_file, "\n", png_file, "\n", sep = "")
