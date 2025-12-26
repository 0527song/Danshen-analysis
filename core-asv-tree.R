#!/usr/bin/env Rscript

# ------------------------------------------------------------
# core_otu_tree.R
# Build a phylogenetic tree from representative OTU/ASV fasta,
# compute grouped mean relative abundance, and plot tree + heatmap.
#
# Refactored from the user's original script (core-otu_tree.R):
# - Tree inference: DECIPHER alignment + phangorn NJ + ML (GTR)
# - Group mean relative abundance from OTU table
# - Plot with ggtree + ggtreeExtra::gheatmap
# - Export shared OTUs across all groups
# ------------------------------------------------------------

# --------------------------
# 0) Helpers: arg parsing
# --------------------------
parse_args <- function(args) {
  # Accept: --key value  OR  --key=value
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
      "Please install them before running this script.\n",
      "Tip: use renv or BiocManager for Bioconductor packages."
    )
  }
}

# --------------------------
# 1) Read args & validate
# --------------------------
args <- parse_args(commandArgs(trailingOnly = TRUE))

fasta    <- args[["fasta"]]
otu_file <- args[["otu"]]
outdir   <- args[["outdir"]]
metadata <- args[["metadata"]]  # optional

if (is.null(fasta) || is.null(otu_file) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/core_otu_tree.R --fasta <selected_otus.fasta> --otu <otu_table.txt> --outdir <results_dir> [--metadata <metadata.tsv>]\n\n",
    "metadata.tsv format (optional): two columns with header: sample<TAB>group\n",
    "Example:\n",
    "sample\tgroup\n",
    "SXRS101\tSXRS\n",
    "YYRS103\tYYRS\n",
    "...\n"
  )
  quit(status = 1)
}

if (!file.exists(fasta))    stop("FASTA not found: ", fasta)
if (!file.exists(otu_file)) stop("OTU table not found: ", otu_file)
if (!is.null(metadata) && !file.exists(metadata)) stop("Metadata not found: ", metadata)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 2) Load packages
# --------------------------
need(c("ape", "phangorn", "ggplot2"))
need(c("ggtree", "ggtreeExtra"))
need(c("dplyr", "tidyr"))

# Bioconductor packages (often)
need(c("Biostrings", "DECIPHER"))

# Attach (keep minimal)
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(tidyr)

# Use explicit namespaces for Bioc pkgs to reduce conflicts
# Biostrings::readDNAStringSet
# DECIPHER::AlignSeqs

# --------------------------
# 3) Build phylogenetic tree
# --------------------------
cat("[1/5] Building phylogenetic tree from FASTA...\n")

seqs <- Biostrings::readDNAStringSet(fasta)
if (length(seqs) < 3) stop("FASTA must contain >= 3 sequences to build a tree.")

aln <- DECIPHER::AlignSeqs(seqs)
phang_align <- phangorn::as.phyDat(aln, type = "DNA")
dm <- phangorn::dist.ml(phang_align)

tree_nj <- ape::NJ(dm)

# ML optimization (GTR)
fit <- phangorn::pml(tree_nj, phang_align)
fit_gtr <- phangorn::optim.pml(
  fit,
  model = "GTR",
  rearrangement = "stochastic",
  control = phangorn::pml.control(trace = 0)
)

tree <- fit_gtr$tree

tree_file <- file.path(outdir, "tree.nwk")
ape::write.tree(tree, file = tree_file)
cat("  Tree saved: ", tree_file, "\n", sep = "")

# --------------------------
# 4) Read OTU table & compute group mean relative abundance
# --------------------------
cat("[2/5] Reading OTU table and computing grouped mean relative abundance...\n")

otu_table <- read.table(
  otu_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (nrow(otu_table) == 0 || ncol(otu_table) == 0) stop("OTU table looks empty: ", otu_file)

# Convert to relative abundance per sample (column-wise)
col_sums <- colSums(otu_table)
if (any(col_sums == 0)) stop("Some samples have zero total counts in OTU table (cannot compute relative abundance).")
otu_rel <- sweep(otu_table, 2, col_sums, "/")

# Grouping: prefer metadata.tsv; otherwise fallback to built-in groups (as in original script)
group_info <- NULL

if (!is.null(metadata)) {
  meta <- read.table(metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  req_cols <- c("sample", "group")
  if (!all(req_cols %in% colnames(meta))) {
    stop("metadata.tsv must contain columns: sample, group")
  }
  group_info <- split(meta$sample, meta$group)
} else {
  # Fallback: your original hard-coded sample lists (edit as needed for your data)
  group_info <- list(
    SXRS = c("SXRS101","SXRS103","SXRS105","SXRS201","SXRS202","SXRS203","SXRS204","SXRS205",
             "SXRS301","SXRS302","SXRS303","SXRS304","SXRS305","SX2RS104","SX2RS105","SX2RS201",
             "SX2RS203","SX2RS204","SX2RS205","SX2RS301","SX2RS302","SX2RS303","SX2RS304","SX2RS305"),
    YYRS = c("YYRS103","YYRS201","YYRS202","YYRS203","YYRS204","YYRS205","YYRS301","YYRS302","YYRS303","YYRS305"),
    YZRS = c("YZRS101","YZRS102","YZRS103","YZRS104","YZRS105","YZRS202","YZRS203","YZRS204","YZRS205",
             "YZRS301","YZRS302","YZRS303","YZRS304","YZRS305")
  )
}

# Compute mean relative abundance per group
group_avg <- sapply(names(group_info), function(g) {
  samples <- intersect(group_info[[g]], colnames(otu_rel))
  if (length(samples) == 0) stop("No sample columns matched for group: ", g, "\nCheck your sample names or provide --metadata.")
  rowMeans(otu_rel[, samples, drop = FALSE])
})

# Ensure OTU IDs match tree tips
tips <- tree$tip.label
common <- intersect(tips, rownames(group_avg))
if (length(common) == 0) {
  stop(
    "No overlap between tree tip labels and OTU table rownames.\n",
    "Tip: ensure FASTA IDs == OTU table rownames."
  )
}

# Subset to common OTUs and order by tree tips
group_avg <- group_avg[common, , drop = FALSE]
group_avg <- group_avg[tips[tips %in% rownames(group_avg)], , drop = FALSE]

# Save group average table
avg_file <- file.path(outdir, "otu_group_avg.tsv")
write.table(group_avg, avg_file, sep = "\t", quote = FALSE, col.names = NA)
cat("  Group mean relative abundance saved: ", avg_file, "\n", sep = "")

# Export shared OTU IDs present in all groups (>0)
present_in_all <- rownames(group_avg)[rowSums(group_avg > 0) == ncol(group_avg)]
shared_file <- file.path(outdir, "shared_otu_ids.txt")
write.table(present_in_all, shared_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  Shared OTU IDs saved: ", shared_file, "\n", sep = "")

# --------------------------
# 5) Plot: tree + heatmap
# --------------------------
cat("[3/5] Plotting tree + heatmap...\n")

p_tree <- ggtree(tree, layout = "rectangular") +
  geom_tiplab(align = TRUE, linesize = 0.5, size = 3)

# gheatmap expects a matrix/data.frame with rownames = tip labels
heat_df <- as.data.frame(group_avg)
heat_df <- heat_df[tree$tip.label[tree$tip.label %in% rownames(heat_df)], , drop = FALSE]

p_final <- ggtreeExtra::gheatmap(
  p_tree,
  heat_df,
  offset = 0.02,
  width = 0.4,
  colnames_position = "top",
  colnames_angle = 45,
  colnames_offset_y = 0.5,
  font.size = 3
) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.position = "right")

pdf_file <- file.path(outdir, "tree_heatmap.pdf")
png_file <- file.path(outdir, "tree_heatmap.png")

ggsave(pdf_file, p_final, width = 10, height = 8, limitsize = FALSE)
ggsave(png_file, p_final, width = 10, height = 8, dpi = 300, limitsize = FALSE)

cat("  Plot saved:\n")
cat("   - ", pdf_file, "\n", sep = "")
cat("   - ", png_file, "\n", sep = "")

cat("[4/5] Done.\n")
