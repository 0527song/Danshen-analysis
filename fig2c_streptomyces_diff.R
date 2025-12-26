#!/usr/bin/env Rscript

# Figure 2C: Streptomyces differential abundance (NOR vs DOR)
# Input: CSV/TSV with an ID column + replicate columns like NOR1.. / DOR1..
# Outputs: stats table + volcano + jitter + bubble plot

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
  if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

input      <- args[["input"]]
outdir     <- args[["outdir"]]
id_col     <- args[["id_col"]]
nor_prefix <- args[["nor_prefix"]]
dor_prefix <- args[["dor_prefix"]]
p_cutoff   <- args[["p_cutoff"]]
lfc_cutoff <- args[["lfc_cutoff"]]
pseudocnt  <- args[["pseudocount"]]
highlight  <- args[["highlight"]]
make_bub   <- args[["make_bubble"]]

if (is.null(input) || is.null(outdir)) {
  cat("Usage:\n",
      "  Rscript scripts/fig2c_streptomyces_diff.R \\\n",
      "    --input <abundance_table.csv/tsv> --outdir <outdir> \\\n",
      "    [--id_col ID] [--nor_prefix NOR] [--dor_prefix DOR] \\\n",
      "    [--p_cutoff 0.05] [--lfc_cutoff 1] [--pseudocount 1e-6] \\\n",
      "    [--highlight Streptomyces_aurantiacus] [--make_bubble TRUE|FALSE]\n", sep = "")
  quit(status = 1)
}
if (!file.exists(input)) stop("Input not found: ", input)

if (is.null(id_col))     id_col <- "ID"
if (is.null(nor_prefix)) nor_prefix <- "NOR"
if (is.null(dor_prefix)) dor_prefix <- "DOR"
if (is.null(p_cutoff))   p_cutoff <- "0.05"
if (is.null(lfc_cutoff)) lfc_cutoff <- "1"
if (is.null(pseudocnt))  pseudocnt <- "1e-6"
if (is.null(highlight))  highlight <- "Streptomyces_aurantiacus"
if (is.null(make_bub))   make_bub <- "TRUE"

p_cutoff   <- as.numeric(p_cutoff)
lfc_cutoff <- as.numeric(lfc_cutoff)
pseudocnt  <- as.numeric(pseudocnt)

need(c("dplyr", "tidyr", "ggplot2", "ggrepel"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read CSV or TSV by extension
is_tsv <- grepl("\\.tsv$|\\.txt$", tolower(input))
df <- if (is_tsv) {
  read.table(input, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
} else {
  read.csv(input, header = TRUE, check.names = FALSE)
}

if (!id_col %in% colnames(df)) stop("ID column not found: ", id_col)

# Identify NOR/DOR columns
nor_cols <- grep(paste0("^", nor_prefix), colnames(df), value = TRUE)
dor_cols <- grep(paste0("^", dor_prefix), colnames(df), value = TRUE)
if (length(nor_cols) < 2 || length(dor_cols) < 2) {
  stop("Need >=2 replicates for each group. Found NOR cols: ",
       paste(nor_cols, collapse = ", "),
       " ; DOR cols: ",
       paste(dor_cols, collapse = ", "))
}

# Ensure numeric
df[, c(nor_cols, dor_cols)] <- lapply(df[, c(nor_cols, dor_cols), drop = FALSE], as.numeric)

# Stats
df_stat <- df %>%
  rowwise() %>%
  mutate(
    NOR_mean = mean(c_across(all_of(nor_cols)), na.rm = TRUE),
    DOR_mean = mean(c_across(all_of(dor_cols)), na.rm = TRUE),
    log2FC   = log2((DOR_mean + pseudocnt) / (NOR_mean + pseudocnt)),
    p_value  = tryCatch(
      t.test(c_across(all_of(dor_cols)), c_across(all_of(nor_cols)))$p.value,
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    negLog10P = -log10(p_value),
    Status = case_when(
      !is.na(p_value) & p_value < p_cutoff & log2FC >=  lfc_cutoff ~ "enriched",
      !is.na(p_value) & p_value < p_cutoff & log2FC <= -lfc_cutoff ~ "depleted",
      TRUE ~ "no_sig"
    )
  )

# Save table
tab_file <- file.path(outdir, "Fig2C_streptomyces_stats.tsv")
write.table(df_stat, tab_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano
volcano <- ggplot(df_stat, aes(x = log2FC, y = negLog10P, color = Status)) +
  geom_point(alpha = 0.85, size = 2) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  scale_color_manual(values = c(enriched = "red", depleted = "blue", no_sig = "grey60")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Streptomyces differential abundance (DOR vs NOR)",
    x = "log2(Fold Change)",
    y = "-log10(P-value)",
    color = NULL
  )

# Highlight one label (optional)
df_stat$label <- ifelse(df_stat[[id_col]] == highlight, df_stat[[id_col]], NA)

volcano <- volcano +
  geom_text_repel(aes(label = label),
                  size = 4,
                  color = "black",
                  segment.color = "black",
                  segment.linetype = "dashed",
                  max.overlaps = Inf,
                  na.rm = TRUE)

vol_pdf <- file.path(outdir, "Fig2C_volcano.pdf")
vol_png <- file.path(outdir, "Fig2C_volcano.png")
ggsave(vol_pdf, volcano, width = 7.2, height = 6.2, limitsize = FALSE)
ggsave(vol_png, volcano, width = 7.2, height = 6.2, dpi = 300, limitsize = FALSE)

# Jitter plot (single-column comparison)
jitter_df <- df_stat %>%
  mutate(Comparison = "DOR vs NOR")

jitter <- ggplot(jitter_df, aes(x = Comparison, y = log2FC, color = Status)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -lfc_cutoff, ymax = lfc_cutoff,
           alpha = 0.25) +
  geom_jitter(width = 0.15, alpha = 0.85, size = 2) +
  scale_color_manual(values = c(enriched = "red", depleted = "blue", no_sig = "grey60")) +
  theme_minimal(base_size = 13) +
  labs(title = "Streptomyces (DOR vs NOR)", x = NULL, y = "log2(Fold Change)", color = NULL)

jit_pdf <- file.path(outdir, "Fig2C_jitter.pdf")
jit_png <- file.path(outdir, "Fig2C_jitter.png")
ggsave(jit_pdf, jitter, width = 5.8, height = 5.5, limitsize = FALSE)
ggsave(jit_png, jitter, width = 5.8, height = 5.5, dpi = 300, limitsize = FALSE)

# Bubble plot (abundance by sample)
if (toupper(make_bub) == "TRUE") {
  long <- df %>%
    pivot_longer(cols = all_of(c(nor_cols, dor_cols)),
                 names_to = "Sample", values_to = "Abundance") %>%
    mutate(Group = ifelse(grepl(paste0("^", nor_prefix), Sample), "NOR", "DOR"))

  bubble <- ggplot(long, aes(x = Sample, y = .data[[id_col]])) +
    geom_point(aes(size = Abundance, color = Group), alpha = 0.75) +
    scale_size(range = c(1, 9)) +
    scale_color_manual(values = c(NOR = "#1f77b4", DOR = "#ff7f0e")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Streptomyces abundance (bubble plot)",
         x = "Sample", y = "Species/ID", size = "Abundance", color = "Group")

  bub_pdf <- file.path(outdir, "Fig2C_bubble.pdf")
  bub_png <- file.path(outdir, "Fig2C_bubble.png")
  ggsave(bub_pdf, bubble, width = 10.5, height = 7.5, limitsize = FALSE)
  ggsave(bub_png, bubble, width = 10.5, height = 7.5, dpi = 300, limitsize = FALSE)
}

cat("Done:\n",
    " - ", tab_file, "\n",
    " - ", vol_pdf, "\n",
    " - ", jit_pdf, "\n", sep = "")
