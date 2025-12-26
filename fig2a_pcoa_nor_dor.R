#!/usr/bin/env Rscript

# Figure 2A: PCoA (NOR vs DOR)
# Input CSV must contain columns: PC1, PC2, group (group in {NOR, DOR})
# Outputs: Fig2A_pcoa.pdf/png and (optional) MANOVA test result.

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

input     <- args[["input"]]
outdir    <- args[["outdir"]]
pc1_label <- args[["pc1_label"]]
pc2_label <- args[["pc2_label"]]
title     <- args[["title"]]
envelope  <- args[["envelope"]]   # hull | ellipse | none
do_manova <- args[["do_manova"]]  # TRUE | FALSE

if (is.null(input) || is.null(outdir)) {
  cat("Usage:\n",
      "  Rscript scripts/fig2a_pcoa_nor_dor.R --input <Data_for_PcoA.csv> --outdir <outdir>\n",
      "  [--pc1_label \"PC1 (..%)\"] [--pc2_label \"PC2 (..%)\"] [--title \"PCoA: NOR vs DOR\"]\n",
      "  [--envelope hull|ellipse|none] [--do_manova TRUE|FALSE]\n", sep = "")
  quit(status = 1)
}
if (!file.exists(input)) stop("Input not found: ", input)

if (is.null(pc1_label)) pc1_label <- "PC1"
if (is.null(pc2_label)) pc2_label <- "PC2"
if (is.null(title))     title     <- "PCoA: NOR vs DOR"
if (is.null(envelope))  envelope  <- "hull"
if (is.null(do_manova)) do_manova <- "TRUE"

need(c("ggplot2", "dplyr"))
library(ggplot2)
library(dplyr)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input, header = TRUE, check.names = FALSE)
req <- c("PC1", "PC2", "group")
if (!all(req %in% colnames(df))) {
  stop("Input must contain columns: PC1, PC2, group")
}

df_sub <- df %>% filter(group %in% c("NOR", "DOR"))
df_sub$group <- factor(df_sub$group, levels = c("NOR", "DOR"))

# Envelope data
hulls <- NULL
if (envelope == "hull") {
  hulls <- df_sub %>%
    group_by(group) %>%
    slice(chull(PC1, PC2))
}

p <- ggplot(df_sub, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 3) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(x = pc1_label, y = pc2_label, title = title)

if (envelope == "hull") {
  p <- p + geom_polygon(data = hulls, alpha = 0.2, aes(fill = group), show.legend = FALSE)
} else if (envelope == "ellipse") {
  p <- p + stat_ellipse(level = 0.95, type = "t", alpha = 0.2, geom = "polygon", show.legend = FALSE)
}

pdf_file <- file.path(outdir, "Fig2A_pcoa.pdf")
png_file <- file.path(outdir, "Fig2A_pcoa.png")
ggsave(pdf_file, p, width = 7.2, height = 6.0, limitsize = FALSE)
ggsave(png_file, p, width = 7.2, height = 6.0, dpi = 300, limitsize = FALSE)

# Optional MANOVA test (Wilks)
if (toupper(do_manova) == "TRUE") {
  fit <- manova(cbind(PC1, PC2) ~ group, data = df_sub)
  res <- summary(fit, test = "Wilks")
  out_txt <- file.path(outdir, "Fig2A_manova_wilks.txt")
  capture.output(res, file = out_txt)
}

cat("Done:\n", pdf_file, "\n", png_file, "\n", sep = "")
