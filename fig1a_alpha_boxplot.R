#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Fig1A: Alpha diversity boxplot + ANOVA/TukeyHSD + LSD letters
# Refactoring of Zhang2019NBT-style workflow:
#   - read alpha table (rows=samples)
#   - read design (rows=samples) + group column
#   - subset groups (optional) + enforce group order
#   - per index: aov -> TukeyHSD table -> LSD.test letters
#   - plot boxplot + jitter + letters
#   - export stats tables + figures
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
      "  install.packages(c('ggplot2','dplyr','readr','tibble','agricolae'))\n"
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

alpha_file <- args[["alpha"]]
design_file <- args[["design"]]
group_col <- args[["group_col"]]
groups_str <- args[["groups"]]     # comma-separated
indices_str <- args[["indices"]]   # comma-separated
outdir <- args[["outdir"]]
p_adjust <- args[["p_adjust"]]     # none | fdr | bonferroni etc (agricolae supports several)
width <- args[["width"]]
height <- args[["height"]]
seed <- args[["seed"]]

if (is.null(alpha_file) || is.null(design_file) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/fig1a_alpha_boxplot.R \\\n",
    "    --alpha <alpha4.txt> --design <design3.txt> --outdir <outdir> \\\n",
    "    [--group_col groupID] \\\n",
    "    [--groups \"SXBS,SXRS,...\"] \\\n",
    "    [--indices \"Shannon,Simpson\"] \\\n",
    "    [--p_adjust none] [--width 5] [--height 3] [--seed 123]\n",
    sep = ""
  )
  quit(status = 1)
}

if (!file.exists(alpha_file)) stop("alpha file not found: ", alpha_file)
if (!file.exists(design_file)) stop("design file not found: ", design_file)

if (is.null(group_col)) group_col <- "groupID"
if (is.null(p_adjust)) p_adjust <- "none"
if (is.null(width)) width <- "5"
if (is.null(height)) height <- "3"
if (is.null(seed)) seed <- "123"

width <- as.numeric(width)
height <- as.numeric(height)
seed <- as.integer(seed)
set.seed(seed)

need(c("ggplot2", "dplyr", "tibble", "agricolae"))
library(ggplot2)
library(dplyr)
library(tibble)
library(agricolae)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------
# 1) Read inputs
# --------------------------
alpha <- read.table(alpha_file, header = TRUE, row.names = 1, sep = "\t",
                    check.names = FALSE, quote = "", comment.char = "")

design <- read.table(design_file, header = TRUE, row.names = 1, sep = "\t",
                     check.names = FALSE, quote = "", comment.char = "")

if (!group_col %in% colnames(design)) {
  stop("design does not contain group column: ", group_col)
}

design$group <- design[[group_col]]

# Optional group selection + order
if (!is.null(groups_str) && nchar(groups_str) > 0) {
  groups <- unlist(strsplit(groups_str, ","))
  design <- subset(design, group %in% groups)
  design$group <- factor(design$group, levels = groups)
} else {
  # Keep as factor but preserve appearance order
  design$group <- factor(design$group)
}

# Cross-filter by sample IDs
common <- intersect(rownames(alpha), rownames(design))
if (length(common) < 3) stop("Too few overlapping samples between alpha and design.")
alpha <- alpha[common, , drop = FALSE]
design <- design[common, , drop = FALSE]

index <- cbind(alpha, design)

# Determine indices to plot
if (!is.null(indices_str) && nchar(indices_str) > 0) {
  indices <- unlist(strsplit(indices_str, ","))
  missing_idx <- indices[!indices %in% colnames(index)]
  if (length(missing_idx) > 0) {
    stop("These indices are not found in alpha table: ", paste(missing_idx, collapse = ", "))
  }
} else {
  # fallback: plot all numeric columns from alpha
  alpha_cols <- colnames(alpha)
  is_num <- vapply(alpha[alpha_cols], is.numeric, FUN.VALUE = logical(1))
  indices <- alpha_cols[is_num]
  if (length(indices) == 0) stop("No numeric alpha indices found in alpha table.")
}

# Default palette for the 9 groups in your design (customize if needed)
default_palette <- c(
  "SXBS"="#1b9e77","SXRS"="#d95f02","SXRT"="#7570b3",
  "YYBS"="#e7298a","YYRS"="#66a61e","YYRT"="#e6ab02",
  "YZBS"="#a6761d","YZRS"="#666666","YZRT"="#ff7f00"
)

# If groups are different, ggplot will fall back unless fully covered
use_manual_colors <- all(levels(index$group) %in% names(default_palette))

main_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

# --------------------------
# 2) Loop each index
# --------------------------
for (m in indices) {
  cat("Processing index: ", m, "\n", sep = "")

  # ANOVA
  model <- aov(index[[m]] ~ group, data = index)

  # TukeyHSD
  tuk <- TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  tuk_tab <- as.data.frame(tuk$group) %>%
    rownames_to_column("comparison")

  tuk_file <- file.path(outdir, paste0("alpha_", m, "_tukey.tsv"))
  write.table(tuk_tab, tuk_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # LSD letters (agricolae)
  lsd <- agricolae::LSD.test(model, "group", p.adj = p_adjust)
  letters_tab <- lsd$groups %>%
    rownames_to_column("group") %>%
    as.data.frame()

  letters_file <- file.path(outdir, paste0("alpha_", m, "_lsd_letters.tsv"))
  write.table(letters_tab, letters_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Join letters to each point for plotting
  plot_df <- index %>%
    mutate(group = factor(group, levels = levels(design$group))) %>%
    left_join(letters_tab[, c("group", "groups")], by = "group") %>%
    rename(stat = groups)

  # One label per group (avoid repeating letters for every point)
  y_max <- max(plot_df[[m]], na.rm = TRUE)
  y_min <- min(plot_df[[m]], na.rm = TRUE)
  offset <- (y_max - y_min) * 0.05

  label_df <- plot_df %>%
    group_by(group) %>%
    summarise(
      y = max(.data[[m]], na.rm = TRUE) + offset,
      stat = dplyr::first(stat),
      .groups = "drop"
    )

  p <- ggplot(plot_df, aes(x = group, y = .data[[m]], color = group)) +
    geom_boxplot(alpha = 1, outlier.size = 0, size = 0.7, width = 0.5, fill = "transparent") +
    geom_jitter(position = position_jitter(width = 0.17), size = 1, alpha = 0.7) +
    geom_text(data = label_df, aes(x = group, y = y, label = stat), inherit.aes = FALSE) +
    labs(x = "Groups", y = m) +
    main_theme

  if (length(levels(plot_df$group)) > 3) {
    p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  }

  if (use_manual_colors) {
    p <- p + scale_color_manual(values = default_palette)
  }

  pdf_file <- file.path(outdir, paste0("alpha_", m, ".pdf"))
  png_file <- file.path(outdir, paste0("alpha_", m, ".png"))
  ggsave(pdf_file, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

cat("Done. Outputs under: ", outdir, "\n", sep = "")
