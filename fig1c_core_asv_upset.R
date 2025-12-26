#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Fig1C: Core ASV (prevalence >= threshold) + UpSet/Venn + shared ASVs export
#
# Input (per group):
#   <group>.txt  (TSV with header, rownames=OTU IDs, columns=samples)
# Core OTU definition:
#   ASV is "core" in a group if it appears (value > 0) in >= threshold * N samples.
#
# Outputs:
#   - core_tables/<group>_core_<pct>.tsv
#   - core_ids/<group>_core_ids_<pct>.txt
#   - shared_core_otus_<pct>.txt
#   - shared_core_otu_table_<pct>.tsv (from --shared_table_source)
#   - upset_core_otus_<pct>.pdf/png
#   - venn_core_otus_<pct>.png (optional; only 2-5 sets)
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
      "  install.packages(c('UpSetR','VennDiagram'))\n"
    )
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

indir  <- args[["indir"]]
groups <- args[["groups"]]          # comma-separated
outdir <- args[["outdir"]]
threshold <- args[["threshold"]]    # e.g. 0.90
suffix <- args[["suffix"]]          # default ".txt"
make_upset <- args[["make_upset"]]  # TRUE/FALSE
make_venn  <- args[["make_venn"]]   # TRUE/FALSE
venn_groups <- args[["venn_groups"]]# comma-separated (2-5 groups)
shared_table_source <- args[["shared_table_source"]]  # group name to export shared table from

if (is.null(indir) || is.null(groups) || is.null(outdir)) {
  cat(
    "Usage:\n",
    "  Rscript scripts/fig1c_core_otu_upset.R \\\n",
    "    --indir <dir_with_group_tables> \\\n",
    "    --groups \"SXBS,SXRS,...\" \\\n",
    "    --threshold 0.90 \\\n",
    "    --outdir <outdir> \\\n",
    "    [--suffix .txt] [--make_upset TRUE] [--make_venn FALSE] \\\n",
    "    [--venn_groups \"SXRT,YYRT,YZRT\"] [--shared_table_source SXRT]\n",
    sep = ""
  )
  quit(status = 1)
}

if (!dir.exists(indir)) stop("indir not found: ", indir)

if (is.null(threshold)) threshold <- "0.90"
threshold <- as.numeric(threshold)
if (threshold <= 0 || threshold > 1) stop("--threshold must be in (0, 1].")

if (is.null(suffix)) suffix <- ".txt"
if (is.null(make_upset)) make_upset <- "TRUE"
if (is.null(make_venn))  make_venn  <- "FALSE"

groups_vec <- unlist(strsplit(groups, ","))
if (length(groups_vec) < 2) stop("Need at least 2 groups.")

if (is.null(shared_table_source)) shared_table_source <- groups_vec[1]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
core_tbl_dir <- file.path(outdir, "core_tables")
core_id_dir  <- file.path(outdir, "core_ids")
dir.create(core_tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(core_id_dir,  recursive = TRUE, showWarnings = FALSE)

pct_tag <- paste0(round(threshold * 100), "percent")

# --------------------------
# 1) Compute core ASVs per group
# --------------------------
otu_list <- list()
core_tables <- list()

for (g in groups_vec) {
  f <- file.path(indir, paste0(g, suffix))
  if (!file.exists(f)) stop("Group file not found: ", f)

  otu <- read.table(f, header = TRUE, row.names = 1, sep = "\t",
                    check.names = FALSE, quote = "", comment.char = "")

  if (nrow(otu) == 0 || ncol(otu) == 0) stop("Empty table: ", f)

  presence <- rowSums(otu > 0)
  thr_n <- threshold * ncol(otu)

  core_otu <- otu[presence >= thr_n, , drop = FALSE]

  # Save core table
  out_tbl <- file.path(core_tbl_dir, paste0(g, "_core_", pct_tag, ".tsv"))
  write.table(core_otu, file = out_tbl, sep = "\t", quote = FALSE, col.names = NA)

  # Save core OTU ids
  out_ids <- file.path(core_id_dir, paste0(g, "_core_ids_", pct_tag, ".txt"))
  write.table(rownames(core_otu), file = out_ids,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  otu_list[[g]] <- rownames(core_otu)
  core_tables[[g]] <- core_otu

  cat(sprintf("[Core] %s: %d core OTUs (threshold=%.2f of %d samples)\n",
              g, nrow(core_otu), threshold, ncol(otu)))
}

# --------------------------
# 2) Shared ASVs across all groups
# --------------------------
shared_otus <- Reduce(intersect, otu_list)
shared_out <- file.path(outdir, paste0("shared_core_otus_", pct_tag, ".txt"))
write.table(shared_otus, file = shared_out,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Shared core OTUs across ALL groups: ", length(shared_otus), "\n", sep = "")

# Export a shared OTU table (from a chosen group's core table or raw table)
if (!shared_table_source %in% names(core_tables)) {
  stop("--shared_table_source must be one of: ", paste(names(core_tables), collapse = ", "))
}
src_tbl <- core_tables[[shared_table_source]]
shared_tbl <- src_tbl[rownames(src_tbl) %in% shared_otus, , drop = FALSE]
shared_tbl_out <- file.path(outdir, paste0("shared_core_otu_table_", pct_tag, "_from_", shared_table_source, ".tsv"))
write.table(shared_tbl, file = shared_tbl_out, sep = "\t", quote = FALSE, col.names = NA)

# --------------------------
# 3) UpSet plot (recommended for 6+ sets)
# --------------------------
if (toupper(make_upset) == "TRUE") {
  need(c("UpSetR"))

  all_otus <- sort(unique(unlist(otu_list)))
  binary_matrix <- sapply(otu_list, function(x) all_otus %in% x)
  binary_df <- as.data.frame(binary_matrix * 1)
  # Keep set order
  binary_df <- binary_df[, groups_vec, drop = FALSE]

  upset_pdf <- file.path(outdir, paste0("upset_core_otus_", pct_tag, ".pdf"))
  upset_png <- file.path(outdir, paste0("upset_core_otus_", pct_tag, ".png"))

  pdf(upset_pdf, width = 10, height = 6)
  UpSetR::upset(binary_df,
                sets = groups_vec,
                order.by = "freq",
                keep.order = TRUE,
                nsets = length(groups_vec))
  dev.off()

  png(upset_png, width = 3000, height = 1800, res = 300)
  UpSetR::upset(binary_df,
                sets = groups_vec,
                order.by = "freq",
                keep.order = TRUE,
                nsets = length(groups_vec))
  dev.off()

  cat("UpSet saved:\n  - ", upset_pdf, "\n  - ", upset_png, "\n", sep = "")
}

# --------------------------
# 4) Venn plot (only 2–5 sets)
# --------------------------
if (toupper(make_venn) == "TRUE") {
  need(c("VennDiagram", "grid"))

  if (is.null(venn_groups) || nchar(venn_groups) == 0) {
    # default: first 3 groups
    venn_vec <- groups_vec[1:min(3, length(groups_vec))]
  } else {
    venn_vec <- unlist(strsplit(venn_groups, ","))
  }

  if (length(venn_vec) < 2 || length(venn_vec) > 5) {
    stop("--venn_groups must contain 2 to 5 groups.")
  }
  if (!all(venn_vec %in% names(otu_list))) {
    stop("Some --venn_groups are not in --groups: ", paste(setdiff(venn_vec, names(otu_list)), collapse = ", "))
  }

  venn_list <- otu_list[venn_vec]

  venn_png <- file.path(outdir, paste0("venn_core_otus_", pct_tag, "_", paste(venn_vec, collapse = "_"), ".png"))

  # VennDiagram draws via grid; use png device
  png(venn_png, width = 2400, height = 2000, res = 300)
  vp <- VennDiagram::venn.diagram(
    x = venn_list,
    category.names = venn_vec,
    filename = NULL,
    imagetype = "png",
    cex = 1.0,
    cat.cex = 1.0,
    fill = grDevices::rainbow(length(venn_vec)),
    alpha = 0.45,
    main = paste0("Core OTUs Venn (", pct_tag, ")")
  )
  grid::grid.newpage()
  grid::grid.draw(vp)
  dev.off()

  cat("Venn saved:\n  - ", venn_png, "\n", sep = "")
}

cat("Done. Outputs under: ", outdir, "\n", sep = "")
