#!/usr/bin/env Rscript
# Regression check against the figures omakase was built to reproduce.
#
# Renders a start-site cohort through omakase and writes the pages next to the
# originals, so the two can be compared one by one. Nothing here is shipped
# with the package, and no data is committed: the source tables stay on the
# machine that produced them, and this script does nothing unless they are
# present. Point it somewhere with:
#
#   OMAKASE_REFERENCE_SRC=/path/to/tables Rscript data-raw/verify_reference.R [gene ...]

devtools::load_all(".", quiet = TRUE)

SRC <- Sys.getenv("OMAKASE_REFERENCE_SRC", unset = "")
# A stable, git-ignored directory rather than tempdir(), so the pages are still
# there to look at after the session that made them has exited.
OUT <- Sys.getenv("OMAKASE_CHECK_DIR", "reference-check")

if (!nzchar(SRC) || !dir.exists(SRC)) {
  message("Set OMAKASE_REFERENCE_SRC to a directory of source tables to run ",
          "this check. Nothing to do.")
  quit(status = 0)
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

genes <- commandArgs(trailingOnly = TRUE)

# The originals are keyed by gene name and use project-specific column names;
# `rename` maps them onto the sashimi_data contract without touching the files.
sd <- sashimi_from_tables(
  loci      = file.path(SRC, "*sashimi_genes.tsv"),
  tracks    = file.path(SRC, "*sashimi_tracks.parquet"),
  junctions = file.path(SRC, "*sashimi_arcs.tsv"),
  models    = file.path(SRC, "*sashimi_models.tsv"),
  psi       = file.path(SRC, "*sashimi_psi.tsv"),
  features  = file.path(SRC, "*sashimi_repeats.tsv"),
  rename = list(
    loci      = c(main_apex = "main_apex", alt_apex = "atss_apex"),
    tracks    = c(group = "stage", value = "rpm"),
    junctions = c(group = "stage", x0 = "apex", x1 = "anchor",
                  count = "activity"),
    models    = c(tx_id = "iso_id", start = "exon_start", end = "exon_end"),
    psi       = c(group = "stage"),
    features  = c(start = "rep_start", end = "rep_end", name = "rep_name",
                  class = "rep_class")
  ),
  locus_col = "gene_name"
)

# The originals draw their conditions in biological order, which is the order
# the tables list them in, not alphabetical order. sashimi_from_tables already
# preserves first-appearance order, so nothing more is needed unless
# OMAKASE_REFERENCE_GROUPS overrides it.
grp <- Sys.getenv("OMAKASE_REFERENCE_GROUPS", unset = "")
if (nzchar(grp)) sd$meta$group_order <- trimws(strsplit(grp, ",")[[1]])

# With no genes named, check the first few the cohort happens to contain.
if (!length(genes)) genes <- utils::head(sd$loci$gene_name, 4)

cat("loaded:\n")
print(sd)

for (g in genes) {
  if (!g %in% sd$loci$gene_name) {
    message("  skipping ", g, ": not in the cohort")
    next
  }
  p <- plot_sashimi(sd, preset = "tss", locus = g)
  f <- file.path(OUT, paste0(g, ".omakase.pdf"))
  save_sashimi(p, f, width = 5.9)
  ref <- file.path(SRC, "sashimi", paste0(g, ".pdf"))
  if (file.exists(ref)) file.copy(ref, file.path(OUT, paste0(g, ".reference.pdf")),
                                  overwrite = TRUE)
  cat("  ", g, "\n")
}

cat("\nwrote to ", OUT, "\n", sep = "")
cat("compare <gene>.omakase.pdf against <gene>.reference.pdf\n")
