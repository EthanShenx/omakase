#!/usr/bin/env Rscript
# Regenerate every figure shown in README.md, at 600 dpi.
#
# Sources: the six public ENCODE libraries over the human ABI1 locus shipped in
# inst/extdata, the real skipped-exon event derived from those junctions, the
# synthetic consequence table, and a synthetic start-site dataset built here.
# Nothing is hand-drawn or retouched afterwards.
#
#   Rscript data-raw/make_figures.R

devtools::load_all(".", quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

ex <- function(...) file.path("inst", "extdata", ...)
FIG <- file.path("man", "figures")
DPI <- 600
# Every figure is set in Arial. Text geoms do not inherit a theme's family, so
# the package threads `base_family` to each of them; the theme defaults below
# cover the few pieces ggplot2 draws on its own (titles, axis text, legends).
FAM <- "Arial"
theme_set(theme_get() + theme(text = element_text(family = FAM)))
update_geom_defaults("text", list(family = FAM))
update_geom_defaults("label", list(family = FAM))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

png_out <- function(p, name, width, height, saver = ggsave) {
  f <- file.path(FIG, paste0(name, ".png"))
  # ragg writes a smaller file than grDevices at the same resolution, which
  # matters when every figure is 600 dpi.
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL
  args <- list(f, p, width = width, height = height, units = "in", dpi = DPI,
               bg = "white")
  if (!is.null(dev)) args$device <- dev
  do.call(saver, args)
  cat(sprintf("  %-26s %5.1f x %-4.1f in  %6.1f KB\n", name, width, height,
              file.size(f) / 1024))
  invisible(f)
}

REGION <- "chr10:27040584-27048100"
GTF <- ex("annotation.gtf")
BAMS <- ex("samples.tsv")

cat("reading alignments...\n")
sd <- sashimi_from_bam(BAMS, REGION, annotation = GTF, min_count = 10)
sd <- add_models(sd, read_annotation(GTF, gene = "ABI1",
                                     feature = c("exon", "CDS")), max_tx = 4)

cat("drawing sashimi:\n")

# ---- 1. the default figure ------------------------------------------------
png_out(
  plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count"),
  "sashimi-default", 7.2, 4.9
)

# ---- 2. one panel per library, fixed scale --------------------------------
png_out(
  plot_sashimi(sd, base_family = FAM, arc_label_format = "count", fix_y_scale = TRUE,
               arc_width_rule = "log", arc_width = 0.3),
  "sashimi-replicates", 7.2, 7.6
)

# ---- 3. intron compression on and off -------------------------------------
lab <- function(t) plot_annotation(
  title = t,
  theme = theme(plot.title = element_text(size = 9, family = FAM)))
plain <- plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
                      show_gene_label = FALSE) + lab("shrink = FALSE")
shrunk <- plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
                       shrink = TRUE, show_gene_label = FALSE) +
  lab("shrink = TRUE  (ticks label the true coordinates)")
png_out(wrap_elements(plain) / wrap_elements(shrunk), "sashimi-shrink",
        7.2, 9.4)

# ---- 4. every arc geometry ------------------------------------------------
png_out(
  wrap_plots(lapply(arc_shapes(), function(s) {
    sashimi_track(sd, base_family = FAM, group = "Endothelial", arc_shape = s,
                  arc_label_format = "count", arc_width_rule = "log",
                  arc_width = 0.3, psi_pad = 0, group_label = NA) +
      ggtitle(sprintf('arc_shape = "%s"', s)) +
      theme(plot.title = element_text(size = 8, hjust = 0, family = FAM))
  }), ncol = 2),
  "arc-shapes", 8.4, 5.2
)

# ---- 5. label placement and panel background ------------------------------
on <- plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
                   arc_label_position = "on", show_gene_label = FALSE) +
  lab('arc_label_position = "on"  -  the box interrupts the curve')
above <- plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
                      arc_label_position = "above", background = "#F2F4F6",
                      show_gene_label = FALSE) +
  lab('arc_label_position = "above",  background = "#F2F4F6"')
png_out(wrap_elements(on) / wrap_elements(above), "sashimi-labels", 7.2, 9.4)

# ---- 6. overlay, log axis, arcs below -------------------------------------
png_out(
  plot_sashimi(sd, base_family = FAM, overlay = list(Endothelial = "Endothelial",
                                  `Epithelial + Mesenchymal` =
                                    c("Epithelial", "Mesenchymal")),
               alpha = 0.55, arc_label_format = "count",
               arc_width_rule = "log", arc_width = 0.3),
  "sashimi-overlay", 7.2, 4.2
)
png_out(
  plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
               log_y = TRUE, palette = "ember"),
  "sashimi-log", 7.2, 4.9
)
png_out(
  plot_sashimi(sd, base_family = FAM, aggregate = "mean", arc_label_format = "count",
               arc_side = "below", arc_height_frac = c(0.5, 0.75),
               palette = "okabe"),
  "sashimi-below", 7.2, 4.9
)

# ---- 7. the rMATS event ---------------------------------------------------
cat("reading the rMATS event...\n")
ev <- sashimi_from_rmats(ex("SE.MATS.JC.txt"), BAMS, min_count = 5)
png_out(
  plot_sashimi(ev, base_family = FAM, aggregate = "mean", arc_label_format = "count",
               arc_width_rule = "log", arc_width = 0.3),
  "sashimi-rmats", 7.2, 4.6
)

# ---- 8. genome-track view -------------------------------------------------
cat("drawing genome tracks:\n")
png_out(
  plot_tracks(aggregate_tracks(sd, "mean"), title = "ABI1",
              base_family = FAM),
  "tracks-default", 8.5, 3.4, saver = ggsave
)

# ---------------------------------------------------------------------------
# A simulated start-site locus, in the shape the `tss` preset was built for:
# 5'-tag coverage (a flat body with a sharp peak at each active start site),
# two isoforms sharing a downstream anchor, and an arc from each start site
# labelled with that site's activity.
#
# Everything here is invented - the gene name, the contig, the coordinates, the
# condition names and every number. It is not a real locus and the values
# describe nothing. The activities are chosen to span the range that
# format_activity() switches precision across, and to tell a clean switch
# story: the main site dominates early and the alternative one takes over.
# ---------------------------------------------------------------------------
set.seed(7)
WIN <- c(1000000, 1024000)
ANCHOR <- 1001500
ALT <- 1018000
MAIN <- 1021500
pos <- seq(WIN[1], WIN[2], by = 25)

peak <- function(at, h, sd_bp = 90) h * exp(-0.5 * ((pos - at) / sd_bp)^2)
stage_track <- function(main_h, alt_h) {
  pmax(0, peak(MAIN, main_h) + peak(ALT, alt_h) +
         stats::runif(length(pos)) * 0.35)
}
stages <- c("Early", "Mid", "Late")
# Per stage: the height of the main peak, of the alternative peak, and the two
# activities the arcs are labelled with.
sim <- list(
  Early = list(h = c(62, 0.3), act = c(main = 182.4, alt = 0.004)),
  Mid   = list(h = c(18, 25),  act = c(main = 41.3,  alt = 57.9)),
  Late  = list(h = c(4, 47),   act = c(main = 6.2,   alt = 133.5))
)

# Two isoforms differing only in where they start: six shared exons, and a
# first exon at the main or the alternative site.
shared <- data.frame(
  start = c(1001500, 1004200, 1006800, 1009100, 1012000, 1015300),
  end   = c(1002010, 1004310, 1006895, 1009190, 1012100, 1015390)
)
mk_model <- function(tx, role, apex) {
  rbind(
    data.frame(locus_id = "g1", tx_id = tx, role = role, feature = "exon",
               start = shared$start, end = shared$end),
    data.frame(locus_id = "g1", tx_id = tx, role = role, feature = "exon",
               start = apex - 190, end = apex)
  )
}

tss <- sashimi_data(
  loci = data.frame(locus_id = "g1", gene_name = "Demo1", chrom = "7",
                    strand = "-", win_lo = WIN[1], win_hi = WIN[2],
                    anchor = ANCHOR, main_apex = MAIN, alt_apex = ALT),
  tracks = do.call(rbind, lapply(stages, function(s) {
    data.frame(locus_id = "g1", group = s, pos = pos,
               value = stage_track(sim[[s]]$h[1], sim[[s]]$h[2]))
  })),
  junctions = do.call(rbind, lapply(stages, function(s) {
    data.frame(locus_id = "g1", group = s, x0 = ANCHOR, x1 = c(MAIN, ALT),
               role = c("main", "ATSS"), count = unname(sim[[s]]$act))
  })),
  models = rbind(mk_model("main isoform", "main", MAIN),
                 mk_model("ATSS isoform", "ATSS", ALT)),
  meta = list(group_order = stages)
)
tss <- compute_psi(tss, main = "main", alt = "ATSS")

png_out(plot_sashimi(tss, base_family = FAM, preset = "tss"), "sashimi-tss", 5.9, 4.85)

# The same simulated locus as a genome-track view, so the two idioms can be
# compared side by side.
stage_cols <- c(Early = "#66C2A5", Mid = "#FC8D62", Late = "#8DA0CB")
png_out(
  plot_tracks(
    tracks = c(
      list(
        track_axis(),
        track_models(tss, title = "main / ATSS isoforms", color = OM_ROLE_FILL,
                     row_by = "role", labels = TRUE),
        track_features(
          data.frame(start = c(MAIN, ALT), end = c(MAIN + 1, ALT + 1),
                     name = c("main TSS", "ATSS"),
                     color = unname(OM_ROLE_FILL[c("main", "ATSS")])),
          title = "TSS apex", color = "bed_rgb", shape = "marker"
        ),
        track_spacer(0.05)
      ),
      lapply(stages, function(s) {
        track_coverage(tss, title = s, group = s, color = stage_cols[[s]],
                       ymax = max(tss$tracks$value))
      })
    ),
    region = sprintf("7:%d-%d", WIN[1], WIN[2]), axis = "none",
    title = "Demo1 - simulated alternative start site", chrom_style = "ucsc",
    base_family = FAM
  ),
  "tracks-tss", 8.5, 3.6, saver = ggsave
)

# ---- 9. the consequence figures -------------------------------------------
cat("drawing consequences:\n")
cons <- utils::read.delim(ex("demo_consequence.tsv"))
png_out(
  plot_consequence(cons, base_family = FAM, filter = both_full_length == 1,
                   title = "What the start-site switch does to the transcript"),
  "consequence-donut", 6.6, 5.0
)
png_out(plot_consequence(cons, base_family = FAM, filter = both_full_length == 1,
                         style = "lollipop"),
        "consequence-lollipop", 6.8, 3.4)
png_out(plot_consequence(cons, base_family = FAM, filter = both_full_length == 1, style = "bar"),
        "consequence-bar", 6.8, 3.4)

# ---- 10. the palettes -----------------------------------------------------
pal <- do.call(rbind, lapply(omakase_palettes(), function(nm) {
  cols <- omakase_palette(nm)
  data.frame(palette = nm, i = seq_along(cols), colour = cols)
}))
pal$palette <- factor(pal$palette, levels = rev(omakase_palettes()))
png_out(
  ggplot(pal, aes(x = i, y = palette, fill = I(colour))) +
    geom_tile(width = 0.92, height = 0.72) +
    scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_omakase_axes(9, FAM) +
    theme(axis.line = element_blank(), axis.ticks = element_blank()),
  "palettes", 5.6, 2.2
)

cat("\ndone.\n")
