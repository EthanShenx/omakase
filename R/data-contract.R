# ---------------------------------------------------------------------------
# The sashimi_data object.
#
# Every reader in the package (BAM, rMATS, tag BED, junction file, tidy tables)
# produces one of these, and every plotting function consumes one. Keeping a
# single documented intermediate representation is what lets a user swap the
# input format without touching the plotting call, and lets them build the
# object by hand when their counts come from somewhere omakase has never heard
# of.
# ---------------------------------------------------------------------------

# Column specifications for each slot. Used both to validate user input and to
# manufacture correctly-typed empty slots.
SLOT_SPEC <- list(
  loci = c(locus_id = "character", gene_name = "character",
           chrom = "character", strand = "character",
           win_lo = "numeric", win_hi = "numeric"),
  tracks = c(locus_id = "character", group = "character",
             pos = "numeric", value = "numeric"),
  junctions = c(locus_id = "character", group = "character",
                x0 = "numeric", x1 = "numeric", count = "numeric"),
  models = c(locus_id = "character", tx_id = "character",
             start = "numeric", end = "numeric"),
  psi = c(locus_id = "character", group = "character", psi = "numeric"),
  features = c(locus_id = "character", start = "numeric", end = "numeric",
               name = "character")
)

# Optional columns each slot understands, with the default used when absent.
SLOT_OPTIONAL <- list(
  loci = list(anchor = NA_real_, main_apex = NA_real_, alt_apex = NA_real_),
  tracks = list(sample = NA_character_, strand = "*", bin = NA_real_),
  junctions = list(role = NA_character_, label = NA_character_,
                   sample = NA_character_, strand = "*"),
  models = list(role = NA_character_, feature = "exon", strand = "*"),
  psi = list(numerator = NA_real_, denominator = NA_real_),
  features = list(role = NA_character_, class = NA_character_,
                  color = NA_character_)
)

#' Construct a sashimi data object
#'
#' The container every `sashimi_from_*()` reader returns and every plotting
#' function accepts. It is a plain list of tidy data frames, so it can be built
#' by hand when counts come from a source omakase does not read directly.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{`loci`}{One row per plotted locus. Required: `locus_id`,
#'     `gene_name`, `chrom`, `strand`, `win_lo`, `win_hi`. Optional: `anchor`
#'     (a shared downstream point arcs may point at), `main_apex`, `alt_apex`
#'     (start sites to mark with a triangle).}
#'   \item{`tracks`}{Binned coverage. Required: `locus_id`, `group`, `pos`,
#'     `value`. Optional: `sample`, `strand`, `bin`.}
#'   \item{`junctions`}{Arcs. Required: `locus_id`, `group`, `x0`, `x1`,
#'     `count`. Optional: `role`, `label`, `sample`, `strand`. `x0`/`x1` are
#'     genomic coordinates; `count` drives the arc label, and optionally its
#'     width and height.}
#'   \item{`models`}{Transcript models drawn beneath the tracks. Required:
#'     `locus_id`, `tx_id`, `start`, `end`. Optional: `role` (used to colour
#'     and to order the rows), `feature` (`"exon"`, `"CDS"`, `"UTR"`),
#'     `strand`.}
#'   \item{`psi`}{Per-group inclusion values printed in the right-hand gutter.
#'     Required: `locus_id`, `group`, `psi`. Optional: `numerator`,
#'     `denominator`.}
#'   \item{`features`}{An extra annotation row under each model - repeats,
#'     motifs, peaks, anything. Required: `locus_id`, `start`, `end`, `name`.
#'     Optional: `role`, `class`, `color`.}
#'   \item{`meta`}{A named list recording how the object was built. Printed by
#'     `summary()` and written to the methods table by [write_sashimi_data()].}
#' }
#'
#' @param loci,tracks,junctions,models,psi,features Data frames as described
#'   above. Any may be `NULL`.
#' @param meta A named list of provenance/parameter values.
#'
#' @return An object of class `sashimi_data`.
#'
#' @examples
#' loci <- data.frame(
#'   locus_id = "demo", gene_name = "DEMO", chrom = "chr1",
#'   strand = "+", win_lo = 1000, win_hi = 2000
#' )
#' tracks <- data.frame(
#'   locus_id = "demo", group = "A",
#'   pos = seq(1000, 2000, by = 10),
#'   value = abs(sin(seq(0, pi, length.out = 101))) * 10
#' )
#' sashimi_data(loci = loci, tracks = tracks)
#'
#' @export
sashimi_data <- function(loci = NULL, tracks = NULL, junctions = NULL,
                         models = NULL, psi = NULL, features = NULL,
                         meta = list()) {
  x <- list(
    loci = normalise_slot(loci, "loci"),
    tracks = normalise_slot(tracks, "tracks"),
    junctions = normalise_slot(junctions, "junctions"),
    models = normalise_slot(models, "models"),
    psi = normalise_slot(psi, "psi"),
    features = normalise_slot(features, "features"),
    meta = as.list(meta)
  )
  structure(x, class = "sashimi_data")
}

# Coerce one slot: fill in optional columns, check required ones, and cast the
# types the plot layer assumes.
#' @noRd
normalise_slot <- function(df, slot) {
  spec <- SLOT_SPEC[[slot]]
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0L)) {
    return(empty_df(spec))
  }
  df <- as_df(df)
  require_cols(df, names(spec), paste0("`", slot, "` table"))

  for (nm in names(spec)) {
    df[[nm]] <- switch(spec[[nm]],
      character = as.character(df[[nm]]),
      numeric = as.numeric(df[[nm]]),
      df[[nm]]
    )
  }
  for (nm in names(SLOT_OPTIONAL[[slot]])) {
    if (is.null(df[[nm]])) df[[nm]] <- SLOT_OPTIONAL[[slot]][[nm]]
  }
  # A junction drawn right-to-left is the same junction; normalising the
  # orientation here means the arc code never has to test which end is which.
  if (identical(slot, "junctions") && nrow(df)) {
    flip <- !is.na(df$x0) & !is.na(df$x1) & df$x0 > df$x1
    if (any(flip)) {
      tmp <- df$x0[flip]
      df$x0[flip] <- df$x1[flip]
      df$x1[flip] <- tmp
    }
  }
  rownames(df) <- NULL
  df
}

#' Validate a sashimi data object
#'
#' Checks slot classes and cross-slot referential integrity: every `locus_id`
#' appearing in `tracks`, `junctions`, `models`, `psi` or `features` must exist
#' in `loci`.
#'
#' @param x A `sashimi_data` object.
#' @param strict If `TRUE`, an unknown `locus_id` is an error; otherwise a
#'   warning and the offending rows are reported but kept.
#' @return `x`, invisibly.
#' @examples
#' validate_sashimi_data(sashimi_data())
#' @export
validate_sashimi_data <- function(x, strict = TRUE) {
  if (!inherits(x, "sashimi_data")) {
    om_abort("Expected a {.cls sashimi_data} object, got {.cls {class(x)[1]}}.")
  }
  known <- unique(x$loci$locus_id)
  for (slot in c("tracks", "junctions", "models", "psi", "features")) {
    d <- x[[slot]]
    if (!nrow(d)) next
    unknown <- setdiff(unique(d$locus_id), known)
    if (length(unknown)) {
      msg <- c("Slot {.field {slot}} references locus_id{?s} absent from {.field loci}: {.val {unknown}}.")
      if (strict) om_abort(msg) else om_warn(msg)
    }
  }
  if (nrow(x$loci) && any(x$loci$win_hi <= x$loci$win_lo, na.rm = TRUE)) {
    bad <- x$loci$locus_id[which(x$loci$win_hi <= x$loci$win_lo)]
    om_abort("Loc{?us/i} {.val {bad}} ha{?s/ve} a window whose end is not after its start.")
  }
  invisible(x)
}

#' @export
print.sashimi_data <- function(x, ...) {
  n_loci <- nrow(x$loci)
  cli::cli_text("{.cls sashimi_data}: {n_loci} loc{?us/i}, {length(unique(x$tracks$group))} group{?s}")
  slots <- c("tracks", "junctions", "models", "psi", "features")
  rows <- vapply(slots, function(s) nrow(x[[s]]), integer(1))
  for (i in seq_along(slots)) {
    if (rows[i]) cli::cli_bullets(c("*" = "{.field {slots[i]}}: {rows[i]} row{?s}"))
  }
  if (n_loci) {
    shown <- utils::head(x$loci$gene_name, 6)
    more <- if (n_loci > 6) paste0(", ... (+", n_loci - 6, ")") else ""
    cli::cli_text("loci: {paste(shown, collapse = ', ')}{more}")
  }
  invisible(x)
}

#' @export
summary.sashimi_data <- function(object, ...) {
  print(object)
  if (length(object$meta)) {
    cli::cli_h3("meta")
    for (nm in names(object$meta)) {
      v <- object$meta[[nm]]
      cli::cli_bullets(c("*" = "{.field {nm}}: {paste(format(v), collapse = ', ')}"))
    }
  }
  invisible(object)
}

#' Subset a sashimi data object by locus
#'
#' @param x A `sashimi_data` object.
#' @param i Locus identifiers, gene names, or a numeric/logical index into
#'   `loci`.
#' @param ... Unused.
#' @return A `sashimi_data` object holding only the selected loci.
#' @examples
#' sd <- sashimi_data(loci = data.frame(
#'   locus_id = c("a", "b"), gene_name = c("A", "B"), chrom = "chr1",
#'   strand = "+", win_lo = 1, win_hi = 10
#' ))
#' sd["a"]
#' @export
`[.sashimi_data` <- function(x, i, ...) {
  loci <- x$loci
  keep <- if (is.character(i)) {
    loci$locus_id %in% i | loci$gene_name %in% i
  } else {
    seq_len(nrow(loci)) %in% seq_len(nrow(loci))[i]
  }
  ids <- loci$locus_id[keep]
  out <- x
  out$loci <- loci[keep, , drop = FALSE]
  for (slot in c("tracks", "junctions", "models", "psi", "features")) {
    d <- x[[slot]]
    out[[slot]] <- d[d$locus_id %in% ids, , drop = FALSE]
  }
  out
}

#' Locus identifiers held by a sashimi data object
#'
#' @param x A `sashimi_data` object.
#' @return A character vector of `locus_id` values.
#' @examples
#' loci(sashimi_data())
#' @export
loci <- function(x) {
  if (!inherits(x, "sashimi_data")) om_abort("Expected a {.cls sashimi_data} object.")
  x$loci$locus_id
}

#' Combine sashimi data objects
#'
#' Row-binds matching slots, so several regions read separately can be plotted
#' or written as one.
#'
#' @param ... `sashimi_data` objects.
#' @return A single `sashimi_data` object.
#' @examples
#' a <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
#'   chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
#' b <- sashimi_data(loci = data.frame(locus_id = "b", gene_name = "B",
#'   chrom = "chr2", strand = "-", win_lo = 5, win_hi = 20))
#' c_ab <- combine_sashimi(a, b)
#' loci(c_ab)
#' @export
combine_sashimi <- function(...) {
  objs <- list(...)
  objs <- Filter(Negate(is.null), objs)
  if (!length(objs)) return(sashimi_data())
  bad <- !vapply(objs, inherits, logical(1), "sashimi_data")
  if (any(bad)) om_abort("All arguments must be {.cls sashimi_data} objects.")

  slot_bind <- function(slot) {
    parts <- lapply(objs, function(o) o[[slot]])
    cols <- unique(unlist(lapply(parts, names)))
    parts <- lapply(parts, function(d) {
      for (cn in setdiff(cols, names(d))) d[[cn]] <- NA
      d[, cols, drop = FALSE]
    })
    rbind_all(parts) %||% empty_df(SLOT_SPEC[[slot]])
  }
  out <- sashimi_data(
    loci = slot_bind("loci"), tracks = slot_bind("tracks"),
    junctions = slot_bind("junctions"), models = slot_bind("models"),
    psi = slot_bind("psi"), features = slot_bind("features"),
    meta = Reduce(function(a, b) utils::modifyList(a, b),
                  lapply(objs, function(o) o$meta), list())
  )
  dup <- duplicated(out$loci$locus_id)
  if (any(dup)) out$loci <- out$loci[!dup, , drop = FALSE]
  out
}

#' Write a sashimi data object to disk
#'
#' Writes one file per non-empty slot, plus a `*_methods.tsv` recording `meta`.
#' The layout mirrors what [sashimi_from_tables()] reads back, so a slow
#' counting step can be run once and re-plotted cheaply.
#'
#' @param x A `sashimi_data` object.
#' @param dir Output directory, created if absent.
#' @param prefix File name prefix.
#' @param parquet If `TRUE` and the arrow package is available, write the
#'   (largest) `tracks` slot as parquet instead of TSV.
#' @return The character vector of written paths, invisibly.
#' @examples
#' sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
#'   chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
#' write_sashimi_data(sd, tempdir(), prefix = "demo")
#' @export
write_sashimi_data <- function(x, dir, prefix = "sashimi", parquet = FALSE) {
  validate_sashimi_data(x)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  written <- character()
  for (slot in c("loci", "tracks", "junctions", "models", "psi", "features")) {
    d <- x[[slot]]
    if (!nrow(d)) next
    if (parquet && identical(slot, "tracks") && has_pkg("arrow")) {
      p <- file.path(dir, sprintf("%s_%s.parquet", prefix, slot))
      arrow::write_parquet(d, p)
    } else {
      p <- file.path(dir, sprintf("%s_%s.tsv", prefix, slot))
      utils::write.table(d, p, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    written <- c(written, p)
  }
  if (length(x$meta)) {
    p <- file.path(dir, sprintf("%s_methods.tsv", prefix))
    md <- data.frame(
      parameter = names(x$meta),
      value = vapply(x$meta, function(v) paste(format(v), collapse = ","),
                     character(1)),
      stringsAsFactors = FALSE
    )
    utils::write.table(md, p, sep = "\t", quote = FALSE, row.names = FALSE)
    written <- c(written, p)
  }
  invisible(written)
}
