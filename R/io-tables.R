# ---------------------------------------------------------------------------
# Building a sashimi_data object from tables the user already has.
#
# This is the escape hatch that keeps the package useful when counts come from
# somewhere omakase does not read: a caller who can produce tidy data frames
# can draw the figure, regardless of how the numbers were made. It is also how
# a slow BAM pass is separated from fast re-plotting.
# ---------------------------------------------------------------------------

#' Build a sashimi data object from tidy tables
#'
#' Accepts data frames or paths to TSV/CSV/parquet files, one per slot, and
#' maps their columns onto the [sashimi_data()] contract. Column names that
#' differ from the contract can be remapped through `rename`, so tables written
#' by an existing pipeline usually need no editing.
#'
#' @param loci,tracks,junctions,models,psi,features Data frames, or paths to
#'   files. Any may be `NULL`.
#' @param rename A named list of named character vectors, one per slot, mapping
#'   contract column names to the names used in your tables. For example
#'   `rename = list(tracks = c(group = "stage", value = "rpm"))` reads a table
#'   whose columns are `stage` and `rpm`.
#' @param locus_col The column identifying the locus, if it is not `locus_id`.
#'   A common case is a gene name column shared by every table.
#' @param meta A named list recorded on the object.
#'
#' @return A `sashimi_data` object.
#'
#' @examples
#' loci <- data.frame(gene_name = "DEMO", chrom = "chr1", strand = "+",
#'                    win_lo = 1000, win_hi = 2000)
#' tracks <- data.frame(gene_name = "DEMO", stage = "early",
#'                      pos = seq(1000, 2000, 100), rpm = runif(11))
#' sashimi_from_tables(
#'   loci = loci, tracks = tracks,
#'   rename = list(tracks = c(group = "stage", value = "rpm")),
#'   locus_col = "gene_name"
#' )
#'
#' @export
sashimi_from_tables <- function(loci = NULL, tracks = NULL, junctions = NULL,
                                models = NULL, psi = NULL, features = NULL,
                                rename = list(), locus_col = NULL,
                                meta = list()) {
  get_slot <- function(v, slot) {
    if (is.null(v)) return(NULL)
    d <- if (is.character(v) && length(v) == 1L) read_table_any(v) else as_df(v)
    if (!is.null(rename[[slot]])) d <- rename_cols(d, as.list(rename[[slot]]))
    if (!is.null(locus_col) && !has_col(d, "locus_id") && has_col(d, locus_col)) {
      d$locus_id <- as.character(d[[locus_col]])
    }
    d
  }

  loci <- get_slot(loci, "loci")
  if (!is.null(loci)) {
    # A table keyed only by gene name still needs a locus_id, and the gene name
    # is the natural one.
    if (!has_col(loci, "locus_id") && has_col(loci, "gene_name")) {
      loci$locus_id <- as.character(loci$gene_name)
    }
    if (!has_col(loci, "gene_name") && has_col(loci, "locus_id")) {
      loci$gene_name <- as.character(loci$locus_id)
    }
  }

  sashimi_data(
    loci = loci,
    tracks = get_slot(tracks, "tracks"),
    junctions = get_slot(junctions, "junctions"),
    models = get_slot(models, "models"),
    psi = get_slot(psi, "psi"),
    features = get_slot(features, "features"),
    meta = c(list(source = "tables"), meta)
  )
}

#' Read a directory of sashimi tables
#'
#' The counterpart to [write_sashimi_data()]: reads back the per-slot files it
#' wrote. Files are matched by the slot name appearing in the file name, so
#' both `sashimi_tracks.tsv` and a prefixed `mystudy_tracks.parquet` are found.
#'
#' @param dir Directory to read.
#' @param prefix Optional file name prefix to require.
#' @return A `sashimi_data` object.
#' @examples
#' d <- tempfile(); dir.create(d)
#' sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
#'   chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
#' write_sashimi_data(sd, d)
#' read_sashimi_dir(d)
#' @export
read_sashimi_dir <- function(dir, prefix = NULL) {
  if (!dir.exists(dir)) om_abort("Directory not found: {.path {dir}}.")
  files <- list.files(dir, full.names = TRUE)
  pick <- function(slot) {
    pat <- paste0(if (is.null(prefix)) "" else paste0("^", prefix, ".*"),
                  slot, "\\.(tsv|csv|txt|parquet)(\\.gz)?$")
    hit <- files[grepl(pat, basename(files))]
    if (!length(hit)) NULL else hit[1]
  }
  args <- lapply(c("loci", "tracks", "junctions", "models", "psi", "features"),
                 pick)
  names(args) <- c("loci", "tracks", "junctions", "models", "psi", "features")

  meta <- list(source = "directory", dir = dir)
  mf <- pick("methods")
  if (!is.null(mf)) {
    md <- read_table_any(mf)
    if (all(c("parameter", "value") %in% names(md))) {
      meta <- c(meta, stats::setNames(as.list(md$value), md$parameter))
    }
  }
  do.call(sashimi_from_tables, c(args, list(meta = meta)))
}

#' Coerce an object to sashimi data
#'
#' Used by the plotting functions so they accept a `sashimi_data` object, a
#' list of slot tables, or a directory path interchangeably.
#'
#' @param x A `sashimi_data` object, a named list of slot data frames, or a
#'   path to a directory of tables.
#' @param ... Passed to [sashimi_from_tables()].
#' @return A `sashimi_data` object.
#' @examples
#' as_sashimi_data(list(loci = data.frame(locus_id = "a", gene_name = "A",
#'   chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10)))
#' @export
as_sashimi_data <- function(x, ...) {
  if (inherits(x, "sashimi_data")) return(x)
  if (is.character(x) && length(x) == 1L && dir.exists(x)) {
    return(read_sashimi_dir(x))
  }
  if (is.list(x) && !is.data.frame(x)) {
    known <- c("loci", "tracks", "junctions", "models", "psi", "features",
               "meta")
    return(do.call(sashimi_from_tables,
                   c(x[intersect(names(x), known)], list(...))))
  }
  om_abort("Cannot interpret {.cls {class(x)[1]}} as {.cls sashimi_data}.")
}

#' Read a sample manifest
#'
#' Reads the tab-separated manifest used to describe a set of alignment files.
#' The format matches the one `ggsashimi` accepts, so an existing manifest
#' works unchanged: column 1 is a sample identifier, column 2 the path to the
#' file, and any further columns are metadata that can be named in `group_col`
#' or `color_col`.
#'
#' A header row is optional. When absent, columns are named `sample`, `path`,
#' `V3`, `V4` and so on.
#'
#' @param path Path to the manifest, or a character vector of file paths, or a
#'   single alignment file.
#' @param group_col Column (name or 1-based index) that assigns samples to
#'   groups. Defaults to column 3 when present, otherwise every sample is its
#'   own group.
#' @param label_col Column carrying the label to print on each panel. Defaults
#'   to the sample identifier.
#' @param base_dir Directory that relative paths in the manifest are resolved
#'   against; defaults to the manifest's own directory.
#'
#' @return A data frame with columns `sample`, `path`, `group`, `label`, plus
#'   any extra manifest columns.
#'
#' @examples
#' m <- tempfile(fileext = ".tsv")
#' writeLines(c("s1\ta.bam\tEndothelial", "s2\tb.bam\tEpithelial"), m)
#' read_manifest(m)
#'
#' @export
read_manifest <- function(path, group_col = NULL, label_col = NULL,
                          base_dir = NULL) {
  # A bare list of files, or one file, is a manifest with no metadata.
  if (length(path) > 1L || (length(path) == 1L && grepl("\\.(bam|sam|cram)$",
                                                        path, ignore.case = TRUE))) {
    s <- sub("\\.(bam|sam|cram)$", "", basename(path), ignore.case = TRUE)
    return(data.frame(sample = s, path = path, group = s, label = s,
                      stringsAsFactors = FALSE))
  }
  if (!file.exists(path)) om_abort("Manifest not found: {.path {path}}.")

  first <- utils::read.delim(path, header = FALSE, nrows = 1, sep = "\t",
                             stringsAsFactors = FALSE, comment.char = "#")
  # A manifest's second column is a path. If it is a number, this is a data
  # file (a junction table, a tag BED) that was handed straight to a reader,
  # and it describes one unnamed sample rather than a list of them.
  if (ncol(first) >= 2 &&
      !is.na(suppressWarnings(as.numeric(as.character(first[[2]][1]))))) {
    s <- tools::file_path_sans_ext(basename(path))
    return(data.frame(sample = s, path = path, group = s, label = s,
                      stringsAsFactors = FALSE))
  }
  # Treat the first row as a header only if it names the path column rather
  # than pointing at a file.
  has_header <- ncol(first) >= 2 &&
    tolower(trimws(as.character(first[[2]][1]))) %in%
      c("path", "file", "bam", "filename", "bam_file")
  d <- utils::read.delim(path, header = has_header, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")
  if (!has_header) {
    names(d)[seq_len(min(2, ncol(d)))] <- c("sample", "path")[seq_len(min(2, ncol(d)))]
  } else {
    names(d)[1:2] <- c("sample", "path")
  }
  d$sample <- as.character(d$sample)
  d$path <- as.character(d$path)

  base_dir <- base_dir %||% dirname(path)
  rel <- !grepl("^([/~]|[A-Za-z]:)", d$path)
  d$path[rel] <- file.path(base_dir, d$path[rel])

  pick_col <- function(spec, default) {
    if (is.null(spec)) return(default)
    if (is.numeric(spec)) {
      if (spec > ncol(d)) om_abort("Manifest has no column {spec}.")
      return(as.character(d[[spec]]))
    }
    if (!spec %in% names(d)) om_abort("Manifest has no column {.field {spec}}.")
    as.character(d[[spec]])
  }
  default_group <- if (ncol(d) >= 3) as.character(d[[3]]) else d$sample
  d$group <- pick_col(group_col, default_group)
  d$label <- pick_col(label_col, d$sample)
  d
}
