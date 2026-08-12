# Internal helpers shared across the package. Nothing here is exported.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Exact column access.
#
# `df$n` partial-matches, so on a table that happens to carry `n_uATG_gained`
# it silently returns that column instead of NULL. Every optional column this
# package probes on user-supplied data goes through here, where `[[` matches
# exactly.
#' @noRd
col <- function(df, name, default = NULL) {
  if (is.null(df) || !name %in% names(df)) return(default)
  df[[name]]
}

#' @noRd
has_col <- function(df, name) !is.null(df) && name %in% names(df)

# +1 when arcs are drawn above the coverage, -1 when below, so an offset can be
# written once and applied in whichever direction the arcs run.
#' @noRd
sign_of <- function(side) if (identical(side, "below")) -1 else 1

# Abort with an omakase-branded condition.
#
# cli evaluates the `{}` expressions in a message against `.envir`, which
# defaults to the caller of cli_abort() - that would be this wrapper's own
# frame, where none of the message's variables exist. Forwarding the calling
# function's environment is what makes `{.val {x}}` resolve to the caller's
# `x`, and the same applies to the warning and message helpers below.
#' @noRd
om_abort <- function(msg, ..., class = NULL, call = rlang::caller_env(),
                     .envir = rlang::caller_env()) {
  cli::cli_abort(msg, ..., class = c(class, "omakase_error"), call = call,
                 .envir = .envir)
}

#' @noRd
om_warn <- function(msg, ..., .envir = rlang::caller_env()) {
  cli::cli_warn(msg, ..., .envir = .envir)
}

#' @noRd
om_inform <- function(msg, ..., .envir = rlang::caller_env()) {
  cli::cli_inform(msg, ..., .envir = .envir)
}

# Is a package available? Used to gate the optional Biostrings/arrow paths.
#' @noRd
has_pkg <- function(pkg) isTRUE(requireNamespace(pkg, quietly = TRUE))

#' @noRd
need_pkg <- function(pkg, what) {
  if (!has_pkg(pkg)) {
    om_abort(c(
      "{.pkg {pkg}} is required for {what}.",
      "i" = 'Install it with {.run BiocManager::install("{pkg}")} or {.run install.packages("{pkg}")}.'
    ))
  }
  invisible(TRUE)
}

# Coerce anything data-frame-ish (tibble, data.table, arrow Table) to a plain
# data.frame with stringsAsFactors off, so downstream code has one shape to
# reason about.
#' @noRd
as_df <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "data.frame") && !inherits(x, "tbl_df") &&
      !inherits(x, "data.table")) {
    return(x)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

# Read a delimited table, transparently handling .gz and picking the separator
# from the extension when `sep` is not given.
#' @noRd
read_table_any <- function(path, sep = NULL, ...) {
  if (!file.exists(path)) om_abort("File not found: {.path {path}}.")
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  if (identical(ext, "parquet")) {
    need_pkg("arrow", "reading parquet files")
    return(as_df(arrow::read_parquet(path)))
  }
  sep <- sep %||% if (ext %in% c("csv")) "," else "\t"
  as_df(utils::read.delim(
    path, sep = sep, header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, comment.char = "", quote = "", ...
  ))
}

# Require a set of columns, naming every missing one at once rather than
# failing on the first.
#' @noRd
require_cols <- function(df, cols, what = "table") {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    om_abort(c(
      "The {what} is missing required column{?s} {.field {missing}}.",
      "i" = "Columns present: {.field {names(df)}}."
    ))
  }
  invisible(df)
}

# Rename columns by a named vector c(new = "old"), silently skipping pairs whose
# source column is absent. Used by the adapters to map project-specific column
# names onto the sashimi_data contract.
#' @noRd
rename_cols <- function(df, map) {
  for (new in names(map)) {
    old <- map[[new]]
    if (!is.null(old) && old %in% names(df) && !identical(old, new)) {
      names(df)[match(old, names(df))] <- new
    }
  }
  df
}

# An empty data.frame with a given column spec, so downstream code can rbind or
# iterate over an absent slot without branching.
#' @noRd
empty_df <- function(spec) {
  out <- lapply(spec, function(type) vector(type, 0L))
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Drop a leading "chr" (or add one) so an event file and a BAM header that
# disagree about contig naming can still be joined.
#' @noRd
harmonise_chrom <- function(x, style = c("keep", "ucsc", "ensembl")) {
  style <- match.arg(style)
  switch(style,
    keep = x,
    ucsc = ifelse(grepl("^chr", x), x, paste0("chr", x)),
    ensembl = sub("^chr", "", x)
  )
}

# Recycle a scalar/short vector to length n, erroring when the length is neither
# 1 nor n (a silent recycle here would mislabel samples).
#' @noRd
recycle_to <- function(x, n, what = "value") {
  if (is.null(x)) return(NULL)
  if (length(x) == n) return(x)
  if (length(x) == 1L) return(rep(x, n))
  om_abort("{what} has length {length(x)}; expected 1 or {n}.")
}

# Bind a list of data frames, tolerating NULLs and zero-row members.
#' @noRd
rbind_all <- function(lst) {
  lst <- Filter(function(d) !is.null(d) && nrow(d) > 0L, lst)
  if (!length(lst)) return(NULL)
  do.call(rbind, c(lst, list(make.row.names = FALSE)))
}

# Order a character vector by a preferred order, appending unknown values in
# their natural order rather than dropping them. With no preference the order
# of first appearance wins: stage and condition names ("Early", "Mid", "Late",
# "8-cell") almost never sort into their biological order, and the order they
# were supplied in is the author's intent.
#' @noRd
order_levels <- function(x, preferred = NULL) {
  u <- unique(x)
  if (is.null(preferred)) return(u)
  c(preferred[preferred %in% u], setdiff(u, preferred))
}
