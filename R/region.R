#' Parse a genomic region string
#'
#' Accepts the usual browser syntax and the colon-separated form used by
#' `rmats2sashimiplot`, and is forgiving about thousands separators and
#' whitespace.
#'
#' Recognised forms:
#' * `"chr10:27035000-27050000"`
#' * `"chr10:27,035,000-27,050,000"`
#' * `"chr10:27035000..27050000"`
#' * `"chr10:+:27035000:27050000"` (rmats2sashimiplot `-c` style, with strand)
#' * `"chr10"` (whole contig; `start`/`end` become `NA`)
#'
#' @param x A region string, or a `list`/`data.frame` already carrying
#'   `chrom`, `start`, `end` (returned unchanged after validation).
#' @param one_based Logical. If `TRUE` (the default) the returned `start` is
#'   treated as 1-based inclusive, matching the way coordinates are written in
#'   genome browsers and GTF files. Set to `FALSE` to interpret the input as
#'   0-based half-open BED coordinates, which are then converted internally.
#'
#' @return A list with elements `chrom`, `start`, `end`, `strand` and
#'   `label`. `strand` is `"*"` unless the input specified one.
#'
#' @examples
#' parse_region("chr10:27,035,000-27,050,000")
#' parse_region("chr10:+:27035000:27050000")
#'
#' @export
parse_region <- function(x, one_based = TRUE) {
  if (inherits(x, "omakase_region")) return(x)

  if (is.list(x) || is.data.frame(x)) {
    require_cols(as_df(x), c("chrom", "start", "end"), "region")
    return(new_region(
      as.character(x$chrom[1]), as.numeric(x$start[1]), as.numeric(x$end[1]),
      as.character(x$strand[1] %||% "*")
    ))
  }

  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    om_abort("{.arg region} must be a single region string, e.g. {.val chr1:1000-2000}.")
  }

  raw <- trimws(x)
  clean <- gsub("[, ]", "", raw)

  # rmats2sashimiplot form: chrom:strand:start:end
  parts <- strsplit(clean, ":", fixed = TRUE)[[1]]
  if (length(parts) == 4L && parts[2] %in% c("+", "-", "*")) {
    return(new_region(parts[1], as.numeric(parts[3]), as.numeric(parts[4]),
                      parts[2], one_based = one_based, label = raw))
  }
  if (length(parts) == 3L) {
    # Either chrom:start:end, or chrom:start-end:strand.
    if (parts[3] %in% c("+", "-", "*")) {
      se <- strsplit(sub("\\.\\.", "-", parts[2]), "-", fixed = TRUE)[[1]]
      if (length(se) == 2L) {
        return(new_region(parts[1], as.numeric(se[1]), as.numeric(se[2]),
                          parts[3], one_based = one_based, label = raw))
      }
    }
    return(new_region(parts[1], as.numeric(parts[2]), as.numeric(parts[3]),
                      "*", one_based = one_based, label = raw))
  }

  if (length(parts) == 1L) {
    # A bare contig name. Anything with whitespace in it is not one, and is
    # far more likely to be a malformed region than a real sequence name.
    if (grepl("\\s", raw) || !nzchar(clean)) {
      om_abort("Could not parse region {.val {raw}}; expected {.val chr:start-end}.")
    }
    return(new_region(parts[1], NA_real_, NA_real_, "*", label = raw))
  }
  if (length(parts) != 2L) {
    om_abort("Could not parse region {.val {raw}}.")
  }

  chrom <- parts[1]
  span <- sub("\\.\\.", "-", parts[2])
  # A trailing strand suffix, e.g. chr1:100-200:+
  strand <- "*"
  m <- regmatches(span, regexec("^(.*?)([+-])$", span))[[1]]
  if (length(m) == 3L) {
    span <- m[2]
    strand <- m[3]
  }
  se <- strsplit(span, "-", fixed = TRUE)[[1]]
  if (length(se) != 2L || anyNA(suppressWarnings(as.numeric(se)))) {
    om_abort("Could not parse region {.val {raw}}; expected {.val chr:start-end}.")
  }
  new_region(chrom, as.numeric(se[1]), as.numeric(se[2]), strand,
             one_based = one_based, label = raw)
}

#' @noRd
new_region <- function(chrom, start, end, strand = "*", one_based = TRUE,
                       label = NULL) {
  if (!is.na(start) && !is.na(end) && end < start) {
    om_abort("Region end ({end}) is before its start ({start}).")
  }
  # Internally omakase is 1-based inclusive throughout; BED-style input is
  # nudged once, here, so nothing downstream has to remember the convention.
  if (!one_based && !is.na(start)) start <- start + 1

  structure(
    list(
      chrom = as.character(chrom),
      start = start,
      end = end,
      strand = strand,
      label = label %||% sprintf("%s:%s-%s", chrom,
                                 format(start, big.mark = ",", trim = TRUE),
                                 format(end, big.mark = ",", trim = TRUE))
    ),
    class = "omakase_region"
  )
}

#' @export
format.omakase_region <- function(x, ...) {
  s <- if (identical(x$strand, "*")) "" else paste0(" (", x$strand, ")")
  if (is.na(x$start)) return(paste0(x$chrom, " [whole contig]", s))
  sprintf("%s:%s-%s%s", x$chrom,
          format(x$start, big.mark = ",", trim = TRUE),
          format(x$end, big.mark = ",", trim = TRUE), s)
}

#' @export
print.omakase_region <- function(x, ...) {
  cat("<omakase region> ", format(x), "\n", sep = "")
  invisible(x)
}

#' @noRd
region_width <- function(r) {
  if (is.na(r$start)) return(NA_real_)
  r$end - r$start + 1
}

# Widen a region by a flank, clamping the low edge at 1.
#' @noRd
expand_region <- function(r, flank = 0) {
  if (is.na(r$start) || flank == 0) return(r)
  new_region(r$chrom, max(1, r$start - flank), r$end + flank, r$strand)
}

#' @noRd
region_to_granges <- function(r) {
  need_pkg("GenomicRanges", "converting a region to a GRanges")
  GenomicRanges::GRanges(
    seqnames = r$chrom,
    ranges = IRanges::IRanges(start = r$start, end = r$end),
    strand = r$strand
  )
}
