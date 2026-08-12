#' Format junction or activity values with adaptive precision
#'
#' A single number of decimal places never suits a column of junction counts or
#' TSS activities: the values span several orders of magnitude, so rounding to
#' whole numbers prints small non-zero values as `"0"`, while a fixed two
#' decimals makes large counts unreadable. `format_activity()` picks the
#' precision from the magnitude of each value, and prints an explicit
#' `"<0.01"` so that a small non-zero value is never shown as absent.
#'
#' The rule is:
#'
#' \deqn{f(v) = \begin{cases}
#'   \texttt{"0"}      & v = 0 \\
#'   \texttt{"<0.01"}  & 0 < v < 0.01 \\
#'   \mathrm{fixed}(v, 2) & 0.01 \le v < 1 \\
#'   \mathrm{fixed}(v, 1) & 1 \le v < 10 \\
#'   \mathrm{comma}(\mathrm{round}(v)) & v \ge 10
#' \end{cases}}
#'
#' @param v Numeric vector.
#' @param big_mark Thousands separator used for values of 10 or more.
#' @param zero String used for exact zeros.
#'
#' @return A character vector the same length as `v`.
#'
#' @examples
#' format_activity(c(0, 0.004, 0.37, 4.2, 103.01, 4218))
#'
#' @export
format_activity <- function(v, big_mark = ",", zero = "0") {
  v <- as.numeric(v)
  out <- ifelse(
    is.na(v), NA_character_,
    ifelse(v == 0, zero,
    ifelse(abs(v) < 0.01, "<0.01",
    ifelse(abs(v) < 1, sprintf("%.2f", v),
    ifelse(abs(v) < 10, sprintf("%.1f", v),
           format(round(v), big.mark = big_mark, trim = TRUE, scientific = FALSE))))))
  out
}

#' Format an integer count
#'
#' Junction read counts are whole numbers, so they get thousands separators and
#' nothing else. Used as the default `arc_label_format` when values come from a
#' BAM rather than from a normalised activity table.
#'
#' @param v Numeric vector.
#' @param big_mark Thousands separator.
#' @return A character vector.
#' @examples
#' format_count(c(3, 45, 1200))
#' @export
format_count <- function(v, big_mark = ",") {
  ifelse(is.na(v), NA_character_,
         format(round(as.numeric(v)), big.mark = big_mark, trim = TRUE,
                scientific = FALSE))
}

#' Format a genomic coordinate
#'
#' @param v Numeric vector of coordinates.
#' @param big_mark Thousands separator.
#' @return A character vector.
#' @examples
#' format_coord(27035000)
#' @export
format_coord <- function(v, big_mark = ",") {
  format(round(as.numeric(v)), big.mark = big_mark, trim = TRUE,
         scientific = FALSE)
}

# Resolve the `arc_label_format` / `psi_format` arguments, which accept a
# function or one of a few shorthand strings.
#' @noRd
resolve_formatter <- function(x, default = format_activity) {
  if (is.null(x)) return(default)
  if (is.function(x)) return(x)
  if (is.character(x) && length(x) == 1L) {
    return(switch(x,
      activity = format_activity,
      count    = format_count,
      integer  = format_count,
      coord    = format_coord,
      none     = function(v) as.character(v),
      om_abort("Unknown label format {.val {x}}; use {.val activity}, {.val count}, {.val coord}, {.val none}, or a function.")
    ))
  }
  om_abort("Label format must be a function or a string.")
}
