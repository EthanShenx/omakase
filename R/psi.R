# ---------------------------------------------------------------------------
# Percent spliced in.
#
# Two flavours are in circulation and they are not the same number:
#   - the plain activity ratio, main / (main + alternative), which is what the
#     start-site activity figures print;
#   - rMATS's length-corrected inclusion level, which divides each count by the
#     number of positions at which a read could have produced it.
# omakase computes both and says which it used.
# ---------------------------------------------------------------------------

#' Compute PSI (percent spliced in) per locus and group
#'
#' @details
#' The `"ratio"` method is the share of activity carried by the reference
#' feature,
#' \deqn{\psi = \frac{A_{\mathrm{main}}}{A_{\mathrm{main}} + A_{\mathrm{alt}}},}
#' which is well defined for TSS activities, junction counts, or anything else
#' additive. Where both terms are zero, \eqn{\psi} is `NA` rather than 0 - no
#' signal is not the same as no inclusion.
#'
#' The `"rmats"` method is the length-corrected inclusion level
#' \deqn{\psi = \frac{I / \ell_I}{I / \ell_I + S / \ell_S},}
#' where \eqn{I} and \eqn{S} are inclusion and skipping counts and
#' \eqn{\ell_I}, \eqn{\ell_S} the effective lengths, defaulting to the values
#' rMATS itself uses when they are not supplied.
#'
#' @param x A `sashimi_data` object whose `junctions` slot carries a `role`
#'   column, or a data frame with columns `locus_id`, `group`, `role` and
#'   `count`.
#' @param main,alt The `role` values treated as the reference and alternative
#'   features.
#' @param method `"ratio"` (the default) or `"rmats"`.
#' @param len_inc,len_skip Effective lengths \eqn{\ell_I}, \eqn{\ell_S} used by
#'   `method = "rmats"`.
#'
#' @return For a `sashimi_data` input, the object with its `psi` slot filled
#'   in. For a data frame, a data frame of `locus_id`, `group`, `psi`,
#'   `numerator`, `denominator`.
#'
#' @examples
#' j <- data.frame(
#'   locus_id = "a", group = c("early", "early", "late", "late"),
#'   role = c("main", "alt", "main", "alt"), count = c(90, 10, 20, 80)
#' )
#' compute_psi(j)
#'
#' @export
compute_psi <- function(x, main = "main", alt = c("alt", "ATSS", "alternative"),
                        method = c("ratio", "rmats"),
                        len_inc = NULL, len_skip = NULL) {
  method <- match.arg(method)
  is_sd <- inherits(x, "sashimi_data")
  j <- if (is_sd) x$junctions else as_df(x)

  if (!nrow(j)) {
    if (is_sd) return(x)
    return(empty_df(SLOT_SPEC$psi))
  }
  require_cols(j, c("locus_id", "group", "count"), "junction table")
  if (!has_col(j, "role") || all(is.na(j$role))) {
    om_abort(c(
      "PSI needs a {.field role} column marking which junctions are inclusion and which are skipping.",
      "i" = "Set {.code role} to {.val {main}} / {.val {alt[1]}} on the junctions slot."
    ))
  }

  key <- paste(j$locus_id, j$group, sep = "\r")
  sum_role <- function(roles) {
    v <- ifelse(j$role %in% roles, j$count, 0)
    tapply(v, key, sum, na.rm = TRUE)
  }
  I <- sum_role(main)
  S <- sum_role(alt)
  keys <- names(I)

  if (identical(method, "rmats")) {
    li <- len_inc %||% 1
    ls <- len_skip %||% 1
    num <- I / li
    den <- num + S / ls
  } else {
    num <- I
    den <- I + S
  }
  psi <- ifelse(den > 0, num / den, NA_real_)

  parts <- strsplit(keys, "\r", fixed = TRUE)
  out <- data.frame(
    locus_id = vapply(parts, `[`, character(1), 1),
    group = vapply(parts, `[`, character(1), 2),
    psi = as.numeric(psi),
    numerator = as.numeric(num),
    denominator = as.numeric(den),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  if (!is_sd) return(out)
  x$psi <- normalise_slot(out, "psi")
  x$meta$psi_method <- method
  x
}

#' Difference in PSI between two groups
#'
#' The quantity rMATS reports as `IncLevelDifference`: \eqn{\Delta\psi =
#' \psi_{g_1} - \psi_{g_2}}, per locus.
#'
#' @param x A `sashimi_data` object with a populated `psi` slot, or a data
#'   frame with `locus_id`, `group`, `psi`.
#' @param group1,group2 The two group labels to compare.
#' @return A data frame with `locus_id`, `psi_1`, `psi_2` and `dpsi`.
#' @examples
#' p <- data.frame(locus_id = c("a", "a"), group = c("early", "late"),
#'                 psi = c(0.9, 0.2))
#' delta_psi(p, "early", "late")
#' @export
delta_psi <- function(x, group1, group2) {
  p <- if (inherits(x, "sashimi_data")) x$psi else as_df(x)
  require_cols(p, c("locus_id", "group", "psi"), "psi table")
  a <- p[p$group == group1, c("locus_id", "psi")]
  b <- p[p$group == group2, c("locus_id", "psi")]
  names(a)[2] <- "psi_1"
  names(b)[2] <- "psi_2"
  out <- merge(a, b, by = "locus_id", all = TRUE)
  out$dpsi <- out$psi_1 - out$psi_2
  out[order(-abs(out$dpsi)), , drop = FALSE]
}
