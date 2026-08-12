#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom ggplot2 .pt
#' @importFrom rlang .data %||% abort arg_match inform warn
#' @importFrom stats median approx setNames quantile
#' @importFrom utils head tail read.delim write.table modifyList
## usethis namespace: end
NULL

# Silence R CMD check notes for columns referenced inside dplyr/ggplot2 verbs
# that are resolved from the data mask rather than the calling environment.
utils::globalVariables(c("."))
