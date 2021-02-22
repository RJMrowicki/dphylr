#' Convert a vector into a factor with sequential levels
#'
#' Convert a character or numeric vector into a factor with sequential
#' pseudo-numeric levels, beginning at a specific 'number'.
#'
#' @param .x A pre-ordered character or numeric vector.
#' @param start_num The 'number' at which the sequence of factor levels should
#'   start (default 0).
#' @return A factor with pseudo-numeric levels beginning at the specified
#'   number.
#' @examples
#' recode_seq()
#' @export
recode_seq <- function(.x, start_num = 0) {
  # convert the vector into a factor (preserve order of levels):
  var_fct <- as.character(.x) %>% forcats::as_factor()

  level_key <-  # create a named vector for recoding factor levels
    levels(var_fct) %>%  # extract factor levels
    # create
    magrittr::set_names(seq_along(.)-1 + start_num)

  # recode the factor using the named vector:
  forcats::fct_recode(var_fct, !!!level_key)
}
