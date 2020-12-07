#' Match posterior probabilities with bootstrap values
#'
#' Find the node labels (= bootstrap values) of a Maximum Likelihood
#' phylogenetic tree that correspond to the node numbers and labels (= posterior
#' probabilities) of a Bayesian tree generated from the same data.
#'
#' @param bi_tree A Bayesian Inference phylogenetic tree (i.e. the tree to be
#'   plotted, defining the nodes onto which the corresponding Maximum Likelihood
#'   bootstrap values are to be mapped), in `ape` phylo format, with node labels
#'   formatted as numeric.
#' @param ml_tree A Maximum Likelihood phylogenetic tree (generated from the
#'   same data as the Bayesian tree), in `ape` phylo format, with node labels
#'   formatted as numeric.
#' @return A table of Bayesian tree node numbers and posterior probabilities,
#'   with corresponding Maximum Likelihood bootstrap values.
#' @examples
#' match_pp_boot()
#' @export
match_pp_boot <- function (bi_tree, ml_tree) {
  # extract subtrees for BI and ML trees:
  bi_st <- ape::subtrees(bi_tree)
  ml_st <- ape::subtrees(ml_tree)
  
  node_pp_boot <-  # create table of node no.s and BI posterior probabilities
    dplyr::tibble(
      node = seq_len(bi_tree$Nnode) + (1 + bi_tree$Nnode),  # node numbers
      pp = bi_tree$node.label  # node posterior probabilities
    )
  
  node_pp_boot$boot <-  # add column for corresponding ML bootstrap values
    purrr::map_dbl(bi_st, ~ {  # for each BI subtree (i.e. node),
      # (NB -- return final value as a numeric vector)
      purrr::imap_dbl(ml_st, function(.m, .n) {  # for each ML subtree (i.e. node),
        # (NB -- uses `seq_along(.x` as second argument; return numeric vector)
        if (setequal(.x$tip.label, .m$tip.label)) {  # if tip labels are equal,
          ml_tree$node.label[.n]  # extract the corresponding ML node label
        } else { NA_real_ }  # else output numerical NA (required for `_dbl`)
      }) %>% .[!is.na(.)] %>%  # extract non-NA values only
        .[1]  # return the value or, where the value is `integer(0)`, NA
    })
  
  # output table of node no.s, posterior probabilities and bootstrap values:
  return(node_pp_boot)
}