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

  bi_node_pp <-  # create table of BI node no.s and posterior probabilities
    ggplot2::fortify(bi_tree) %>%  # extract data table
    # remove rows corresponding to tip 'nodes' (i.e. retain internal nodes):
    dplyr::filter(isTip == FALSE) %>%
    # select `node` and `label` (rename to `pp`) columns only:
    dplyr::select(node, pp = label)

  ml_node_boot <-  # create table of ML node no.s and bootstrap values
    ggplot2::fortify(ml_tree) %>%  # extract data table
    # remove rows corresponding to tip 'nodes' (i.e. retain internal nodes):
    dplyr::filter(isTip == FALSE) %>%
    # select `node` and `label` (rename to `boot`) columns only:
    dplyr::select(node, boot = label)

  bi_node_boot <-  # create table of BI node no.s and ML bootstrap values
    purrr::map_dfr(bi_st, function(.b) {  # for each BI subtree (i.e. node),
      # (NB -- return final value as a data frame by binding rows)
      # determine the corresponding node number on the original BI tree:
      bi_node <- dphylr::get_ognode(bi_tree, .b$tip.label)
      # determine which ML subtrees (if any) contain matching tip labels:
      match_ml <- purrr::map_lgl(ml_st, ~ { setequal(.b$tip.label, .$tip.label) })
      # (NB -- returns logical string to use as predicate in `purrr::map_if()`)

      # for each ML subtree (i.e. node), if tip labels are equal,
      boot_val <- purrr::map_if(ml_st, .p = match_ml, function(.m) {
        # determine the corresponding node no. on the original ML tree:
        ml_node <- dphylr::get_ognode(ml_tree, .m$tip.label)
        # extract the corresponding ML bootstrap value:
        ml_node_boot$boot[ml_node_boot$node == ml_node]
      }, .else = NA_character_) %>%  # else output character NA (results in NULL)
        flatten_chr %>% .[!is.null(.)]  # extract non-NULL values only

      # output combined BI node no. and ML bootstrap value:
      tibble::tibble(node = bi_node, boot = boot_val)
    })

  suppressWarnings(
    bi_node_pp_boot <-
      # merge BI posterior probabilities and ML bootstrap values:
      dplyr::left_join(bi_node_pp, bi_node_boot, by = "node") %>%
      # format columns as numerical (generates warnings if char. values present):
      dplyr::mutate(dplyr::across(c(pp, boot), as.numeric))
  )

  # output table of node no.s, posterior probabilities and bootstrap values:
  return(bi_node_pp_boot)
}
