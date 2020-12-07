#' Get outgroup(s) node
#'
#' Determine the number of the node corresponding to one or more outgroup
#' sequences in a phylogenetic tree.
#'
#' @param tree_obj A phylogenetic tree, in `ape` phylo format.
#' @param og A character vector comprising the names of one or more sequences
#'   designated as outgroup(s).
#' @return An integer.
#' @examples
#' get_ognode()
#' @export
get_ognode <- function (tree_obj, og) {
  if (length(og) > 1) {  # if multiple outgroups are specified,
    # determine the MRCA node of the outgroups:
    phytools::findMRCA(tree_obj, tips = og, type = "node")
  } else { # else, if a single outgroup is specified,
    # extract the number of the corresponding tip label:
    which(tree_obj$tip.label == og[1])
  }
}