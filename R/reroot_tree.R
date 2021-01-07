#' Re-root a tree
#'
#' Re-root a phylogenetic tree at a position half-way along the outgroup branch.
#'
#' @param tree_obj A phylogenetic tree, in `ape` phylo format.
#' @param og A character vector comprising the names of one or more sequences
#'   designated as outgroup(s).
#' @param pos The position along the target edge at which the tree is to be
#'   re-rooted, expressed as a number between 0 and 1.
#' @return A phylogenetic tree, in `ape` phylo format.
#' @examples
#' reroot_tree()
#' @export
reroot_tree <- function (tree_obj, og, pos = 0.5) {
  og_node <-  # create object for storing outgroup(s) node number
    get_ognode(tree_obj, og)
  og_edge_length <-  # create object for storing outgroup(s) edge length
    # extract `edge.length` corresponding to node in `edge` column 2:
    with(tree_obj, edge.length[edge[,2] == og_node])

  # re-root tree, at a position half-way (pos = 0.5) along the target edge:
  # (NB -- following same approach as old version of `ggtree::reroot()`)
  phytools::reroot(tree_obj, position = pos * og_edge_length, node = og_node)
}
