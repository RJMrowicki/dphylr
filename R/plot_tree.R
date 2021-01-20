#' Plot a tree
#'
#' Plot a labelled and formatted phylogenetic tree.
#'
#' @param tree_obj A phylogenetic tree, in `ape` phylo format.
#' @param og A character vector comprising the names of one or more sequences
#'   designated as outgroup(s).
#' @param shorten_og Should the outgroup(s) branch(es) be shortened?
#' @param node_vals_lab A table containing node numbers (`node`) and
#'   corresponding character string labels (`lab`). The node numbers must relate
#'   to those in `tree_obj`.
#' @param samp_seq_lab A table of sequence codes (`seq_code`; must be the first
#'   column) and corresponding tip labels (`tiplab`) as parsable character
#'   strings.
#' @param xext Factor by which to increase the plot width relative to its
#'   original size (a value of 1 will result in an increase of 100%). This value
#'   also determines the size and position of the species delimitation matrix,
#'   if plotted.
#' @param lht Line height (in cm) for tree tips (and the species delimitation
#'   heatmap, if plotted).
#' @return A `ggtree` plot object, including associated data table.
#' @examples
#' plot_tree()
#' @export
plot_tree <- function (
  tree_obj, og, shorten_og, node_vals_lab, samp_seq_lab,
  xext = 0.5, lht = 0.5
) {
  p0 <-  # initiate the plot and assign to object
    # (NB -- reduce line thickness; plot ladderised tree in correct orientation)
    ggtree::ggtree(tree_obj, size = 0.4, right = TRUE)

  if (shorten_og) {  # if the outgroup(s) branch is to be shortened:
    og_node <- get_ognode(tree_obj, og)  # extract the outgroup(s) node number
    og_node_x <- dplyr::filter(p0$data, node==og_node)$x  # original og node xpos
    og_node_x_new <- og_node_x / 2  # shorten the node branch by 50%

    # replace original value in plot data with new value:
    p0$data %<>% dplyr::mutate(x = replace(x, node == og_node, og_node_x_new))

    if (length(og) > 1) {  # if there are multiple outgroups,
      og_x <- dplyr::filter(p0$data, label %in% og)$x  # original og branch xpos(s)
      og_bl <- og_x - og_node_x  # original outgroup branch length(s)
      non_og_x <- dplyr::filter(p0$data, !label %in% og)$x  # non-og branch xpos(s)
      # shorten og branch length(s), relative to mean xpos for other branches:
      og_bl_new <- og_bl * (mean(non_og_x) / mean(og_x))
      og_x_new <- og_node_x_new + og_bl_new  # convert into xpos

      # replace original values in plot data with new values:
      p0$data %<>% dplyr::mutate(x = replace(x, label %in% og, og_x_new))
    }

    # apply outgroup clade grouping structure to tree, in order to use different
    # line styles for outgroup(s) (NB -- `groupOTU()` instead of `groupClade()`?):
    p0 %<>% tidytree::groupClade(get_ognode(tree_obj, og))

    # apply different line style for outgroup (specify manually):
    p0 <- p0 +
      ggplot2::aes(linetype = group, show.legend = FALSE) +
      ggplot2::scale_linetype_manual(values = c("solid", "longdash")) +
      ggplot2::guides(linetype = FALSE)  # (NB -- avoid legend in `gheatmap`)
  }


  # add combined ML boot and BI pp node labels to plot data:
  p0$data %<>% dplyr::left_join(
    dplyr::select(node_vals_lab, node, nodelab), by = "node")


  # extract current x limits:
  xlim0 <- ggplot2::ggplot_build(p0)$layout$panel_params[[1]]$x.range
  xdist0 <- diff(xlim0)  # calculate x distance
  # extract current y limits:
  ylim0 <- ggplot2::ggplot_build(p0)$layout$panel_params[[1]]$y.range
  ydist0 <- diff(ylim0)  # calculate y distance


  p <-  # extend the plot object
    # map data frame to tree for annotation
    # (NB -- requires sequence codes [= tip labels] to be in first column):
    p0 %<+% dplyr::select(samp_seq_lab, seq_code, tidyselect::everything()) +

    # disable axis limit extension, expand b and r margins:
    ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
    ggplot2::theme(plot.margin = margin(b = 0.75, unit = "cm")) +
    ggplot2::xlim(NA, xlim0[2] + xdist0 * xext) +  # extend x limits by `xext`

    ggtree::geom_treescale(  # add scale bar
      x = 0, width = signif(0.15*xdist0, 1), linesize = 1) +

    ggtree::geom_tiplab(  # add tip labels
      ggplot2::aes(  # `subset = grp != 0` for differential formatting
        subset = grp != 0, label = tiplab),  # non group0 samples
      parse = TRUE, align = FALSE, size = 3) +
    ggtree::geom_tiplab(  # add tip labels
      ggplot2::aes(
        subset = grp == 0, label = tiplab),  # group0 samples
      parse = TRUE, align = FALSE, size = 3, colour = "black") +

    geom_nodelab(  # add customised node labels
      aes(label = nodelab),  # use 'x = branch' for label on branch
      hjust = 1.11, vjust = -0.5, size = 2.5)
    # ggrepel::geom_text_repel(  # add non-overlapping node labels
    #   ggplot2::aes(label = nodelab),  # use 'x = branch' for label on branch
    #   na.rm = TRUE,
    #   box.padding = 0.01,  # default 0.25
    #   # # for `ggrepel::geom_label_repel()`:
    #   # label.padding = 0.1, label.size = 0, fill = adjustcolor("white", alpha = 0.75),
    #   nudge_x = -0.005*xdist0,  # nudge_y = 0.0025*ydist0,
    #   direction = "y", hjust = 0.25, vjust = -0.25, segment.size = 0.1, size = 2.5)


  p +  # plot tree only
    ggplot2::theme(plot.margin = margin(t = lht/2, unit = "cm"))
}
