#' Plot a tree
#'
#' Plot a labelled and formatted phylogenetic tree, optionally with a matrix or
#' bars summarising the results of species delimitation analyses.
#'
#' @param tree_obj A phylogenetic tree, in `ape` phylo format.
#' @param og A character vector comprising the names of one or more sequences
#'   designated as outgroup(s).
#' @param shorten_og Should the outgroup(s) branch(es) be shortened?
#' @param node_vals_lab A table containing node numbers (`node`) and
#'   corresponding character string labels (`lab`). The node numbers must relate
#'   to those in `tree_obj`.
#' @param samp_seq_lab A table of either sequence codes (`seq_code`) or sample
#'   codes (`samp_code`, for concatenated sequence data) and corresponding tip
#'   labels (`tiplab`) as parsable character strings.
#' @param match_by Either "seq_code", for matching labels to `tree_obj` tips by
#'   sequence code (default), or "samp_code", for matching labels by sample code
#'   (for concatenated sequence data).
#' @param p_width Plot width (in cm) used to determine the relative width of the
#'   species delimitation matrix, if plotted. For square matrix cells, specify
#'   the same width when using graphics output devices.
#' @param xext Factor by which to increase the plot width relative to its
#'   original size (a value of 1 will result in an increase of 100%). This value
#'   also determines the size and position of the species delimitation matrix,
#'   if plotted.
#' @param lht Line height (in cm) for tree tips (and the species delimitation
#'   matrix, if plotted).
#' @param plot_delim Plot corresponding species delimitation results as either a
#'   group association matrix ("matrix") or clade bars ("bars") to the right of
#'   the tree; otherwise, do not include species delimitation results ("none").
#' @param spp_delim A data frame containing the combined species-sample
#'   association results from one or more delimitation analyses, as output by
#'   `comb_delim()`.
#' @return A `ggtree` plot object, including associated data table.
#' @examples
#' plot_tree()
#' @export
plot_tree <- function (
  tree_obj, og, shorten_og, node_vals_lab,
  samp_seq_lab, match_by = c("seq_code", "samp_code"),
  p_width = 20, xext = 0.5, lht = 0.5,
  plot_delim = c("matrix", "bars", "none"), spp_delim
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
    p0 %<>% tidytree::groupClade(dphylr::get_ognode(tree_obj, og))

    # apply different line style for outgroup (specify manually):
    p0 <- p0 +
      ggplot2::aes(linetype = group, show.legend = FALSE) +
      ggplot2::scale_linetype_manual(values = c("solid", "longdash")) +
      ggplot2::guides(linetype = "none")  # (NB -- avoid legend in `gheatmap`)
  }


  # add combined ML boot and BI pp node labels to plot data:
  p0$data %<>% dplyr::left_join(
    dplyr::select(node_vals_lab, node, nodelab), by = "node"
  )


  # extract current x limits:
  xlim0 <- ggplot2::ggplot_build(p0)$layout$panel_params[[1]]$x.range
  xdist0 <- diff(xlim0)  # calculate x distance
  # extract current y limits:
  ylim0 <- ggplot2::ggplot_build(p0)$layout$panel_params[[1]]$y.range
  ydist0 <- diff(ylim0)  # calculate y distance


  p <-  # extend the plot object
    # map data frame to tree for annotation
    # (NB -- requires seq/samp codes [= tip labels] to be in first column):
    p0 %<+% dplyr::select(samp_seq_lab, tidyselect::all_of(match_by), tidyselect::everything()) +

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

    ggtree::geom_nodelab(  # add customised node labels
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


  # if plotting a species delimitation matrix (heatmap),
  if (plot_delim == "matrix") {
    # extract group association variables from delimitation results table and
    # **convert `seq_code` col to rownames** (required by `ggtree::gheatmap()`)
    sp_grps <- spp_delim %>% tibble::column_to_rownames("seq_code")
    # extract overall number of levels (should be the same for all variables,
    # representing the maximum number of groups delimited by any analysis):
    n_sp_grps <- length(levels(dplyr::pull(sp_grps, 1)))

    sp_colours <-  # create custom colour palette for heatmap
      colorRampPalette(  # define custom colour sequence
        rep(  # 2 colours, alternating between each adjacent group
          c(grey(0.45), grey(0.85)), times = ceiling( n_sp_grps / 2 )
        )[1:n_sp_grps],
        alpha = TRUE, bias = 1  # return alpha values, alter colour spread
      )

    # specify total width of species matrix *relative* to tree width:
    hm_width <- (lht/p_width) * (1 + xext) * ncol(sp_grps)
    # (desired single column width [cm; = lht] / plot width [cm] *
    # tree width [1] + x extension width [relative to tree width]) *
    # no. of heatmap columns)

    # specify species matrix offset in terms of *x scale*:
    hm_offset <- xdist0*xext - (hm_width*xdist0)
    # (max x limit - matrix width [converted to x units])

    p %>%  # plot tree and associated species delimitation matrix
      ggtree::gheatmap(
        data = sp_grps,
        width = hm_width, offset = hm_offset, color = "white",
        colnames_angle = -45, hjust = 0, font.size = 3) +
      ggplot2::scale_fill_manual(  # set matrix colours
        breaks = as.character(seq_len(n_sp_grps)),
        # breaks = ggplot2::waiver(),  # use default breaks (scale limits)
        values = sp_colours(n_sp_grps),  # colour values from custom palette
        na.value = NA,  # ensure NA values are blank (default "grey50")
        guide = "none")
  }

  # otherwise, if plotting species delimitation bars (clade groups),
  if (plot_delim == "bars") {
    # make a list of tree nodes representing species delimitation groups:
    node_grps <- spp_delim %>%
      # for all variables except `seq_code`,
      purrr::map_at(., dplyr::vars(-seq_code), ~ {
        # combine `seq_code` and variable (= delimitation group) into tibble:
        tibble::tibble(seq_code = spp_delim$seq_code, delim_grp = .) %>%
          # remove NAs (i.e. outgroup rows):
          tidyr::drop_na() %>%
          # split `seq_code` into list by delimitation group:
          dplyr::group_by(delim_grp) %>%
          # for each group of seq_codes,
          dplyr::group_map(., ~ {
            # obtain the tree node corresponding to that group:
            dphylr::get_ognode(tree_obj, dplyr::pull(.))
          }) %>%
          # flatten list and convert into tibble column:
          unlist %>% dplyr::as_tibble() %>%
          # convert rownames (= group) to column, and rename columns:
          tibble::rownames_to_column() %>% rlang::set_names("grp", "nnode") %>%
          # add dummy column for blank `ggtree::geom_cladelab()` labels:
          tibble::add_column(grp_lab = NA)
      }) %>%
      # remove element `seq_code` from the list:
      purrr::list_modify("seq_code" = NULL)

    # obtain the taxon name corresponding to the last (bottom) tip:
    last_taxon <- dplyr::last(ggtree::get_taxa_name(p))
    # obtain the corresponding tree node:
    lab_node <- dphylr::get_ognode(tree_obj, last_taxon)

    # set bar offset increment (= distance between bar groups):
    bar_offset_increment <- 0.02 * xdist0*(1 + xext)
    # (2% * total plot width (in terms of *x units*))
    # set initial offset for clade bars:
    bar_offset <- xdist0*xext - ((length(node_grps)-1) * bar_offset_increment)
    # (max x limit - ((no. sets of delimitation groups - 1) * bar offset increment))

    for (i in names(node_grps)) {  # for each set of delimitation groups,
      p <- p +
        ggtree::geom_cladelab(  # add clade labels (bar only)
          data = node_grps[[i]],
          mapping = aes(node = nnode, label = grp_lab),
          barcolor = grey(0.25), barsize = 3, extend = 0.35,
          textcolour = "transparent", align = TRUE, offset = bar_offset
        ) +
        ggtree::geom_cladelab(  # add delimitation method label
          node = lab_node, label = paste0("     ", i), angle = -90,
          barsize = 0, align = TRUE, offset = bar_offset,
          fontsize = 2.5, offset.text = 0, vjust = 1
        )
      # increase clade bar offset by increment:
      bar_offset <- bar_offset + bar_offset_increment
    }

    p +  # plot tree (set top and bottom margins)
      ggplot2::theme(plot.margin = margin(t = lht/2, b = lht, unit = "cm"))
  }

  # otherwise, if not including species delimitation results:
  if (plot_delim == "none") {
    p +  # plot tree only (set top margin)
      ggplot2::theme(plot.margin = margin(t = lht/2, unit = "cm"))
  }
}
