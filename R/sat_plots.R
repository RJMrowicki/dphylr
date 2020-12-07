#' Create saturation plots
#'
#' Create nucleotide substitution saturation plots of raw (uncorrected) vs.
#' model-corrected pairwise distances and the number of transitions and
#' transversions (expressed as a proportion of sequence length) for sequences in
#' a multiple alignment, based on one or more models of sequence evolution.
#'
#' @param aln_obj A `phangorn` phyDat object containing the multiple sequence
#'   alignment.
#' @param sub_mods A character vector containing the names of one or more
#'   nucleotide substitution model(s) for which plots are to be generated. These
#'   must correspond with model options used by `ape::dist.dna()`.
#' @param plot_title The main title for the plot.
#' @return ...
#' @examples
#' sat_plots()
#' @export
sat_plots <- function (
  aln_obj, sub_mods = c("JC69", "F84", "GG95"), plot_title = NULL
) {
  # convert alignment to `ape` DNAbin format:
  aln_dat <- ape::as.DNAbin(aln_obj)
  
  suppressMessages(  # suppress message from `.name_repair`
    plot_dat <-  # create data frame to store plot data
      purrr::map_dfc(c("raw", "TS", "TV"), ~ {  # for each 'distance' model
        # ('raw' = uncorrected dist., 'TS' = #transitions, 'TV' = #transversions)
        ape::dist.dna(aln_dat, .) %>% as.numeric  # calculate pairwise 'distances'
      }) %>%  # (`purrr::map_dfc()` combines vectors into a table)
      purrr::set_names(c("dist_raw", "TS", "TV")) %>%  # define column names
      # convert no.s of transitions/transversions into proportions:
      dplyr::mutate(dplyr::across(  # (NB -- use `across()`, not `mutate_at()`)
        c("TS", "TV"), ~ { ./ncol(aln_dat) }, .names = "{.col}_prop"
      ), .keep = "unused")  # keep existing variables not used in the function
  )
  
  plot_list <-  # create list to store plots for each substitution model
    purrr::map(sub_mods, ~ {  # for each substitution model,
      # calculate model-corrected pairwise distances:
      dist_mod <- ape::dist.dna(aln_dat, .) %>% as.vector
      
      # create plot of uncorrected (raw) vs. corrected distance:
      plot_raw <-  # assign to object
        tibble::tibble(plot_dat, dist_mod) %>%  # combine with plot data table
        ggplot2::ggplot(aes(dist_mod, dist_raw)) +  # create plot
        # set x and y axis lower limits to 0, expand upper limit by 2.5%:
        ggplot2::scale_x_continuous(
          expand = expansion(c(0, 0.025)), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(
          expand = expansion(c(0, 0.025)), limits = c(0, NA)) +
        ggplot2::labs(x = NULL, y = "Uncorrected distance") +  # label y axis
        ggplot2::geom_point(colour = grey(0.4)) +  # add points (grey)
        ggplot2::geom_abline(  # add 1:1 indicator line (grey, dashed)
          intercept = 0, slope = 1, colour = "grey", linetype = 2) +
        ggplot2::geom_smooth(  # add linear regression line (grey, solid)
          method = "lm", formula = y ~ x, se = FALSE, na.rm = TRUE,
          colour = grey(0.4), weight = 1.5, fullrange = TRUE) +
        theme_custom() +  # apply custom theme
        theme(  # specify plot margins
          plot.margin = ggplot2::margin(1, 1, 4, 4))
      
      # create plot of transitions/transversions vs. distance:
      plot_tstv <-  # assign to object
        tibble::tibble(plot_dat, dist_mod) %>%  # combine with plot data table
        tidyr::pivot_longer(  # pivot into long table according to TS/TV vars
          ends_with("_prop"), names_to = "var", values_to = "y") %>%
        ggplot2::ggplot(aes(dist_mod, y, colour = var)) +  # create plot
        scale_colour_manual(  # set group labels and colours
          labels = c("TS", "TV"), values = c("red", "blue")) +
        # set x and y axis lower limits to 0, expand upper limit by 2.5%:
        ggplot2::scale_x_continuous(
          expand = expansion(c(0, 0.025)), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(
          expand = expansion(c(0, 0.025)), limits = c(0, NA)) +
        ggplot2::labs(x = ., y = "Ts and Tv") +  # label x and y axes
        ggplot2::geom_point() +  # add points (group colours)
        ggplot2::geom_smooth(  # add linear regression lines (group colours)
          method = "lm", formula = y ~ x, se = FALSE, na.rm = TRUE,
          weight = 1.5, fullrange = TRUE) +
        theme_custom()  +  # apply custom theme
        theme(  # position legend in top left, specify plot margins
          legend.position = c(0.025, 1), legend.justification = c(0, 1),
          plot.margin = ggplot2::margin(1, 1, 4, 4))
      
      # return plots bound together as a 'grob' (aligns left axes):
      gridExtra::gtable_rbind(ggplotGrob(plot_raw), ggplotGrob(plot_tstv))
    })

  # arrange plots together in a grid:
  gridExtra::grid.arrange(
    grobs = plot_list, nrow = 1, as.table = FALSE, top = plot_title
  )
}