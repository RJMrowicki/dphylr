#' Generate diagnostics for BEAST2 results
#'
#' Calculate potential scale reduction factor (PSRF) values and plot average
#' standard deviation of split frequency (ASDSF) values for one or more BEAST2
#' MCMC runs, in order to diagnose convergence.
#'
#' @param beast_dir Path to the directory containing the results of the BEAST2
#'   analysis (`.log` and `.trees` files).
#' @param burnin A numeric value specifying the proportion of MCMC trace values
#'   or trees to be discarded as burnin (default = 0.25).
#' @param plot_asdsf If TRUE, a plot of ASDSF values (calculated cumulatively
#'   across all MCMC chains) will be produced.
#' @return A list containing the following elements: `beast_logs`, a
#'   `coda::mcmc.list` object containing parameter trace values (minus burnin)
#'   for each MCMC chain; `beast_trees`, a `rwty::rwty.chain` object containing
#'   the trees (as `ape::multiPhylo` objects) and associated trace log data for
#'   each MCMC chain; and `beast_psrf`, the multivariate PSRF values calculated
#'   for each parameter.
#' @examples
#' diag_beast()
#' @export
diag_beast <- function(
  beast_dir, burnin = 0.25, plot_asdsf = TRUE
) {
  # ~ PSRF:
  beast_logs <-  # create list of BEAST2 MCMC chain logs
    # generate vector of log file names, including paths:
    fs::dir_ls(beast_dir, glob = "*.log") %>%
    purrr::map(., ~ {  # parse each log file and convert to `coda::mcmc` object:
      dat <- tracerer::parse_beast_tracelog_file(.)[, -1]  # NB -- remove `Sample`
      tracerer::remove_burn_ins(dat, burnin) %>%  # discard % burnin
        coda::mcmc(., start = burnin*nrow(dat), end = nrow(dat))  # specify limits
    }) %>%
    coda::mcmc.list()  # convert list to `coda::mcmc.list` object

  # calculate combined PSRF values (Gelman-Rubin convergence diagnostic):
  beast_psrf <- coda::gelman.diag(  # (NB - 'autoburnin = FALSE'?)
    beast_logs, autoburnin = FALSE, multivariate = FALSE
  )  # (should all be ~1)

  # ~ ASDSF:
  rwty.processors <<- 4  # assign no. of processors (`rwty::`)

  beast_trees <-  # create list of BEAST2 MCMC trees + logs
    # generate vector of tree file names, including paths:
    fs::dir_ls(beast_dir, glob = "*.trees") %>%
    # append "END;" final line to all tree files, otherwise cannot be read:
    purrr::walk(., ~ { readr::write_lines("END;", ., append = TRUE) }) %>%
    # load each tree file (and associated log file) into `rwty.chain` object:
    purrr::map2(., fs::dir_ls(beast_dir, glob = "*.log"), ~ {
      rwty::load.trees(., type = "nexus", trim = 10, logfile = .y)
    })  # (NB -- thin by retaining every 10 [or 100] trees???)

  if (plot_asdsf) {  # if specified, create ASDSF plot
    beast_trees %>%  # (NB -- discard % burnin)
      rwty::makeplot.asdsf(., burnin = round(burnin * length(.[[1]]$trees)))
    # (should reach asymptote of <0.01)
  }

  # return list containing the BEAST2 logs and trees:
  return(tibble::lst(beast_logs, beast_trees, beast_psrf))
}
