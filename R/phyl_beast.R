#' Conduct ultrametric phylogenetic analysis (BEAST2)
#'
#' Create a bash script for conducting ultrametric (time-calibrated)
#' phylogenetic analysis via BEAST2.
#'
#' @details An existing BEAST2 config (`.xml`) file is required. This config
#'   file can be generated using the `beautier` package (part of the [`babette`
#'   suite](https://docs.ropensci.org/babette/)); however, `beautier` does NOT
#'   support >1 site, clock or tree models. Therefore, for partitioned data, the
#'   config file must be generated via the `BEAUti.exe` programme.
#'
#' @param beast_dir Path to the target directory for the BEAST2 analysis, which
#'   should contain the relevant config file generated via `BEAUTi.exe`.
#' @param xml_file Name of the BEAST2 config (`.xml`) file.
#' @param trees_file Name of the BEAST2 treelog (`.trees`) file specified in the
#'   config file.
#' @param bscr_dir Path to the directory in which the BEAST2 bash script will be
#'   created.
#' @param run If TRUE, the bash script is run in Windows Subsystem for Linux
#'   (WSL) as a native R process.
#' @return No objects are returned. A bash script named `beast.sh` is created in
#'   the specified script directory. If run, this script performs the BEAST2
#'   analysis, with results saved in the specified BEAST2 directory.
#' @examples
#' phyl_beast()
#' @export
phyl_beast <- function (
  beast_dir, xml_file, trees_file, bscr_dir, run = FALSE
) {
  # create bash script for running BEAST2:
  readr::write_lines(  # write the following lines...
    c(
      "#!/bin/bash", # bash script header
      "# Time-calibrated Bayesian phylogenetic analysis (BEAST2)",
      paste0('cd "',beast_dir,'"'),  # change directory
      "rm *.trees", "rm *.log", "rm *state",  # remove files from previous analysis
      # set BEAGLE library path variable (home directory),
      # (otherwise BEAST2 cannot locate BEAGLE library & produces errors):
      "export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH",

      paste(  # run BEAST2: (`-beagle_SSE` fastest)
        "~/beast/bin/beast -beagle_SSE -instances 2 -threads 4 -overwrite",
        xml_file),

      # # NB -- may use 'logcombiner' (same directory as 'beast' &
      # # 'treeannotator') to combine trace/tree logs from multiple runs; for
      # # coupled MCMC, use trace/tree logs from **cold** (non-heated) chain
      # # only, i.e. chain 0.
      # paste0(  # combine traces
      #   "~/beast/bin/logcombiner -log run1_aln.log -log run2_aln.log ",
      #   "-b 25 -o aln.log"),  # (NB -- 25% burnin; but check traces in Tracer first)
      # paste0(  # combine trees
      #   "~/beast/bin/logcombiner -log run1_aln.trees -log run2_aln.trees ",
      #   "-b 25 -o aln.trees"),  # (NB -- 25% burnin)

      paste(  # summarise trees (TreeAnnotator):
        # NB -- 0 burnin, if using trees combined from separate runs (otherwise 25):
        "~/beast/bin/treeannotator -b 0 -lowMem",  # avoid exceeding memory
        # posterior probability limit 0, keep target node heights:
        "-limit 0 -heights keep", trees_file, "aln_mcc.tre")
      # (NB -- 'aln.trees' may be called 'tree.trees', if treelog filename
      # was not changed from default value in BEAUti)
    ),
    # ... to the following file:
    paste0(bscr_dir,"beast.sh")
  )

  if (run) {  # if specified, run the script via WSL:
    system(paste0('bash -c "',bscr_dir,'raxmlng.sh"'))
  }

  # check output to ensure convergence:
  # - open 'cold' chain log files in Tracer, and for each one (and combined),
  # examine trace and estimated sample size (ESS; should be >200) for every
  # parameter;
  # - examine standard deviation of split frequencies (ASDSF; should be <0.01)
  # and potential scale reduction factor (PSRF; should be ~1) in R.
}
