#' Conduct maximum likelihood phylogenetic analysis (RAxML-NG)
#'
#' Create a bash script for conducting maximum likelihood phylogenetic analysis
#' via RAxML-NG.
#'
#' @param raxmlng_dir Path to the target directory for the RAxML-NG analysis,
#'   which should contain the multiple sequence alignment as a FASTA file, in
#'   addition to the RAxML-NG partition definitions file.
#' @param aln_file Name of the FASTA file containing the multiple sequence
#'   alignment.
#' @param part_file Name of the file containing the RAxML-NG partition
#'   definitions.
#' @param rseed An integer to use as the random seed for analyses. If
#'   unspecified, a random seed will be generated.
#' @param nboot An integer specifying the number of bootstrap replicates for the
#'   analysis (default 1e+03).
#' @param bscr_dir Path to the directory in which the RAxML-NG bash script will
#'   be created.
#' @param run If TRUE, the bash script is run in Windows Subsystem for Linux
#'   (WSL) as a native R process.
#' @return No objects are returned. A bash script named `raxmlng.sh` is created
#'   in the specified script directory. If run, this script performs the
#'   RAxML-NG analysis, with results saved in the specified RAxML-NG directory.
#' @examples
#' phyl_raxmlng()
#' @export
phyl_raxmlng <- function (
  raxmlng_dir, aln_file, part_file, rseed = NULL, nboot = 1e+03, bscr_dir, run = FALSE
) {
  # if random seed is unspecified, generate random seed:
  if (is.null(rseed)) { rseed <- round(runif(1, max = 10^4)) }

  # create bash script for running RAxML-NG:
  readr::write_lines(  # write the following lines...
    c(
      "#!/bin/bash",  # bash script header
      "# Maximum likelihood tree searching (RAxML-NG)",
      paste0('cd "',raxmlng_dir,'"'),  # change directory to alignment
      "rm RAxML-NG_*",  # remove files from previous analysis

      paste0(  # check alignment for formatting errors (read 'raxmlng_0check.raxml.log'):
        'raxml-ng-mpi ',  # (MPI-enabled version > PTHREADS version?)
        '--msa ',aln_file,' --prefix RAxML-NG_0check ',  # input file, output prefix
        '--parse ',
        # (NB -- using 'parse' instead of 'check' estimates optimal no. of cores)
        # substitution model (if concatenated, specify partition file):
        paste0('--model ',part_file,' ')), # ifelse(marker == "conc", "raxmlng_parts", use_mod_raxmlng)

      # # A. Combined ML search + bootstrapping: (convenient for small datasets)
      # paste0(  # ML tree search and bootstrapping:
      #   'raxml-ng-mpi ',
      #   '--all --msa ',aln_file,' --prefix RAxML-NG_1all ',  # input file and output prefix
      #   # substitution model (if concatenated, specify partition file):
      #   paste0('--model ',part_file,' '), # ifelse(marker == "conc", part_file, use_mod_raxmlng)
      #   paste0('--threads 4 --force --seed ',rseed,' '),  # no. threads (forced), random seed
      #   '--tree rand{50},pars{50} ',  # (default rand{10},pars{10})
      #   # perform 100 tree searches using {n} random and {n} parsimony-based starting trees
      #   # (NB -- for partitioned analyses, RAxML-NG estimates 'scaled' branch lengths
      #   # by default, which is considered most appropriate, compared to 'linked' or 'unlinked'
      #   # [RAxML 8.x can only estimate 'linked' or 'unlinked', with 'linked' as default].)
      #   # bootstrapping (autoMRE{n} / 5000?):
      #   paste0('--bs-trees ',format(nboot, scientific = FALSE),' --bs-cutoff 0.03')
      # ),
      #
      # paste0(  # test for bootstrap convergence (read 'raxmlng_converge.raxml.log'):
      #   '# ',  # NOT RUN
      #   'raxml-ng-mpi ',
      #   # specify input (bootstrap trees from previous analysis) and output prefix:
      #   '--bsconverge --bs-trees RAxML-NG_1all.raxml.bootstraps --prefix RAxML-NG_2bsconverge ',
      #   paste0('--threads 4 --force --seed ',rseed,' '),  # no. threads (forced), random seed
      #   '--bs-cutoff 0.03'  # make cutoff more stringent? (default 0.03, i.e. 3%)
      # )

      # B. Separate ML search, bootstrapping & support mapping: (large datasets?)
      # (NB -- allows opportunity to re-run individual parts of the analysis)
      paste0(  # perform ML tree search:
        'raxml-ng-mpi ',
        '--msa ',aln_file,' --prefix RAxML-NG_1run ',  # input file and output prefix
        # substitution model (if concatenated, specify partition file):
        paste0('--model ',part_file,' '), # ifelse(marker == "conc", part_file, use_mod_raxmlng)
        paste0('--threads 4 --force --seed ',rseed,' '),  # no. threads (forced), random seed
        '--tree rand{50},pars{50}'  # (default rand{10},pars{10})
      ),

      paste0(  # perform standard bootstrapping:
        'raxml-ng-mpi ',
        '--msa aln.fasta --prefix RAxML-NG_2boot ',  # input file and new output prefix
        # substitution model (if concatenated, specify partition file):
        paste0('--model ',"raxmlng_parts",' '), # ifelse(i=="conc", "raxmlng_parts", use_mod_raxmlng)
        paste0('--threads 4 --force --seed ',rseed,' '),  # no. threads (forced), random seed
        # bootstrapping (autoMRE{n} / nboot?):
        paste0('--bootstrap --bs-trees ',format(nboot, scientific = FALSE))
      ),

      paste0(  # test for bootstrap convergence (read 'raxmlng_converge.raxml.log'):
        'raxml-ng-mpi ',
        # specify input (bootstrap trees from previous analysis) and output prefix:
        '--bsconverge --bs-trees RAxML-NG_2boot.raxml.bootstraps --prefix RAxML-NG_3converge ',
        paste0('--threads 4 --force --seed ',rseed,' '),  # no. threads (forced), random seed
        '--bs-cutoff 0.03'  # make cutoff more stringent? (default 0.03, i.e. 3%)
      ),

      paste0(  # map bootstrap support values to best-scoring ML tree:
        'raxml-ng-mpi ',
        # specify best-scoring ML tree and bootstrap tree files:
        '--support --tree RAxML-NG_1run.raxml.bestTree --bs-trees RAxML-NG_2boot.raxml.bootstraps ',
        '--prefix RAxML-NG_4support'  # output prefix
      )
    ),
    # ... to the following file:
    paste0(bscr_dir,"raxmlng.sh")
  )

  if (run) {  # if specified, run the script via WSL:
    system(paste0('bash -c "',bscr_dir,'raxmlng.sh"'))
  }
}
