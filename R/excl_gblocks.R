#' Exclude ambiguous alignment sites (Gblocks)
#'
#' Create a bash script for masking ambiguous alignment sites via Gblocks.
#'
#' @param bscr_dir Path to the directory in which the Gblocks bash script will
#'   be created.
#' @param aln_dir Path to the directory containing the multiple sequence
#'   alignment, as a single FASTA file.
#' @param aln_file Name of the FASTA file containing the multiple sequence
#'   alignment.
#' @param run If TRUE, the bash script is run via Windows Subsystem for Linux
#'   (WSL) as a native R process.
#' @return No objects are returned. A bash script named `gblocks.sh` is created
#'   in the specified script directory. If run, this script performs the Gblocks
#'   analysis, with results saved in the specified alignment directory.
#' @examples
#' excl_gblocks()
#' @export
excl_gblocks <- function (bscr_dir, aln_dir, aln_file, run = FALSE) {
  readr::write_lines(
    c(  # output the following lines...:
      "#!/bin/bash", # bash script header
      "# Exclude ambiguous alignment data (Gblocks)",  # script title
      paste0('cd "',aln_dir,'"'),  # change to alignment directory
      paste0(  # Gblocks command:
        "Gblocks ", aln_file,
        " -t=d -b2=",ceiling((length(cons_seqs)/2)+1)," -b4=5 -b5=n")
      # '-b4=5' allows for smaller final blocks (default=10)
      # '-b5=h' allows gap positions in final blocks (default=n)
      # '-b2=...' allows less strict flanking positions (default=85%)
    ),
    # ...to the following file:
    paste0(bscr_dir,"gblocks.sh")
  )
  if (run) {  # if specified, run the script via WSL:
    system(paste0('bash -c "',bscr_dir,'gblocks.sh"'))
  }
}