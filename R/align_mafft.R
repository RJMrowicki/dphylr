#' Align DNA sequences (MAFFT)
#'
#' Create a bash script for aligning DNA sequences via MAFFT.
#'
#' @param bscr_dir Path to the directory in which the MAFFT bash script will be
#'   created.
#' @param aln_dir Path to the directory containing the sequences to be aligned,
#'   as a single FASTA file.
#' @param aln_file Name of the FASTA file containing sequences to be aligned.
#' @param run If TRUE, the bash script is run via Windows Subsystem for Linux
#'   (WSL) as a native R process.
#' @return No objects are returned. A bash script named `mafft.sh` is created in
#'   the specified script directory. If run, this script generates an alignment
#'   as a FASTA file in the specified alignment directory.
#' @examples
#' align_mafft()
#' @export
align_mafft <- function (bscr_dir, aln_dir, aln_file, run = FALSE) {
  readr::write_lines(
    c(  # output the following lines...
      "#!/bin/bash", # bash script header
      "# Perform alignment of consensus sequence data (MAFFT)",  # script title
      paste0('cd "',aln_dir,'"'),  # change to alignment directory
      "rm mafft*",  # remove files from previous analysis
      # MAFFT alignment command:
      paste("mafft --localpair --maxiterate 10000", aln_file, "> mafft.fasta")
      # '--auto' automatically selects appropriate strategy
      # '--localpair' and '--maxiterate 1000' for L-INS-i (most accurate?)
      # '--globalpair' and '--maxiterate 1000' for G-INS-i
    ),
    # ...to the following file:
    paste0(bscr_dir,"mafft.sh")
  )
  if (run) {  # if specified, run the script via WSL:
    system(paste0('bash -c "',bscr_dir,'mafft.sh"'))
  }
}