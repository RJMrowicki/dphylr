#' Read DNA sequences
#'
#' Import all DNA sequences (FASTA format) contained in a specified directory
#' into a single `Biostrings` DNAStringSet object.
#'
#' @param seq_dir Path to the directory containing sequences, as individual
#'   FASTA files.
#' @return A `Biostrings` DNAStringSet object.
#' @examples
#' read_seqs()
#' @export
read_seqs <- function (seq_dir) {
  # generate vector of file names, including paths:
  fs::dir_ls(seq_dir, glob = "*.fa*") %>%
    # import each file as a DNAStringSet object:
    purrr::map(~ { Biostrings::readDNAStringSet(.)[[1]] }) %>%
    # rename list entries according to file name (= sequence code):
    set_names(purrr::map_chr(names(.), ~ {
      stringr::str_remove(., seq_dir) %>%  # remove directory path
        stringr::str_split(., "[.]") %>% map_chr(1)  # remove file extension
    } )) %>%
    # convert list to a DNAStringSet object:
    Biostrings::DNAStringSet()
}