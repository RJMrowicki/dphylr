#' Create a Gblocks mask for a multiple sequence alignment
#'
#' Create an object for masking multiple sequence alignments based on the start
#' and end positions ('flanks') of Gblocks 'chunks'.
#'
#' @param aln_dir Path to the directory containing the Gblocks output.
#' @param gb_file Name of the Gblocks `.html` output file.
#' @return An `IRanges` NormalRanges object containing start and end positions
#'   of Gblocks chunks.
#' @examples
#' get_gbmask()
#' @export
get_gbmask <- function (aln_dir, gb_file) {
  gb_out <-  # create object to store Gblocks flank positions output
    readr::read_lines(paste0(aln_dir, gb_file)) %>%  # read Gblocks output
    stringr::str_subset(., "Flanks")  # extract flank positions line
  
  gb_out_flanks <-  # create list to store flank start and end positions
    tibble::lst(
      # extract all patterns matching "[(>=1 digit)", i.e. start position(s):
      starts = unlist(stringr::str_extract_all(gb_out, "\\[(\\d)+")),
      # extract all patterns matching "(>=1 digit)]", i.e. end position(s):
      ends = unlist(stringr::str_extract_all(gb_out, "(\\d)+\\]"))
    ) %>%
    purrr::map(., ~ {  # for each list element,
      as.numeric(stringr::str_remove(., "\\[|\\]"))  # remove [ or ]
    })
  
  # return `IRanges` NormalRanges object representing mask, based on flanks:
  IRanges::IRanges(start = gb_out_flanks$starts, end = gb_out_flanks$ends)
}