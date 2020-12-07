#' Shift a codon mask
#'
#' Shift a codon mask according to a Gblocks mask
#'
#' @param cmask An `IRanges` NormalRanges object specifying the positions of
#'   either codon 1, 2 or 3 in a multiple sequence alignment.
#' @param gbmask An `IRanges` NormalRanges object specifying the start and end
#'   position(s) of Gblocks 'chunks' in a multiple sequence alignment (from
#'   `get_gbmask`).
#' @return An `IRanges` NormalRanges object specifying codon positions relative
#'   to their position within Gblocks chunks.
#' @examples
#' shift_cmask()
#' @export
shift_cmask <- function (cmask, gbmask) {
  # extract Gblocks mask start position(s) and width(s):
  gbmask_starts <- BiocGenerics::start(gbmask)
  gbmask_widths <- c(0, BiocGenerics::width(gbmask)) + 1  # first value 1
  
  c_gb_shifts <-  # create list of results for each Gblocks "chunk"
    purrr::map(seq_along(gbmask), ~ {  # for each chunk (= `IRanges` range),
      # subset codon mask by Gblocks mask:
      c_gb <- IRanges::subsetByOverlaps(cmask, gbmask[.])
      # shift subsetted codon mask according to (1) position within Gblocks chunk
      # and (2) width (+ 1) of previous Gblocks chunk (if multiple chunks)
      IRanges::shift(c_gb, (-gbmask_starts[.] + gbmask_widths[.]))
    })
  
  # combine into single `IRanges` IRangesList object and return result:
  unlist(IRanges::IRangesList(c_gb_shifts))
}