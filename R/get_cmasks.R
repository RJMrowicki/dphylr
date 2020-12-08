#' Create codon masks for a multiple sequence alignment
#'
#' Create objects for masking a complete (i.e. prior to exclusion of ambiguous
#' sites via Gblocks) multiple sequence alignment according to codon position
#' within Gblocks 'chunks', and vectors describing codon positions.
#'
#' @param msa_obj A `Biostrings` DNAMultipleAlignment object containing the
#'   aligned sequences.
#' @param full_seq Name of the 'full' sequence within the multiple sequence
#'   alignment used to calculate the first codon 1 position.
#' @param gbmask An `IRanges` NormalRanges object specifying the start and end
#'   position(s) of Gblocks 'chunks' in the multiple sequence alignment (from
#'   `get_gbmask`).
#' @return A list containing the following elements: `cmasks`, an `IRanges`
#'   IRangesList of NormalRanges objects specifying the positions of codons 1, 2
#'   and 3 within the post-Gblocks alignment; `cmsas`, a list of `Biostrings`
#'   DNAMultipleAlignment objects for each codon position separately (used for
#'   plotting saturation plots); and `cpos`, a list of character vectors
#'   representing codon positions for use in subsequent codon partitioning
#'   phylogenetic analyses.
#' @examples
#' get_cmasks()
#' @export
get_cmasks <- function (msa_obj, full_seq, gbmask) {
  seq_starts <-  # create matrix of sequence start positions
    # locate the position of the first nucleotide base for each sequence:
    stringr::str_locate(msa_obj, "[A,C,G,T]") %>%
    # set rownames as sequence names:
    magrittr::set_rownames(BiocGenerics::rownames(msa_obj))
  full_seq_start <- seq_starts[full_seq, "start"]  # subset for 'full' sequence

  # the 'full' sequence may not start at position 1-3, therefore:
  if (full_seq_start > 3) {  # if it starts beyond position 3,
    # determine sequence of codon numbers before full sequence start:
    left_cod_pos <- rev(rep_len(3:1, full_seq_start)[-full_seq_start])
    c1_start <- which(left_cod_pos == 1)[1]  # extract position of first '1'
  } else {
    left_cod_pos <- NULL
    c1_start <- full_seq_start  # else use full sequence start position
  }

  # define the alignment width (total number of characters):
  msa_width <- BiocGenerics::ncol(msa_obj)
  # (NB -- Biostrings::maskedncol() for masked, nchar() for unmasked)

  # define `IRanges` NormalRanges objects representing masks,
  # based on positions of codons 1, 2 and 3 in full alignment:
  # (NB -- defined using first codon1 position;
  # otherwise these DO NOT represent true codon numbers,
  # UNLESS alignment is trimmed to start of coding region)
  mask_c1 <- IRanges::IRanges(start = seq(c1_start, msa_width, 3), width = 1)
  mask_c2 <- IRanges::IRanges(start = seq(c1_start+1, msa_width, 3), width = 1)
  mask_c1_c2 <- IRanges::IRanges(start = seq(c1_start, msa_width, 3), width = 2)
  mask_c3 <- IRanges::IRanges(start = seq(c1_start+2, msa_width, 3), width = 1)

  # apply codon position masks to alignments:
  # create copies of alignment, for each codon position/combination:
  msa_c1 <- msa_c2 <- msa_c1_c2 <- msa_c3 <- msa_obj
  # apply masks:
  Biostrings::colmask(msa_c1, invert = TRUE, append = "union") <- mask_c1
  Biostrings::colmask(msa_c2, invert = TRUE, append = "union") <- mask_c2
  Biostrings::colmask(msa_c1_c2, invert = TRUE, append = "union") <- mask_c1_c2
  Biostrings::colmask(msa_c3, invert = TRUE, append = "union") <- mask_c3
  # (NB -- use `append = "union"` to combine with existing masks)
  # combine masked alignments into list:
  cmsas <- tibble::lst(c1 = msa_c1, c2 = msa_c2, c1_c2 = msa_c1_c2, c3 = msa_c3)

  # # subset codon masks by Gblocks mask:
  # # (use these 'full' masks to define 'data blocks' in PartitionFinder,
  # # if using full (unmasked, pre-Gblocks) alignment as input)
  # mask_c1_gblocks <- IRanges::subsetByOverlaps(mask_c1, gbmask)
  # mask_c2_gblocks <- IRanges::subsetByOverlaps(mask_c2, gbmask)
  # mask_c1_c2_gblocks <- IRanges::subsetByOverlaps(mask_c1_c2, gbmask)
  # mask_c3_gblocks <- IRanges::subsetByOverlaps(mask_c3, gbmask)

  # subset codon masks by Gblocks mask,
  # and shift codon masks according to position in Gblocks 'chunks':
  # (use these 'collapsed' masks to define 'data blocks' in PartitionFinder,
  # if using masked (post-Gblocks) alignment as input)
  mask_c1_gblocks_col <- shift_cmask(mask_c1, gbmask)
  mask_c2_gblocks_col <- shift_cmask(mask_c2, gbmask)
  mask_c1_c2_gblocks_col <- shift_cmask(mask_c1_c2, gbmask)
  mask_c3_gblocks_col <- shift_cmask(mask_c3, gbmask)

  # create overall `IRanges` IRangesList object of codon masks:
  cmasks <- IRanges::IRangesList(
    c1 = mask_c1_gblocks_col, c2 = mask_c2_gblocks_col, c3 = mask_c3_gblocks_col
  )
  # (NB -- exclude combined 'c1_c2' mask)

  # extract positions of 1st, 2nd & 3rd codons in masked (Gblocks) alignment
  # (as vectors of single-position ranges), for later partitioning analyses:
  cpos <- IRanges::paste(
    BiocGenerics::start(cmasks), BiocGenerics::end(cmasks), sep = "-"
  ) %>% IRanges::as.list()  # (NB -- output as generic [not IRanges] list)

  # return list containing codon masks, masked alignments, and position vectors:
  return(tibble::lst(cmasks, cmsas, cpos))
}
