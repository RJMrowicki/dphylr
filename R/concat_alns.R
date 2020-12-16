#' Concatenate alignments and create alignment/codon masks
#'
#' Concatenate two or more single-marker multiple sequence alignments, which
#' have already been aligned via MAFFT and had ambiguous sites excluded via
#' Gblocks. Also, create objects for masking the concatenated alignment by
#' marker and codon position, and vectors describing marker/codon positions.
#'
#' @param aln_main_dir Path to the main directory containing all single-marker
#'   alignments to be concatenated, each contained in a separate folder named
#'   according to the marker, in FASTA format (named `mafft.fasta-gb`).
#' @param samp_seq_obj A table of sample codes (`samp_code`) and corresponding
#'   sequence codes (`seq_code`), accounting for all sequences contained in the
#'   single-marker alignments.
#' @param cpart Are the single-marker alignments to be partitioned by codon
#'   position for subsequent phylogenetic analyses?
#' @return A list containing the following elements: `aln_conc`, a `BioStrings`
#'   DNAMultipleAlignment object containing the concatenated alignment;
#'   `aln_conc_samp_codes`, a table of sample codes in the concatenated
#'   alignment and corresponding sequence codes for each marker; `alns_masks`,
#'   an `IRanges` NormalRanges object specifying the start and end positions of
#'   single-marker alignments within the concatenated alignment, also
#'   represented by `alns_pos`, a list of character vectors, for use in
#'   subsequent partitioning phylogenetic analyses; and (if partitioning by
#'   codon) `alns_cmasks`, a list of `IRanges` IRangesLists of NormalRanges
#'   objects specifying the positions of codons 1, 2 and 3 within each
#'   single-marker alignment, also represented by `alns_cpos`, a list of
#'   character vectors for use in subsequent partitioning phylogenetic analyses.
#' @examples
#' concat_alns()
#' @export
concat_alns <- function (aln_main_dir, samp_seq_obj, cpart = FALSE) {
  alns <-  # create list of **post-Gblocks** single-marker alignments
    fs::dir_ls(aln_main_dir, type = "directory") %>%  # get aln. folder paths
    stringr::str_subset("conc", negate = TRUE) %>%  # remove "conc" directory
    # name vector of folder paths (remove path, leaving folder name):
    purrr::set_names(., ~ { stringr::str_remove(., aln_main_dir) }) %>%
    {suppressWarnings(  # suppress warnings related to spaces in input
      purrr::map(., ~ {  # for each folder, import the sequence alignment
        Biostrings::readDNAStringSet(paste0(.,"/mafft.fasta-gb"), format = "fasta")
      })
    )}

  aln_samp_codes <-  # create list of corresponding alignment sample codes
    purrr::map(alns, ~ {  # for each single-marker alignment,
      # subset samp./seq. code table according to alignment sequence codes:
      dplyr::slice(samp_seq_obj, match(names(.x), seq_code)) %>%
        dplyr::select(., samp_code, seq_code)  # retain samp. and seq. codes
    })

  shared_samp_codes <-  # create vector of sample codes shared among alignments
    # (NB -- ensure none of the sample codes in `samp_seq` are 'NA')
    purrr::map(aln_samp_codes, "samp_code") %>%  # extract samp. codes from each table
    purrr::reduce(intersect)  # reduce list to vector via intersect operation

  aln_conc <-  # create concatenated alignment
    purrr::map2(aln_samp_codes, alns, ~ {  # for each single-marker alignment,
      shared_seq_codes <-  # create vector of shared sequence codes
        # subset samp./seq. code table according to shared samp. codes
        # (NB -- this preserves the sequence order when subsetting below):
        dplyr::slice(.x, match(shared_samp_codes, samp_code)) %>%
        dplyr::select(seq_code) %>% dplyr::pull()  # extract seq. codes vector
      # subset the alignment according to shared seq. codes:
      .y[shared_seq_codes]
    }) %>%
    purrr::reduce(Biostrings::xscat) %>%  # concatenate alignments
    stats::setNames(shared_samp_codes) %>%  # set names as samp. codes
    Biostrings::DNAMultipleAlignment()  # convert to DNAMultipleAlignment object

  aln_conc_samp_codes <-  # create summary table of samp. and seq. codes
    aln_samp_codes %>%  # using list of alignment samp./seq. codes,
    # join tables by samp. code, name separate seq. code columns by marker:
    purrr::reduce(full_join, by = "samp_code", suffix = paste0(".",names(.))) %>%
    # retain rows corresponding to shared sample codes only:
    dplyr::filter(., samp_code %in% shared_samp_codes)

  aln_widths <-  # create vector of alignment widths (= number of characters)
    purrr::map_int(alns, ~ { BiocGenerics::width(.)[1] })
  # calculate start and end positions from cumulative sums of aln. widths:
  aln_starts <- (cumsum(aln_widths) + 1) - aln_widths
  aln_ends <- cumsum(aln_widths)

  alns_masks <-  # create list of single-marker alignment masks
    purrr::map2(aln_starts, aln_ends, ~ {  # for each aln. start & end position,
      # define `IRanges` NormalRanges object representing alignment mask:
      IRanges::IRanges(start = .x, end = .y)
    }) %>% IRanges::IRangesList()  # convert into `IRanges` IRangesList

  # alns_masks <-  # create masks representing start and end positions of aln.s
  #   aln_widths %>% IRanges::IRanges(  # define `IRanges` NormalRanges object
  #     # start and end positions calculated from cumulative sums of aln. widths:
  #     start = (cumsum(.) + 1) - ., end = cumsum(.), names = names(.)
  #   )

  alns_pos <-  # create vector of single-marker aln. start and end positions
    IRanges::paste(
      BiocGenerics::start(alns_masks), BiocGenerics::end(alns_masks), sep = "-"
    ) %>% IRanges::as.list()  # (NB -- output generic [not IRanges] list)

  if (cpart) {  # if partitioning by codon position,
    # (NB -- what to do if combining protein-coding and non-protein-coding
    # markers? distinguish by e.g. `alns_cmasks == NULL`, below...???)

    alns_cshifts <-  # calculate corresponding 'shift' values for codon masks
      cumsum(aln_widths) - aln_widths

    alns_cmasks <-  # create list of shifted alignment codon masks
      purrr::map2(names(alns), alns_cshifts, ~ {  # for each single-marker aln.,
        # if there is an existing 'cmasks' file (i.e. marker has been
        # partitioned by codon in previous single-marker analysis),
        if (fs::file_exists(paste0(aln_main_dir, .x ,"/cmasks"))) {
          # load previously created 'codon_masks' from alignment folder:
          load(paste0(aln_main_dir, .x ,"/cmasks"))
          # 'shift' codon masks by corresponding amount:
          IRanges::shift(cmasks, .y)
        } else { NULL }  # otherwise, output NULL
      }) %>% set_names(names(alns))  # name according to markers

    alns_cpos <-  # create list of codon positions (as vectors of pos. ranges)
      purrr::map2(alns_cmasks, alns_pos,  ~ {  # for each single-marker aln.,
        if (!is.null(.x)) {  # if codon masks have been specified,
          IRanges::paste(
            BiocGenerics::start(.x), BiocGenerics::end(.x), sep = "-"
          ) %>% IRanges::as.list()  # (NB -- output generic [not IRanges] list)
        } else {  # otherwise, if codon masks have not specified,
          .y  # use marker start and end position
        }
      }) %>%
      unlist(., recursive = FALSE) %>%  # flatten list
      # replace "." in combined <marker.codon> list element names:
      purrr::set_names(., ~ { stringr::str_replace(., "\\.", "_") })

  } else {  # otherwise, partitions are defined by constituent markers only:
    alns_cmasks <- alns_cpos <-  # create empty lists
      vector("list", length(alns)*3) %>%
      purrr::set_names(  # name by all combinations of marker and codon no.
        paste(rep(names(alns), each = 3), c("c1","c2","c3"), sep = "_")
      )
  }

  # return list containing the concatenated alignment, summary table of samp.
  # and seq. codes, single-marker alignment masks, and codon masks & positions:
  return(tibble::lst(
    aln_conc, aln_conc_samp_codes, alns_masks, alns_pos, alns_cmasks, alns_cpos
  ))
}
