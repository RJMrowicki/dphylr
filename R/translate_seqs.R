#' Translate DNA sequences
#'
#' Translate protein-coding DNA sequences and determine which is the valid
#' reading frame (i.e. resulting in no stop codons).
#'
#' @param seqs_obj `Biostrings` DNAStringSet object containing the sequence(s)
#'   to translate.
#' @return A list containing the following elements: `aa_seqs`, translated amino
#'   acid sequences in all three forward reading frames for each DNA sequence;
#'   `use_rfs`, the identity(ies) of valid reading frame(s) for each DNA
#'   sequence; and `no_rfs`, the number of valid reading frames for each DNA
#'   sequence.
#' @examples
#' translate_seqs()
#' @export
translate_seqs <- function (seqs_obj) {
  aa_seqs <-  # create object containing translated AA sequences
    as.list(seqs_obj) %>%  # convert DNAStringSet object to a list
    purrr::map(., ~ {  # for each list element (= sequence),
      seq_use <- .  # extract sequence to be translated
      purrr::map(1:3, ~ {  # for each codon position (= reading frame),
        suppressWarnings(  # (suppress warnings when length != 3)
          Biostrings::translate(  # (NB -- not `seqinr::translate()`)
            Biostrings::subseq(seq_use, start = .),  # translate sub-sequence
            if.fuzzy.codon = "solve"  # solve non-ambiguous fuzzy codons
          )
        )
      }) %>% Biostrings::AAStringSet()  # convert result to AAStringSet object
    })
  
  valid_rfs <-  # which reading frame(s) contain no stop codons?
    purrr::map(aa_seqs, ~ { BiocGenerics::grepl("\\*", .) == FALSE })
  
  valid_aa_seqs <-  # create object containing valid AA sequences only
    purrr::map2(aa_seqs, valid_rfs, ~ {  # for each set of AA sequences,
      # only if the corresponding no. of valid reading frames = 1,
      if (length(which(.y == TRUE)) == 1) {
        # extract the relevant AA sequence & convert to AAString object:
        .x[.y] %>% as.character() %>% Biostrings::AAString()
      } else { NULL }
    }) %>% unlist %>% Biostrings::AAStringSet()  # convert result to AAStringSet
  
  no_rfs <-  # create list containing no.s of valid reading frames (diagnostic)
    purrr::map(valid_rfs, ~ { length(which(. == TRUE)) })
  
  # output list of AA sequences and reading frame information:
  return(tibble::lst(aa_seqs, valid_rfs, valid_aa_seqs, no_rfs))
}