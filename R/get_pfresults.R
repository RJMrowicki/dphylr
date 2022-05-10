#' Import PartitionFinder results
#'
#' Import results from a PartitionFinder analysis and create partition
#' definitions for RAxML-NG and MrBayes based on the best partitioning scheme.
#'
#' @param pf_dir Path to the root PartitionFinder directory, in which is found
#'   the 'analysis' sub-directory containing the 'best_scheme.txt' results file.
#' @return A list containing the following elements: `pf_results`, a table of
#'   the main PartitionFinder results; `pf_raxmlngparts` and `pf_mbparts`,
#'   vectors of lines specifying partition definitions for subsequent
#'   phylogenetic analyses using RAxML-NG and MrBayes, respectively; and
#'   `pf_aln_sets`, a vector of lines defining character sets to be appended to
#'   an alignment NEXUS file (used in subsequent analysis via MrBayes).
#' @examples
#' get_pfresults()
#' @export
get_pfresults <- function (pf_dir) {
  pf_out <-  # import raw PartitionFinder results
    readr::read_lines(  # (NB -- skip empty rows)
      paste0(pf_dir,"analysis/best_scheme.txt"), skip_empty_rows = TRUE
    )

  # ~ main results table:
  pf_nsubsets <-  # create object to store number of PF subsets
    stringr::str_subset(pf_out, "Number of subsets") %>%  # extract subsets line
    stringr::str_extract("(\\d)+") %>%  # extract the value itself (>=1 digit)
    as.numeric  # convert from character to numeric
  pf_out_tab0 <-  # create object to store results table header line number
    grep("Subset \\| Best Model \\|", pf_out)  # determine the line number

  pf_results <-  # create table containing main results
    pf_out[pf_out_tab0:(pf_out_tab0 + pf_nsubsets)] %>%  # extract the lines
    I %>%  # inhibit object conversion (see `readr::read_delim()`)
    readr::read_delim(delim = "|", trim_ws = TRUE) %>%  # treat as delimited
    # change PartitionFinder model names to RAxML-NG model names (see
    # <https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model>
    # and <http://www.robertlanfear.com/partitionfinder/faq/#toc-beast>):
    dplyr::mutate(dplyr::across(  # (NB -- use `across()`, not `mutate_at()`)
      "Best Model", ~ {  # replace strings within 'Best Model' column
        str_replace_all(., c(
          "T[Rr]N" = "TN93", "T[Rr]Nef" = "TN93ef", "TIM" = "TIM1"
        ))  # "TIMef" = "TIMuf" ???
      })) %>%
    # remove comma+space between multiple partition names:
    dplyr::mutate(dplyr::across(
      "Partition names", ~ {  # replace strings within 'Partition names' column
        stringr::str_replace_all(., ", ", "_")
      }))

  # ~ RAxML-NG partition definitions:
  pf_out_raxmlparts1 <-  # create object to store RAxML def.s start line number
    grep("RaxML-style partition definitions", pf_out) + 2
  pf_out_raxmlparts <-  # create vector of RAxML partition definitions
    # extract the relevant lines:
    pf_out[pf_out_raxmlparts1:(pf_out_raxmlparts1 + pf_nsubsets-1)] %>%
    # for each string, remove everything up to " = " (retain codon positions):
    purrr::map(., ~ { stringr::str_remove(., "^[^=]* = ") }) %>%
    purrr::flatten_chr()  # remove list hierarchy (convert into vector)

  pf_raxmlngparts <-  # create vector of named RAxML-NG partition definitions
    paste0(  # for each partition, paste "model, part_name = positions"
      pf_results$'Best Model', ", ",
      pf_results$'Partition names', " = ",
      stringr::str_remove_all(pf_out_raxmlparts, " ")  # remove codon spaces
    )

  # ~ MrBayes partition definitions:
  pf_out_mbparts1 <-  # create object to store MrBayes def.s start line number
    grep("MrBayes block for partition definitions", pf_out) + 3
  pf_out_mbparts <-  # create vector of MrBayes partition definitions
    # extract the relevant lines:
    pf_out[pf_out_mbparts1:(pf_out_mbparts1 + pf_nsubsets-1)] %>%
    # for each string, remove everything up to " = " (retain codon positions):
    purrr::map(., ~ { stringr::str_remove(., "^[^=]* = ") }) %>%
    purrr::flatten_chr()  # remove list hierarchy (convert into vector)

  pf_mbparts <-  # create vector of named MrBayes partition definitions
    paste0(  # for each partition, paste "  CHARSET part_name = positions;"
      "  CHARSET ", pf_results$'Partition names', " = ", pf_out_mbparts
    )
  chunks <-  # create vector of 'chunk' definitions
    purrr::map2_chr(  # for each subset number AND partition name,
      seq_len(pf_nsubsets), pf_results$'Partition names', ~ {
        paste0("CHUNK", .x, ":", .y)  # paste "CHUNK<number>:<part_name>"
      })

  pf_aln_sets <-  # create vector of NEXUS 'sets' block lines
    c(
      "BEGIN SETS;",
      pf_mbparts,  # MrBayes 'charset' lines
      # 'charpartition' line specifying 'chunks':
      paste0("  CHARPARTITION BYGENECODON = ",paste(chunks, collapse = ", "),";"),
      "END;"
    )

  # return list containing results table, RAxML-NG & MrBayes partition def.s,
  # and NEXUS 'sets' block lines:
  return(tibble::lst(pf_results, pf_raxmlngparts, pf_mbparts, pf_aln_sets))
}
