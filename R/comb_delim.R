#' Combine species delimitation results
#'
#' Combine species-sample association results from multiple species delimitation
#' analyses into a single data frame.
#'
#' @param analyses A character vector containing one or more of "abgd", "gmyc"
#'   and "bptp" (default is all three), corresponding to the individual species
#'   delimitation analyses for which results are to be combined.
#' @param delim_dir Path to the main directory containing the output for ABGD
#'   and/or bPTP analyses, saved directly from the web server(s) in `.html`
#'   and/or `.txt` format, respectively.
#' @param abgd_file (If `analyses` includes "abgd") The name of the
#'   corresponding ABGD results file, in `.html` format.
#' @param gmyc_obj (If `analyses` includes "gmyc") An object describing
#'   species-sample associations determined from GMYC analysis, produced via
#'   `splits::spec.list()`.
#' @param bptp_file (If `analyses` includes "bptp") The name of the
#'   corresponding bPTP results file, in `.txt` format.
#' @param map_tree A rooted and ladderized phylogenetic tree, in ape `phylo`
#'   format, onto which species groups are to be mapped.
#' @param og A character vector containing the name(s) of the outgroup(s) used
#'   to root the reference tree.
#' @return A data frame containing combined species-sample associations with
#'   group associations expressed as sequential pseudo-numeric factors.
#' @examples
#' comb_delim()
#' @export
comb_delim <- function(
  analyses = c("abgd", "gmyc", "bptp"),
  delim_dir, abgd_file, gmyc_obj, bptp_file, map_tree, og
) {
  #
  if ("abgd" %in% analyses) {  # if analyses include ABGD,
    spp_abgd <-  # create data frame of ABGD group-seq. associations
      readr::read_lines(paste0(delim_dir,abgd_file)) %>%  # read ABGD `.html` results
      stringr::str_subset("Group\\[") %>%  # extract lines containing 'Groups' results
      stringr::str_replace(".*id: (.+)<br>", "\\1") %>%  # extract sequence codes
      # (`\\1` = string within [1st] brackets, i.e. between "...id: " & "<br>")
      stringr::str_split(" ") %>%   # split by spaces (list of seq.s per group)
      purrr::imap_dfr(., ~ {  # for each group, create table of seq. codes
        tibble::tibble("ABGD" = .y, "seq_code" = .x)
      })  # (`_dfr` combines separate tables into single data frame via row-binding)
  } else { spp_abgd <- NULL }  # else create empty object

  if ("gmyc" %in% analyses) {  # if analyses include GMYC,
    spp_gmyc <-  # create data frame of GMYC group-seq. associations
      gmyc_obj %>%  # group-seq. associations from `splits::spec.list()`
      magrittr::set_colnames(c("GMYC", "seq_code"))  # change column names
  } else { spp_gmyc <- NULL }

  if ("bptp" %in% analyses) {  # if analyses include bPTP,
    spp_bptp <-  # create data frame of bPTP group-seq. associations
      readr::read_lines(  # read bPTP `.txt` results (skip header & empty rows)
        paste0(delim_dir,bptp_file), skip = 1, skip_empty_rows = TRUE) %>%
      # extract lines containing sequence codes (exclude "Species..." headers):
      stringr::str_subset("Species ", negate = TRUE) %>%
      stringr::str_squish() %>%  # remove whitespace
      stringr::str_split(",") %>%  # split by commas (list of seq.s per group)
      purrr::imap_dfr(., ~ {  # for each group, create table of seq. codes
        tibble::tibble("bPTP" = .y, "seq_code" = .x)
      })  # (`_dfr` combines separate tables into single data frame via row-binding)
  } else { spp_bptp <- NULL }  # else create empty object

  # merge results from separate analyses into single data frame:
  spp_delim <-
    tibble::lst(spp_abgd, spp_gmyc, spp_bptp) %>%  # combine tables into a list
    purrr::discard(is.null) %>%  # remove null elements
    purrr::reduce(dplyr::left_join, by = "seq_code") %>%  # left join by seq. code
    dplyr::relocate("seq_code") %>%  # move seq. code column to beginning
    # re-order rows according to seq. positions on corresponding tree:
    dplyr::slice(match(spider::tiporder(map_tree), seq_code)) %>%
    # recode group numbers to ascending numerical order:
    dplyr::mutate(dplyr::across(-seq_code, ~ {  # for all columns except seq. code,
      ifelse(seq_code %in% og, 0, .) %>%  # replace outgroup value(s) with 0
        # convert into factor with sequentially numbered levels:
        dphylr::recode_seq() %>%  # (NB -- outgroup [i.e. first] level = 0)
        dplyr::na_if("0") %>%  # change outgroup value(s) to NA
        forcats::fct_drop()  # drop unused factor levels (i.e. 0)
    })) %>%
    # ensure same factor levels apply to all columns:
    dplyr::mutate(dplyr::across(-seq_code, ~ {  # for all columns except seq. code,
      forcats::fct_expand(  # add combined levels from the other factors
        ., dplyr::select(dplyr::cur_data_all(), -seq_code) %>% forcats::lvls_union()
      )
    }))

  return(spp_delim)  # output combined results data frame
}
