#' obtain_dada2_amplicon_sizes(mergers_list = NULL, sample_name = NULL, fun_to_apply = "median")
#'
#' This is a JAMS16 function.
#' This function will return the statistical transformation set by fun_to_apply on a given sample's merged F-R pairs within a DADA2 "mergers" object.
#'
#' @export

obtain_dada2_amplicon_sizes <- function(mergers_list = NULL, sample_name = NULL, fun_to_apply = "median"){

    seqs <- mergers_list[[sample_name]]$sequence
    lengthstat <- get(fun_to_apply)(unname(sapply(seqs, function (x) { nchar(x) })))

    return(lengthstat)
}
