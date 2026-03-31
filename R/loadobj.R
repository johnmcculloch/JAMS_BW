#' loadobj(objfn = NULL, never_used_columns = c("AntiFam", "CDD", "Coils", "FunFam", "Gene3D", "PANTHER", "Phobius", "ProSitePatterns", "SFLD", "SignalP_EUK", "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE", "SMART", "TIGRFAM", "MetaCyc")
#' JAMSbeta helper function used by load_jamsfiles_from_system.
#'
#' Do not fiddle with this.
#' @export

loadobj <- function(objfn = NULL, never_used_columns = c("AntiFam", "CDD", "Coils", "FunFam", "Gene3D", "PANTHER", "Phobius", "ProSitePatterns", "SFLD", "SignalP_EUK", "SignalP_GRAM_NEGATIVE", "SignalP_GRAM_POSITIVE", "SMART", "TIGRFAM", "MetaCyc")) {

    #All files within a .jams file should be rds these days, but check.
    is_rds <- tools::file_ext(objfn) == "rds"

    #Load the jams object without dying
    jamsobj <- tryCatch({
        if (is_rds) {
            readRDS(objfn)
        } else {
            # Assuming tsv for any non-rds file
            data.table::fread(data.table = FALSE, file = objfn, sep = "\t", header = TRUE, 
                              quote = "", fill = FALSE, integer64 = "numeric", 
                              logical01 = FALSE, stringsAsFactors = FALSE, nThread = 1)
        }
    }, error = function(e) {
        flog.warn(paste("Unable to read file", objfn, "Error:", e$message))
        return(NULL)
    })

    #Call it a day if there's nothing.
    if (is.null(jamsobj)) { 
        return(NULL)
    }

    #only attempt to prune columns in correct objects
    if (class(jamsobj) != "list"){
        # Ok, jamsobj is a df. Prune columns which won't ever be used 
        cols_to_keep <- setdiff(colnames(jamsobj), never_used_columns)
        if (length(cols_to_keep) < ncol(jamsobj)) {
            jamsobj <- jamsobj[ , cols_to_keep, drop = FALSE]
        }

        if ("Analysis" %in% colnames(jamsobj)) {
            jamsobj <- jamsobj[!(jamsobj$Analysis %in% never_used_columns), ]
        }
    } else {
        #This is a list object, so abundances list. Cycle through the elements and do the same.
        if ("functional" %in% names(jamsobj)){
            avail_taxspaces <- names(jamsobj$functional)
            for (taxsp in avail_taxspaces){
                if ("Analysis" %in% colnames(jamsobj[["functional"]][[taxsp]])) {
                    jamsobj[["functional"]][[taxsp]] <- jamsobj[["functional"]][[taxsp]][!(jamsobj[["functional"]][[taxsp]]$Analysis %in% never_used_columns), ]
                }
            }
        }
    }

    return(jamsobj)
}
