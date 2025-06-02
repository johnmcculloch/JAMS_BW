#' reformat_completeness_to_CheckM_style(completeness_df = NULL)
#'
#' JAMSalpha function
#' @export

reformat_completeness_to_CheckM_style <- function(completeness_df = NULL){
    #I know a loop in inefficient, but keeping it safe for now.
    for (rn in 1:nrow(completeness_df)){
        #Check if Completeness isn't NA, as it may be for unclassifieds
        if (!is.na(completeness_df[rn, "Completeness"])){
            if (completeness_df[rn, "Completeness"] > 100){
                #Cap to 100 and put the remainder in contamination
                completeness_df[rn, "Contamination"] <- completeness_df[rn, "Completeness"] - 100
                completeness_df[rn, "Completeness"] <- 100
            } else {
                #No contamination, so change NA to 0, if "Contamination" is NA
                if (is.na(completeness_df[rn, "Contamination"])){
                    completeness_df[rn, "Contamination"] <- 0
                }
            }
        }
    }

    return(completeness_df)
}
