#' reformat_completeness_to_CheckM_style(completeness_df = NULL)
#'
#' JAMSalpha function
#' @export

reformat_completeness_to_CheckM_style <- function(completeness_df = NULL){
    #I know a loop in inefficient, but keeping it safe for now.
    for (rn in 1:nrow(completeness_df)){

        #If Completeness is NA, set to 0, as completenesses must be evaluated for being greater or smaller than 100
        if (is.na(completeness_df[rn, "Completeness"])){
            completeness_df[rn, "Completeness"] <- 0
        }

        #If Contamination is NA, set to 0, as completenesses exceeding 100 will be added to that.
        if (is.na(completeness_df[rn, "Contamination"])){
            completeness_df[rn, "Contamination"] <- 0
        }

        if (completeness_df[rn, "Completeness"] > 100){
            #Cap to 100 and put the remainder in contamination
            completeness_df[rn, "Contamination"] <-  completeness_df[rn, "Contamination"] + (completeness_df[rn, "Completeness"] - 100)
            completeness_df[rn, "Completeness"] <- 100
        } else {
            #No contamination, so change NA to 0, if "Contamination" is NA
            if (is.na(completeness_df[rn, "Contamination"])){
                completeness_df[rn, "Contamination"] <- 0
            }
        }
    }

    return(completeness_df)
}