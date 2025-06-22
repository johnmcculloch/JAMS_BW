#' rate_bin_quality(completeness_df = NULL, HQ_completeness_threshold = 90, HQ_contamination_threshold = 5, MHQ_completeness_threshold = 70, MHQ_contamination_threshold = 10, MQ_completeness_threshold = 50, MQ_contamination_threshold = 10, LQ_completeness_threshold = 0, LQ_contamination_threshold = 10, High_contamination_threshold = 20)
#' From the CheckM2 paper, https://doi.org/10.1038/s41592-023-01940-w
#' MIMAG quality bins, modified by McCulloch in JAMS2, adding Medium-High quality bins.
#'
#' JAMSalpha function
#' @export


rate_bin_quality <- function(completeness_df = NULL, HQ_completeness_threshold = 90, HQ_contamination_threshold = 5, MHQ_completeness_threshold = 70, MHQ_contamination_threshold = 10, MQ_completeness_threshold = 50, MQ_contamination_threshold = 10, LQ_completeness_threshold = 0, LQ_contamination_threshold = 10, High_contamination_threshold = 20){

    if (all(c("Completeness", "Contamination") %in% colnames(completeness_df))){
        completeness_df$Quality <- NA
        #Start off by censoring high contamination bins
        completeness_df[which(completeness_df$Contamination > High_contamination_threshold) , "Quality"] <- "Contaminated"

        #Mark High quality (HQ) bins, default is Completeness > 90 and contamination < 5
        completeness_df[which((completeness_df$Completeness > HQ_completeness_threshold) & (completeness_df$Contamination < HQ_contamination_threshold)), "Quality"] <- "HQ"

        #Mark medium-high qual (MHQ) bins, default is Completeness between 70 and 90 and contamination < 10
        completeness_df[which((is.na(completeness_df$Quality)) & (completeness_df$Completeness > MHQ_completeness_threshold) & (completeness_df$Contamination < MHQ_contamination_threshold)), "Quality"] <- "MHQ"

        #Mark Medium quality (MQ) bins, default is Completeness between 50 and 70 and contamination < 10
        completeness_df[which((is.na(completeness_df$Quality)) & (completeness_df$Completeness > MQ_completeness_threshold) & (completeness_df$Contamination < MQ_contamination_threshold)), "Quality"] <- "MQ"

        #Mark the remainder as Low Quality (LQ)
        completeness_df[which(is.na(completeness_df$Quality)), "Quality"] <- "LQ"

    } else {
        flog.warn("The completeness data frame must contain a Completeness and a Contamination column. Nothing was done.")
    }

    return(completeness_df)
}
