#' glom_completeness(completeness_df = NULL, taxlvl = NULL)
#'
#' JAMSalpha function
#' @export

glom_completeness <- function(completeness_df = NULL, taxlvl = NULL,  HQ_completeness_threshold = 90, HQ_contamination_threshold = 5, MHQ_completeness_threshold = 70, MHQ_contamination_threshold = 10, MQ_completeness_threshold = 50, MQ_contamination_threshold = 10, LQ_completeness_threshold = 0, LQ_contamination_threshold = 10, High_contamination_threshold = 20){
    aggregate_completeness_df <- completeness_df %>% group_by_at(taxlvl) %>% summarise(Completeness = sum(Completeness), Contamination = sum(Contamination))
    aggregate_completeness_df <- as.data.frame(aggregate_completeness_df)
    rownames(aggregate_completeness_df) <- aggregate_completeness_df[ , taxlvl]
    aggregate_completeness_df <- reformat_completeness_to_CheckM_style(completeness_df = aggregate_completeness_df)
    #Rate my bin
    aggregate_completeness_df <- rate_bin_quality(completeness_df = aggregate_completeness_df,  HQ_completeness_threshold = HQ_completeness_threshold, HQ_contamination_threshold = HQ_contamination_threshold, MHQ_completeness_threshold = MHQ_completeness_threshold, MHQ_contamination_threshold = MHQ_contamination_threshold, MQ_completeness_threshold = MQ_completeness_threshold, MQ_contamination_threshold = MQ_contamination_threshold, LQ_completeness_threshold = LQ_completeness_threshold, LQ_contamination_threshold = LQ_contamination_threshold, High_contamination_threshold = High_contamination_threshold)
    aggregate_completeness_df$Estimated_num_entities <- (aggregate_completeness_df$Completeness + aggregate_completeness_df$Contamination) / 100

    return(aggregate_completeness_df)
}
