#' evaluate_consolidated_Quality_split_df(consolidated_df = NULL)
#'
#' JAMSalpha function
#' @export

evaluate_consolidated_Quality_split_df <- function(consolidated_df = NULL){
    Quality_levels_pecking_order <- c("HQ", "MHQ", "MQ", "LQ", "Contaminated")
    Quality_levels_present <- Quality_levels_pecking_order[Quality_levels_pecking_order %in% unique(consolidated_df$Quality)]
    consolidated_Quality_split_df <- as.data.frame(consolidated_df %>% group_by(Quality) %>% summarise(Pct=round(sum(PPM)/10000, 2)))
    consolidated_Quality_split_df <- consolidated_Quality_split_df[match(Quality_levels_present, consolidated_Quality_split_df$Quality) , ]
    consolidated_Quality_split_df$CumPct <- cumsum(consolidated_Quality_split_df$Pct) 
    rownames(consolidated_Quality_split_df) <- consolidated_Quality_split_df$Quality

    return(consolidated_Quality_split_df)
}