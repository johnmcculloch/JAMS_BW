#' multiple_subsetting_sample_selector(SEobj = NULL, phenotable = NULL, subsetby = NULL, compareby = NULL, cats_to_ignore = "N_A")
#'
#' Application for subsetting a phenotable into up to three tiers of subsets within subsets. Returns a named list. The first element, named "Subsets_stats" is a data frame of the compareby stats within each specific subset. Subsequent list named items are vectors of the sample names within each subset specified by the "Subset_Tier_Class_Name" column of the "Subsets_stats" data frame.
#' @export

multiple_subsetting_sample_selector <- function(SEobj = NULL, phenotable = NULL, subsetby = NULL, compareby = NULL, cats_to_ignore = "N_A"){

    #Ensure subsetby has length up to 3
    subsetby <- subsetby[1:3]
    subsetby <- subsetby[!is.na(subsetby)]
    subset_level_1 <- subsetby[1]
    subset_level_2 <- subsetby[2]
    subset_level_3 <- subsetby[3]

    if (is.null(phenotable)){
        curr_pheno <- as.data.frame(colData(SEobj))
    } else {
        curr_pheno <- as.data.frame(phenotable)
    }

    variable_list <- define_kinds_of_variables(phenotable = curr_pheno, verbose = FALSE)
    #Check if subset variables are in phenotable and are discrete
    if (!any(c(subset_level_1, subset_level_2, subset_level_3) %in% variable_list$discrete)){
        flog.warn(paste("ABORTING:", paste0(c(subset_level_1, subset_level_2, subset_level_3)[!(c(subset_level_1, subset_level_2, subset_level_3) %in% variable_list$discrete)], collapse = ", "), "are not discrete variables within the SummarizedExperiment object phenotable."))

        return(NULL)
    }

    valid_subsets <- c(subset_level_1, subset_level_2, subset_level_3)[c(subset_level_1, subset_level_2, subset_level_3) %in% variable_list$discrete]
    flog.info(paste("Valid subsets are:", paste0(paste(c("Tier 1 subset =", "Tier 2 subset =", "Tier 3 subset =")[1:length(valid_subsets)], valid_subsets), collapse = "; ")))

    #Test for compareby
    if (!is.null(compareby)){
        #Test for silly stuff. Is it in pheno?
        if (!(compareby %in% colnames(curr_pheno))){
            flog.warn(paste0("Variable ", compareby, ", passed to the compareby argument was not found on phenotable. Check your metadata. Setting compareby to NULL."))
            compareby <- NULL
        }
    }

    #Curtail phenotable to contain only relevant information
    curr_pheno <- curr_pheno[ , c(valid_subsets, compareby), drop = FALSE]
    #Eliminate any unwanted cats from phenotable
    for (colm in colnames(curr_pheno)){
        curr_pheno <- curr_pheno[which(!(curr_pheno[ , colm] %in% cats_to_ignore)), , drop = FALSE]
    }

    calculate_pielou_j <- function(class_vec = NULL){
        H_prime <- vegan::diversity(table(class_vec))
        max_H <- log(length(unique(class_vec)))
        pielou_j <- H_prime / max_H

        return(pielou_j)
    }

    apply_appropriate_compareby_QC_test <- function(class_vec = NULL, class_type = NULL){
        if (!all(is.na(class_vec))) {
            if (class_type == "discrete"){
                #Test for balance (as gauged by evenness) of classes
                class_stat <- calculate_pielou_j(class_vec = class_vec)
             } else if (class_type == "continuous"){
                #Test for normalcy of continuous variable
                class_stat <- shapiro.test(as.numeric(class_vec))$statistic
            }
            #Round to 3 SF
            class_stat <- round(class_stat, 3)
        } else {
            #no info passed, avoid calculation and return NA
            class_stat <- NA
        }

        return(class_stat)
    }

    apply_appropriate_compareby_samp_num_test <- function(class_vec = NULL, class_type = NULL){
        if (!all(is.na(class_vec))) {
            if (class_type == "discrete"){
                class_samp_minmax <- c(min(table(class_vec)), max(table(class_vec)))
             } else if (class_type == "continuous"){
                #Test for normalcy of continuous variable
                class_samp_minmax <- NA
            }
       } else {
            #no info passed, avoid calculation and return NA
            class_samp_minmax <- NA
        }

        return(class_samp_minmax)
    }

    #Make data frame of non-unique combis from highest to lowest levels
    colnamedict <- data.frame(VarName = c(compareby, valid_subsets), AttribName = c("Compareby"[length(compareby)], c("Cats_lvl_1", "Cats_lvl_2", "Cats_lvl_3")[1:length(valid_subsets)]))
    tmp_pheno <- curr_pheno
    colnames(tmp_pheno) <- colnamedict$AttribName[match(colnames(tmp_pheno), colnamedict$VarName)]
    if (!("Compareby" %in% colnames(tmp_pheno))){
        tmp_pheno$Compareby <- NA
        compareby <- NA
    }

    #What kind of variable is it?
    CVT <- c("discrete", "continuous", NA)[c((compareby %in% variable_list$discrete), (compareby %in% variable_list$continuous), is.na(compareby))]
    CTN <- c("Pielou_J", "Shapiro-Wilk", NA)[c((compareby %in% variable_list$discrete), (compareby %in% variable_list$continuous), is.na(compareby))]

    #Level 0
    Subsets_stats_df <- tmp_pheno %>% summarise(Subset_Tier_Level = 0, Subset_Tier_Variable_Name = "No subset", Subset_Tier_Class_Name = NA, Num_samples_in_subset = length(Compareby), Compareby_Variable_Name = compareby, Compareby_Variable_Type = CVT, Num_cats_compareby_in_subset = length(unique(Compareby)), Compareby_QC_in_subset = apply_appropriate_compareby_QC_test(class_vec = Compareby, class_type = CVT), Compareby_Test_Name_in_subset = CTN, Compareby_min_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[1], Compareby_max_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[2]) %>% as.data.frame()

    #Level 1
    Subsets_stats_df_suppl <- tmp_pheno %>% group_by(Cats_lvl_1) %>% summarise(Subset_Tier_Level = 1, Subset_Tier_Variable_Name = valid_subsets[1], Subset_Tier_Class_Name = NA, Num_samples_in_subset = length(Compareby), Compareby_Variable_Name = compareby, Compareby_Variable_Type = CVT, Num_cats_compareby_in_subset = length(unique(Compareby)),  Compareby_QC_in_subset = apply_appropriate_compareby_QC_test(class_vec = Compareby, class_type = CVT), Compareby_Test_Name_in_subset = CTN, Compareby_min_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[1], Compareby_max_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[2]) %>% as.data.frame()
    colnames(Subsets_stats_df_suppl)[which(colnames(Subsets_stats_df_suppl) == "Cats_lvl_1")] <- "Subset_Tier_Class_Name"
    Subsets_stats_df_suppl <- Subsets_stats_df_suppl[ , colnames(Subsets_stats_df)] 
    Subsets_stats_df <- rbind(Subsets_stats_df, Subsets_stats_df_suppl)

    #Level 2
    if (length(valid_subsets) > 1) {
        Subsets_stats_df_suppl <- tmp_pheno %>% group_by(Cats_lvl_1, Cats_lvl_2) %>% summarise(Subset_Tier_Level = 2, Subset_Tier_Variable_Name = valid_subsets[2], Num_samples_in_subset = length(Compareby), Compareby_Variable_Name = compareby, Compareby_Variable_Type = CVT, Num_cats_compareby_in_subset = length(unique(Compareby)),  Compareby_QC_in_subset = apply_appropriate_compareby_QC_test(class_vec = Compareby, class_type = CVT), Compareby_Test_Name_in_subset = CTN, Compareby_min_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[1], Compareby_max_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[2]) %>% as.data.frame()
        Subsets_stats_df_suppl$Subset_Tier_Class_Name <- paste(Subsets_stats_df_suppl$Cats_lvl_1, Subsets_stats_df_suppl$Cats_lvl_2, sep = "§")
        Subsets_stats_df_suppl <- Subsets_stats_df_suppl[ , colnames(Subsets_stats_df)] 
        Subsets_stats_df <- rbind(Subsets_stats_df, Subsets_stats_df_suppl)
    }

    #Level 3
    if (length(valid_subsets) > 2) {
        Subsets_stats_df_suppl <- tmp_pheno %>% group_by(Cats_lvl_1, Cats_lvl_2, Cats_lvl_3) %>% summarise(Subset_Tier_Level = 3, Subset_Tier_Variable_Name = valid_subsets[3], Num_samples_in_subset = length(Compareby), Compareby_Variable_Name = compareby, Compareby_Variable_Type = CVT, Num_cats_compareby_in_subset = length(unique(Compareby)), Compareby_QC_in_subset = apply_appropriate_compareby_QC_test(class_vec = Compareby, class_type = CVT), Compareby_Test_Name_in_subset = CTN, Compareby_min_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[1], Compareby_max_num_samples_in_class = apply_appropriate_compareby_samp_num_test(class_vec = Compareby, class_type = CVT)[2]) %>% as.data.frame()
        Subsets_stats_df_suppl$Subset_Tier_Class_Name <- paste(Subsets_stats_df_suppl$Cats_lvl_1, Subsets_stats_df_suppl$Cats_lvl_2, Subsets_stats_df_suppl$Cats_lvl_3, sep = "§")
        Subsets_stats_df_suppl <- Subsets_stats_df_suppl[ , colnames(Subsets_stats_df)] 
        Subsets_stats_df <- rbind(Subsets_stats_df, Subsets_stats_df_suppl)
    }

    #Clean up
    if (!(CVT %in% "discrete")) {
        Subsets_stats_df$Compareby_min_num_samples_in_class <- NULL
        Subsets_stats_df$Compareby_max_num_samples_in_class <- NULL
    }

    sample_subset_list <- list()
    sample_subset_list[["Subsets_stats"]] <- Subsets_stats_df
    #gather vectors of sample names into list
    wanted_subsets <- Subsets_stats_df$Subset_Tier_Class_Name[!is.na(Subsets_stats_df$Subset_Tier_Class_Name)]
    for (subUID in wanted_subsets){
        wanted_pheno_rn <- 1:nrow(curr_pheno)
        split_subset_path <- unlist(strsplit(subUID, split = "§"))
        for (slvl in 1:length(split_subset_path)){
            wanted_pheno_rn <- intersect(wanted_pheno_rn, which(curr_pheno[ , valid_subsets[slvl]] == split_subset_path[slvl]))
            sample_subset_list[[subUID]] <- rownames(curr_pheno[wanted_pheno_rn, , drop = FALSE])
        }
    }

    return(sample_subset_list)
} 

