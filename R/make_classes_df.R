#' make_classes_df(curr_pt = NULL, compareby = NULL, wilcox_paired_by = NULL)
#'
#' This function makes a data frame with sample to class mapping for use by calculate_matrix_stats.
#' @export

make_classes_df <- function(curr_pt = NULL, compareby = NULL, wilcox_paired_by = NULL){

    #Create a standard data frame mapping samples to classes
    classesdf <- data.frame(Sample = rownames(curr_pt), cl = curr_pt[ , which(colnames(curr_pt) == compareby)])
    rownames(classesdf) <- classesdf$Sample
    discretenames <- unique(classesdf$cl)

    #See if we can pair it, if required
    if (!is.null(wilcox_paired_by) & length(discretenames == 2)){
        flog.info(paste("Checking if samples are paired by", wilcox_paired_by, "for Mann-Whitney-Wilcoxon test"))

        classesdf$wilcox_pairs <-  curr_pt[ , which(colnames(curr_pt) == wilcox_paired_by)]

        #Let's see if they actually form pairs
        if (length(unique(table(classesdf$wilcox_pairs))) != 1){
            flog.warn(paste("Unable to find exactly two samples per", wilcox_paired_by))
            flog.warn("Reverting to UNPAIRED Mann-Whitney-Wilcoxon test")
            classesdf$wilcox_pairs <- NULL
            return(classesdf)
        } else {
            if (unique(table(classesdf$wilcox_pairs)) != 2){
                flog.warn(paste("Unable to find exactly two samples per", wilcox_paired_by))
                flog.warn("Reverting to UNPAIRED Mann-Whitney-Wilcoxon test")
                classesdf$wilcox_pairs <- NULL
                return(classesdf)
            }
        }
        #Ensure correct ordering of pairs within each group
        discretenames <- unique(classesdf$cl)
        classesdf_1 <- subset(classesdf, cl == discretenames[1])
        pair_cat_order <- sort(classesdf_1[ , "wilcox_pairs"])
        classesdf_1 <- classesdf_1[base::match(pair_cat_order, classesdf_1$wilcox_pairs), ]
        classesdf_2 <- subset(classesdf, cl == discretenames[2])
        classesdf_2 <- classesdf_2[base::match(pair_cat_order, classesdf_2$wilcox_pairs), ]

        classesdf <- rbind(classesdf_1, classesdf_2)

    }

    return(classesdf)

}