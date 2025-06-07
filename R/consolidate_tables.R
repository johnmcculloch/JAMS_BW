#' consolidate_tables(opt=opt)
#'
#' JAMSalpha function
#' @export

consolidate_tables <- function(opt = opt, elements = NULL, taxonomic_spaces = c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin")){


    valid_taxonomic_spaces <- taxonomic_spaces[taxonomic_spaces %in% colnames(opt$featuredata)]
    #Consolidate resfinder table
    if ("resfinder" %in% colnames(opt$featuredata)){
       resfinder_cons <- subset(opt$featuredata, resfinder != "none")
       if(nrow(resfinder_cons) > 1){
           resfinder_cons <- resfinder_cons[, c("Feature", "resfinder", "Product", valid_taxonomic_spaces)]
           opt$resfinder_cons <- resfinder_cons
       }
    }

    if ("abricate" %in% colnames(opt$featuredata)){
       abricate_cons <- subset(opt$featuredata, abricate != "none")
       if(nrow(abricate_cons) > 1){
           abricate_cons <- abricate_cons[, c("Feature", "abricate", "Product", valid_taxonomic_spaces)]
           opt$abricate_cons <- abricate_cons
       }
    }

    if ("plasmidfinder" %in% colnames(opt$featuredata)){
       plasmidfinder_cons <- subset(opt$featuredata, plasmidfinder != "none")
       if(nrow(plasmidfinder_cons) > 1){
           plasmidfinder_cons<-plasmidfinder_cons[, c("Feature", "plasmidfinder", "Product", valid_taxonomic_spaces)]
           opt$plasmidfinder_cons<-plasmidfinder_cons
       }
    }

    return(opt)

}
