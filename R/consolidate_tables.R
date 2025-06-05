#' consolidate_tables(opt=opt)
#'
#' JAMSalpha function
#' @export

consolidate_tables <- function(opt = opt, elements = NULL){

    #Consolidate resfinder table
    if ("resfinder" %in% colnames(opt$featuredata)){
       resfinder_cons <- subset(opt$featuredata, resfinder != "none")
       if(nrow(resfinder_cons) > 1){
           resfinder_cons <- resfinder_cons[, c("Feature", "resfinder", "Product", "LKT")]
           opt$resfinder_cons <- resfinder_cons
       }
    }

    if ("abricate" %in% colnames(opt$featuredata)){
       abricate_cons <- subset(opt$featuredata, abricate != "none")
       if(nrow(abricate_cons) > 1){
           abricate_cons <- abricate_cons[, c("Feature", "abricate", "Product", "LKT")]
           opt$abricate_cons <- abricate_cons
       }
    }

    if ("plasmidfinder" %in% colnames(opt$featuredata)){
       plasmidfinder_cons <- subset(opt$featuredata, plasmidfinder != "none")
       if(nrow(plasmidfinder_cons) > 1){
           plasmidfinder_cons<-plasmidfinder_cons[, c("Feature", "plasmidfinder", "Product", "LKT")]
           opt$plasmidfinder_cons<-plasmidfinder_cons
       }
    }

    return(opt)

}
