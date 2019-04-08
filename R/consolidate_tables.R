#' consolidate_tables(opt=opt)
#'
#' JAMSalpha function
#' @export

consolidate_tables<-function(opt=opt, elements=NULL){

    #Consolidate 16S table
    if("SixteenSid" %in% colnames(opt$featuredata)){
       taxa_16S_cons<-subset(opt$featuredata, SixteenSid != "none")
       if(nrow(taxa_16S_cons) > 1){
           taxa_16S_cons<-taxa_16S_cons[, c("Feature", "LengthDNA", "Contig", "SixteenSid")]
           taxa_16S_cons<-left_join(taxa_16S_cons, opt$contigsdata)
           taxa_16S_cons<-taxa_16S_cons[, c("Feature", "LengthDNA", "Contig", "Length", "SixteenSid", "LKT")]
           colnames(taxa_16S_cons)<-c("Gene", "Length_16S", "Contig", "Contig_Length", "LKT_dada2", "LKT_kraken")
           opt$taxa_16S_cons<-taxa_16S_cons
       }
    }

    #Consolidate resfinder table
    if("resfinder" %in% colnames(opt$featuredata)){
       resfinder_cons<-subset(opt$featuredata, resfinder != "none")
       if(nrow(resfinder_cons) > 1){
           resfinder_cons<-resfinder_cons[, c("Feature", "resfinder", "Product", "LKT")]
           opt$resfinder_cons<-resfinder_cons
       }
    }

    if("plasmidfinder" %in% colnames(opt$featuredata)){
       plasmidfinder_cons<-subset(opt$featuredata, plasmidfinder != "none")
       if(nrow(plasmidfinder_cons) > 1){
           plasmidfinder_cons<-plasmidfinder_cons[, c("Feature", "plasmidfinder", "Product", "LKT")]
           opt$plasmidfinder_cons<-plasmidfinder_cons
       }
    }

    return(opt)

}
