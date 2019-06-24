#' update_metadata
#'
#' Updates the metadata in from an MetagenomeSeq experiments present in expvec from an Excel spreadsheet with phenotable and phenolabels.
#' @export

update_metadata<-function(metadataXL=NULL, phenotable=NULL, phenolabels=NULL, expvec=NULL, maxclass=10, minGbNAHS=0, minPctAssembly=0, readdata=NULL){

    if(!(is.null(metadataXL))){
        metadata<-load_metadata_from_xl(xlsxFile=metadataXL)
        phenotable<-metadata[[1]]
        phenolabels<-metadata[[2]]
    }

    validcols<-as.character(phenolabels$Var_label)
    print(paste("Phenolabels contains", length(validcols), "categories."))

    #Stop if you are asking for more than you have.
    if(length(validcols) > ncol(phenotable)){
        stop("Phenolabels has more categories than the metadata in phenotable. Review and try again.")
    }

    ptsampcol<-as.character(phenolabels$Var_label[which(phenolabels$Var_type=="Sample")])
    if(length(ptsampcol)!=1){
        stop("You must define exactly one column as being the Sample column in your metadata. Please see documenation.")
    } else {
        print(paste("Samples are in the", ptsampcol, "column in the metadata."))
    }

    print(paste("Metadata classes specified in phenolabels are:", paste0(validcols, collapse=", ")))
    print("Adjusting metadata to contain only these columns.")

    #Stop if the classes you want in phenolabels do not exist in phenotable.
    if(all(validcols %in% colnames(phenotable))){
        phenotable<-phenotable[ , validcols]
    } else {
        stop("One or more classes specified in phenolabels is missing as (a) column(s) in phenotable. Please check files.")
    }

    #Check if samples from old metadata are the same as the new metadata, else stop.
    oldphenotable<-pData(expvec[[1]])
    Samples_have<-rownames(oldphenotable)
    if(!(all(Samples_have %in% phenotable[, ptsampcol]))){
        stop("Cannot drop samples which are already in experiments. Check your new metadata.")
    } else {
        #Redefine phenotable to include only samples I have in list.data
        phenotable<-phenotable[(phenotable[, ptsampcol] %in% Samples_have), ]
    }

    if(!(is.null(readdata))){
        #Make sure that pheno contains only the samples there is data for.
        #phenotable<-phenotable[readdata$Sample, ]
        #Add sequencing info
        newmetadata <- adjust_phenotable(phenotable=phenotable, phenolabels=phenolabels, readdata=readdata)
        #Get rid of samples with not enough data
        pheno<-newmetadata[[1]]
        phenolabels<-newmetadata[[2]]
        pheno<-subset(pheno, GbNAHS > minGbNAHS)
        pheno<-subset(pheno, PctAss > minPctAssembly)
    } else {
        pheno<-adjust_phenotable(phenotable=phenotable, phenolabels=phenolabels)
    }

    for(e in 1:length(expvec)){
        pData(expvec[[e]])<-pheno
    }

    print("Do not forget to re-generate variable_list and the colour dictionary, if you are using any of these.")

    return(expvec)
}
