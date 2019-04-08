#' fix_kraken(taxtable)
#'
#' Fixes the mpa style output of a kaken.labels file. To be used only by the YAMSintegrate.R script.
#' @export

fix_kraken<-function(taxtable=NULL, keepstrain=FALSE){
    #Make sure we have 9 columns.
    if(ncol(taxtable) < 9){
        colstoadd<-paste0("V", (ncol(taxtable)+1):9)
        taxtable[colstoadd]<-NA
    }
    colnames(taxtable)<-c("Sequence","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    #make sure columns are vectors and not factors
    taxtable[] <- lapply(taxtable, as.vector)
    tags<-c("Sequence","^d__","^k__","^p__","^c__","^o__","^f__","^g__","^s__")
    for(l in 9:2){
        iswanted<-as.matrix(apply(taxtable, 2, function(x) grepl(tags[l], x)))
        if(nrow(taxtable)<2){
            iswanted<-t(iswanted)
        }
        for (c in (l-1):2){
            taxtable[which(iswanted[,c]==TRUE), l]<-taxtable[which(iswanted[,c]==TRUE), c]
            if (l != 2){
                    taxtable[which(iswanted[,c]==TRUE), c]<-"Unclassified"                
            }
        }
    }
    #Fill in blanks if they exist
    for (n in 2:ncol(taxtable)){
        taxtable[, n]<-sub("^$", "Unclassified", taxtable[, n])
    }
    taxtable[is.na(taxtable)] <- "Unclassified"
    taxtable$Domain<-sub("root", "Unclassified", taxtable$Domain)
    species<-sapply(strsplit(taxtable$Species, split='_', fixed=TRUE), function(x) (x[4]))
    species[is.na(species)]<-"Unclassified"
    species<-sub("sp\\.", "Unclassified", species)

    #If missing genus from NCBI taxid, fish it out from the species binomial name
    rowsmissinggenus<-which((taxtable$Genus == "Unclassified") & (taxtable$Species != "Unclassified"))
    if(length(rowsmissinggenus)>0){
        taxtable$Genus[rowsmissinggenus]<-paste0("g__", as.character(sapply(strsplit(taxtable$Species, split='_', fixed=TRUE), function(x) (x[3]))[rowsmissinggenus]))
    }
    strain<-sapply(strsplit(taxtable$Species, split='_', fixed=TRUE), function(x) (x[5]))
    strain[is.na(strain)]<-"Unclassified"
    taxtable$Species<-species
    taxtable$Kingdom<-NULL

    #Get rid of pesky characters
    taxtable[] <- lapply(taxtable, gsub, pattern='\\[', replacement='')
    taxtable[] <- lapply(taxtable, gsub, pattern='\\]', replacement='')
    taxtable[] <- lapply(taxtable, gsub, pattern='\'', replacement='')
    
    taxtable<-infer_LKT(taxtable)
    if(keepstrain==TRUE){
        taxtable$Strain<-strain
    }

   #Tag unclassifieds with appropriate level
    for (n in 2:7){
        taxtags<-c("d__","p__","c__","o__","f__","g__")
        taxtable[(which(taxtable[,n]=="Unclassified")), n]<-paste0(taxtags[n-1], "Unclassified")
    }

    return(taxtable)
}
