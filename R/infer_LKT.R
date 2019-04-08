#Infer the Last Known Taxon (LKT)
#' infer_LKT
#'
#' Infers the Last Known Taxon (LKT) from a taxonomic table or dataframe containing taxonomic information and adds an LKT column to the dataframe. 
#' @export

infer_LKT<-function(taxtable){
    dunno<-c("Unknown", "Unclassified")
    taxtable$LKT<-rep("Unclassified", length(taxtable$Domain))
    for (n in 1:length(taxtable$Domain)){
        #Start by Genus
        if (taxtable[n, which(colnames(taxtable)=="Genus")] %in% dunno){
            #Define by Family
            if (taxtable[n, which(colnames(taxtable)=="Family")] %in% dunno){
                #Define by Order
                if (taxtable[n, which(colnames(taxtable)=="Order")] %in% dunno){
                    #Define by Class
                    if (taxtable[n, which(colnames(taxtable)=="Class")] %in% dunno){
                        #Define by Phylum
                        if (taxtable[n, which(colnames(taxtable)=="Phylum")] %in% dunno){
                            #Define by Domain
                            if (taxtable[n, which(colnames(taxtable)=="Domain")] %in% dunno){
                                #nothing can be done at this point...
                            }else{
                                #LKT is the Domain
                                taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Domain")]
                            }
                        }else{
                            #LKT is the Phylum
                            taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Phylum")]
                        }
                    }else{
                         #LKT is the Class
                        taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Class")]
                    }   
                }else{
                    #LKT is the Order
                    taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Order")]
                }
            }else{
                #LKT is the Family
                taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Family")]
            }
        }else{
            taxtable[n, which(colnames(taxtable)=="LKT")]<-paste((gsub("g__", "", (taxtable[n, which(colnames(taxtable)=="Genus")]))), taxtable[n, which(colnames(taxtable)=="Species")], sep="_")
        }
    }

    taxtable$LKT<-paste("LKT", taxtable$LKT, sep="__")
    taxtable$Species<-NULL

    return(taxtable)
}




 
