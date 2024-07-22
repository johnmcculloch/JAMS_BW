#Infer the Last Known Taxon (LKT)
#' infer_LKT_pr2
#'
#' Infers the Last Known Taxon (LKT) from a taxonomic table or dataframe containing taxonomic information and adds an LKT column to the dataframe.
#' @export

infer_LKT_pr2 <- function(taxtable){
    #pr2 has spelling error whereby Unclassified are incorrectly called Unclassifed within the IDtaxa database.
    dunno <- c("Unknown", "Unclassified", "Unclassifed")
    taxtable$LKT <- rep("Unclassified", length(taxtable$Domain))
    for (n in 1:length(taxtable$Domain)){
        #Start by Species
        if (taxtable[n, which(colnames(taxtable)=="Species")] %in% dunno){
           #Define by Genus
            if (taxtable[n, which(colnames(taxtable)=="Genus")] %in% dunno){
            #Define by Family
            if (taxtable[n, which(colnames(taxtable)=="Family")] %in% dunno){
                #Define by Order
                if (taxtable[n, which(colnames(taxtable)=="Order")] %in% dunno){
                    #Define by Class
                    if (taxtable[n, which(colnames(taxtable)=="Class")] %in% dunno){
                        #Define by Subdivision
                        if (taxtable[n, which(colnames(taxtable)=="Subdivision")] %in% dunno){
                            #Define by Division
                            if (taxtable[n, which(colnames(taxtable)=="Division")] %in% dunno){
                                #Define by Supergroup
                                if (taxtable[n, which(colnames(taxtable)=="Supergroup")] %in% dunno){ 
                                    #Define by Domain
                                    if (taxtable[n, which(colnames(taxtable)=="Domain")] %in% dunno){
                                #nothing can be done at this point...
                                    }else{
                                    #LKT is the Domain
                                        taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Domain")]
                                    }
                                    }else{
                                #LKT is the Supergroup
                                    taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Supergroup")]
                                }
                             }else{
                            #LKT is the Division
                                taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Division")]
                        }
                         }else{
                            #LKT is the Subdivision
                               taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Subdivision")]
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
                #LKT is the Genus
                taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Genus")]
            }
        }else{
            #LKT is the Species
                taxtable[n, which(colnames(taxtable)=="LKT")]<-taxtable[n, which(colnames(taxtable)=="Species")]
        }
    }

    taxtable$LKT<-paste("LKT", taxtable$LKT, sep="__")
    taxtable$Species<-NULL

    return(taxtable)
}
