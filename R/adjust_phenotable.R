#Prepare pheno object for metagenomeSeq package
#' adjust_phenotable(phenotable=NULL, phenolabels=NULL, readdata=NULL, list.data=list.data)
#'
#' Adjusts the metadata table into a format which can be used by the MetagenomeSeq package.
#' @export 

adjust_phenotable<-function(phenotable=NULL, phenolabels=NULL, readdata=NULL, list.data=NULL, addtaxlevelstoisolates=NULL){

    print("Adjusting phenotable classes by type of variable.")
    varlist<-define_kinds_of_variables(phenolabels=phenolabels, phenotable=phenotable, maxclass=10000, maxsubclass=10000)

    pheno <- phenotable
    #Adjust discrete and subsettable variables as class character
    charvar <- unique(c(varlist$sample, varlist$discrete, varlist$subsettable)) 
    for(cv in 1: length(charvar)){
        pheno[,charvar[cv]]<-as.character(pheno[,charvar[cv]])
    }
    #Adjust continuous variables as class numeric
    numvar <- unique(varlist$continuous)
    if(length(numvar) > 0){
        for(n in 1:length(numvar)){
            pheno[,numvar[n]]<-as.numeric(pheno[,numvar[n]])
        }
    }
    rownames(pheno) <- pheno[,varlist$sample]
    Samples<-rownames(pheno)

    #add information regarding sample type
    if(!is.null(list.data)){
        projdata<-as.data.frame(matrix(data="unknown", nrow=(length(Samples)), ncol=3))
        projdata[] <- lapply(projdata, as.character)
        colnames(projdata)<-c("Sample", "JAMS_Run_type", "JAMS_Process")
        projdata$Sample<-Samples

        #Fetch data pertaining to each sample
        for (s in 1:length(Samples)){
            projstats<-NULL
            projstats<-list.data[[paste(Samples[s], "projinfo", sep="_")]]
            rownames(projstats)<-projstats$Run_info
            projstats<-projstats[c("Run_type", "Process"), ]
            projstats[] <- lapply(projstats, as.character)
            projdata[which(projdata$Sample == Samples[s]), which(colnames(projdata)=="JAMS_Run_type")]<-projstats["Run_type", "Run_value"]
            projdata[which(projdata$Sample == Samples[s]), which(colnames(projdata)=="JAMS_Process")]<-projstats["Process", "Run_value"]
        }
    
        pheno<-left_join(pheno, projdata, by="Sample")
        phenolabels_projinfo<-data.frame(Var_label=c("JAMS_Run_type", "JAMS_Process"), Var_type=c("discrete", "discrete"), stringsAsFactors =FALSE)
        phenolabels<-rbind(phenolabels, phenolabels_projinfo)

        #add taxonomic information to isolates, if there are any
        if((!is.null(addtaxlevelstoisolates)) && length(which(projdata$JAMS_Run_type == "isolate") > 0)){
            taxdata<-as.data.frame(matrix(data="not_isolate", nrow=(length(Samples)), ncol=(length(addtaxlevelstoisolates) +1 )))
            taxdata[] <- lapply(taxdata, as.character)
            newtaxlvlnames<-paste("Isolate", addtaxlevelstoisolates, sep="_")
            colnames(taxdata)<-c("Sample", newtaxlvlnames)
            taxdata$Sample<-Samples
            isolatesamples<-taxdata$Sample[which(projdata$JAMS_Run_type == "isolate")]
            for (i in 1:length(isolatesamples)){
                taxstats<-NULL
                taxstats<-list.data[[paste(isolatesamples[i], "LKTdose", sep="_")]]
                taxstats<-taxstats[, c(addtaxlevelstoisolates, "NumBases")]
                taxstats[] <- lapply(taxstats, as.character)
                taxstats$NumBases<-as.numeric(taxstats$NumBases)
                taxdata[which(taxdata$Sample == isolatesamples[i]), newtaxlvlnames]<-as.character(as.matrix(unname(taxstats[which(taxstats$NumBases == max(taxstats$NumBases)), addtaxlevelstoisolates])))
            }

            pheno<-left_join(pheno, taxdata, by="Sample")
            phenolabels_taxon<-data.frame(Var_label=newtaxlvlnames, Var_type=rep("discrete", length(newtaxlvlnames)), stringsAsFactors =FALSE)
            phenolabels<-rbind(phenolabels, phenolabels_taxon)
        }
    }

    #Add data from readdata if desired
    if(!(is.null(readdata))){
        dfr<-readdata
        dfr$GbNAHS<-round((dfr$NonHost_bases/1000000000),2)
        dfr$GbTrim<-round((dfr$Trim_bases/1000000000),2)
        pheno$GbNAHS<-rep(0, nrow(pheno)) #Account for the fact that pheno Samples may be missing in readdata
        pheno$PctAss<-rep(0, nrow(pheno)) #Account for the fact that pheno Samples may be missing in readdata
        pheno$GbNAHS<-dfr$GbNAHS[match(pheno[,varlist$sample], dfr$Sample)]
        pheno$PctAss<-dfr$PctAss[match(pheno[,varlist$sample], dfr$Sample)]
        #phenolabels_reads<-data.frame(Var_label=c("GbNAHS", "PctAss"), Var_type=c("continuous","continuous"), stringsAsFactors =FALSE)
        #phenolabels<-rbind(phenolabels, phenolabels_reads)
    }

    metadata<-list()
    rownames(pheno) <- pheno[,varlist$sample]
    metadata[[1]]<-pheno
    metadata[[2]]<-phenolabels
    names(metadata)<-c("phenotable", "phenolabels")

    return(metadata)
}
