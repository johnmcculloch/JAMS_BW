#' make_metagenomeSeq_experiments_bytaxon(pheno=NULL, list.data=NULL, restricttotaxa=NULL)
#'
#' Returns a list vector of functional MetagenomeSeq MRexperiments for every analysis that is possible to make, restricted only to certain taxa. 
#' @export

make_metagenomeSeq_experiments_bytaxon<-function(pheno=NULL, list.data=NULL, restricttotaxa=NULL, taxlevel="LKT", onlyanalyses=NULL){
 
    #Get data for features
    Samples<-rownames(pheno)
    featureobjects<-paste(Samples, "featuredose", sep="_")
    featuredoses<-list.data[featureobjects]
    names(featuredoses)<-Samples

    #Find out which functional analyses can be made
    anallist<-lapply(1:length(featuredoses), function (x) { as.character(unique(featuredoses[[x]][]$Analysis)) })
    possibleanalyses<-Reduce(intersect, anallist)

    if(!is.null(onlyanalyses)){
        possibleanalyses<-possibleanalyses[possibleanalyses %in% onlyanalyses]
    }

    #If restricttotaxa is NOT a list of LKTs, then get LKTs according to taxlevel wanted.
    if(taxlevel != "LKT"){
        print(paste("Getting a taxtable to look up LKTs for", paste0(restricttotaxa, collapse=", ")))

        LKTobjects<-paste(Samples, "LKTdose", sep="_")
        LKTdoses<-list.data[LKTobjects]
        names(LKTdoses)<-Samples
        #Data should be redundant. But if for some reason it is not, make LKTs unique in a safe manner.
        for(s in 1:length(Samples)){
            #get rid of bogus LKT dupes due to faulty NCBI taxonomy.
            LKTdoses[[s]]<-LKTdoses[[s]][!(duplicated(LKTdoses[[s]][]$LKT)), ]
            LKTdoses[[s]][]$LKT<-as.character(LKTdoses[[s]][]$LKT)
        }
        LKTdosesall<-bind_rows(LKTdoses, .id = "id")
        colnames(LKTdosesall)[which(colnames(LKTdosesall)=="id")]<-"Sample"

        #Make tax table
        taxlvlspresent<-colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]

        Ltt<-LKTdosesall[ , taxlvlspresent]
        Ltt<-Ltt[!(duplicated(Ltt)),]
        #get rid of bogus LKT dupes due to faulty NCBI taxonomy.
        Ltt<-Ltt[!(duplicated(Ltt$LKT)), ]
        rownames(Ltt)<-Ltt$LKT
        
        WantedLKTs<-unique(rownames(Ltt[which(Ltt[,taxlevel] %in% restricttotaxa) , ]))

    } else {
        WantedLKTs<-restricttotaxa
    }

    print(paste("Will build experiments only for:", paste0(WantedLKTs, collapse=", ")))

    #Make a vector for holding experiment list 
    expvectx<-NULL
    expvectx<-vector("list",length=(length(possibleanalyses)))
    e=1

    #Make only functional analyses
    for(a in 1:length(possibleanalyses)){
        analysis<-possibleanalyses[a]
        print(paste("Making", analysis, "MRexperiment for the requested taxa"))

        #subset doses to contain only the analysis wanted
        analysisdoses<-NULL
        analysisdoses<-list()
        for(fd in 1:length(featuredoses)){
            analysisofinterest<-subset(featuredoses[[fd]], Analysis == analysis)
            if(!(all(WantedLKTs %in% colnames(analysisofinterest)))){
                #Find out which taxa are missing
                missingbugs<-WantedLKTs[!(WantedLKTs %in% colnames(analysisofinterest))]
                #Add columns for that taxon containing 0
                supplementalcols<-as.data.frame(matrix(data=0, nrow=nrow(analysisofinterest), ncol=length(missingbugs)))
                colnames(supplementalcols)<-missingbugs
                rownames(supplementalcols)<-rownames(analysisofinterest)
                analysisofinterest<-cbind(analysisofinterest, supplementalcols)
            }

            analysisdoses[[fd]]<-analysisofinterest[, c("Accession", "Description", WantedLKTs)]
            names(analysisdoses)[fd]<-names(featuredoses[fd])
        }

        featureall<-bind_rows(analysisdoses, .id = "id")
        featureall[is.na(featureall)]<-0
        colnames(featureall)[1]<-"Sample"

        featureall$NumBases<-sapply(1:nrow(featureall), function(x) { sum(featureall[x , 4:ncol(featureall)]) })
        featureall<-featureall[, c("Sample", "Accession", "Description", "NumBases")]
        cts <- featureall %>% spread(Sample, NumBases)
        cts[is.na(cts)]<-0
        rownames(cts)<-cts$Accession
        cts$Accession<-NULL
        cts$Description<-NULL

        tt<-featureall[,c("Accession", "Description")]
        tt<-tt[!(duplicated(tt)), ]
        rownames(tt)<-tt$Accession

        sampleorder<-rownames(pheno)
        tt<-tt[rownames(cts), ]
        cts<-cts[, sampleorder]

        ##Create metagenomeSeq MRexperiment
        phenotypeData = AnnotatedDataFrame(pheno)
        ttdata = AnnotatedDataFrame(tt)
        mgseq = newMRexperiment(cts, phenoData=phenotypeData, featureData=ttdata)

        attr(mgseq, "analysis")<-analysis

        expvectx[[e]]<-mgseq
        names(expvectx)[e]<-analysis
        e = e + 1
    }

    return(expvectx)
}
