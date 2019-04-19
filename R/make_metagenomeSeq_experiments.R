#' make_metagenomeSeq_experiments(pheno=NULL, list.data=NULL, onlyanalyses=NULL, )
#'
#' Returns a list vector of MetagenomeSeq MRexperiments for every analysis that is possible to make.
#' @export

make_metagenomeSeq_experiments<-function(pheno=NULL, list.data=NULL, onlyanalyses=NULL, minnumsampanalysis=NULL, minpropsampanalysis=0.1){

    #Get data for features
    Samples <- rownames(pheno)
    featureobjects <- paste(Samples, "featuredose", sep="_")
    featuredoses <- list.data[featureobjects]
    names(featuredoses) <- Samples

    #Find out which functional analyses can be made
    anallist <- lapply(1:length(featuredoses), function (x) { as.character(unique(featuredoses[[x]][]$Analysis)) })
    names(anallist) <- Samples

    allanalyses <- Reduce(union, anallist)
    numsampwithanalysis <- NULL
    for (pa in 1:length(allanalyses)){
        possibleanalysis<-possibleanalyses[pa]
        numsampwithanalysis[pa] <- length(which(sapply(1:length(anallist), function(x) { (possibleanalysis %in% anallist[[x]]) } ) == TRUE))
    }

    names(numsampwithanalysis) <- allanalyses
    propsampleswithanalysis <- numsampwithanalysis / length(anallist)

    if (!(is.null(minpropsampanalysis))) {
        possibleanalyses <- names(propsampleswithanalysis[propsampleswithanalysis >= minpropsampanalysis])
    }

    if (!(is.null(minnumsampanalysis))) {
        possibleanalyses <- names(numsampwithanalysis[numsampwithanalysis >= minnumsampanalysis])
    }

    if(!is.null(onlyanalyses)){
        possibleanalyses <- allanalyses[allanalyses %in% onlyanalyses]
    }

    #Stop if unreasonable
    if (length(possibleanalyses) < 1){
        stop("There are no analyses fitting the criteria to make.")
    }

    #Make a vector for holding experiment list
    expvec <- NULL
    expvec <- vector("list",length=(length(possibleanalyses) + 1))
    e=1

    #Start by making LKT experiment
    print("Making LKT MRexperiment")
    LKTobjects <- paste(Samples, "LKTdose", sep="_")
    LKTdoses <- list.data[LKTobjects]
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
    taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]

    tt<-LKTdosesall[ , taxlvlspresent]
    tt<-tt[!(duplicated(tt)),]
    #get rid of bogus LKT dupes due to faulty NCBI taxonomy.
    tt<-tt[!(duplicated(tt$LKT)), ]
    rownames(tt)<-tt$LKT
    ttm<-as.matrix(tt)
    taxtypes<-sapply(1:nrow(ttm), function(x) { paste(ttm[x,], collapse="-") })

    #Ensure that LKT doses conform to the taxtable
    for(s in 1:length(Samples)){
        intt<-NULL
        taxinfo<-NULL
        taxinfo<-as.matrix(LKTdoses[[s]][ , taxlvlspresent])
        taxtypessample<-sapply(1:nrow(taxinfo), function(x) { paste(taxinfo[x,], collapse="-") })
        #Which ones are represented
        intt<-which(taxtypessample %in% taxtypes)
        LKTdoses[[s]]<-LKTdoses[[s]][intt, ]
    }

    #bind them up again
    LKTdosesall<-bind_rows(LKTdoses, .id = "id")
    colnames(LKTdosesall)[which(colnames(LKTdosesall)=="id")]<-"Sample"

    #Make counts table
    LKTallcounts<-LKTdosesall[, c("Sample", "LKT", "NumBases")]
    #Be sure that data is not empty or redundant
    LKTallcounts[is.na(LKTallcounts)]<-0
    cts <- LKTallcounts %>% spread(Sample, NumBases)
    cts[is.na(cts)]<-0
    rownames(cts)<-cts$LKT
    cts$LKT<-NULL

    sampleorder<-rownames(pheno)
    tt<-tt[rownames(cts), ]
    cts<-cts[, sampleorder]

    ##Create metagenomeSeq MRexperiment
    phenotypeData = AnnotatedDataFrame(pheno)
    ttdata = AnnotatedDataFrame(tt)
    mgseq = newMRexperiment(cts, phenoData=phenotypeData, featureData=ttdata)

    attr(mgseq, "analysis")<-"LKT"

    expvec[[e]]<-mgseq
    names(expvec)[e]<-"LKT"
    e = e + 1

    #Now, for the functional analyses
    for (a in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[a]
        print(paste("Making", analysis, "MRexperiment"))

        #Get names of Samples which have data for the analysis
        SamplesofInterest <- names(anallist)[which(sapply(1:length(anallist), function(x) { (analysis %in% anallist[[x]]) } ) == TRUE)]
        for (ad in SamplesofInterest){
            analysisdoses[[ad]] <- subset(featuredoses[[ad]], Analysis == analysis)[, c("Accession", "NumBases", "Description")]
        }
        #Take into account that the phenotable might contain less samples
        phenoanal <- pheno[SamplesofInterest, ]

        #subset doses to contain only the analysis wanted
        analysisdoses <- NULL
        analysisdoses <- list()

        featureall <- bind_rows(analysisdoses, .id = "id")
        featureall[is.na(featureall)] <- 0
        colnames(featureall)[1] <- "Sample"

        cts <- featureall %>% spread(Sample, NumBases)
        cts[is.na(cts)] <- 0
        #Just double check that there are no dupes
        if (length(which(duplicated(cts$Accession) == TRUE)) > 0){
            print(paste("Found", length(which(duplicated(cts$Accession) == TRUE)), "duplicated accessions. Keeping the first on in each case."))
            cts <- cts[!(duplicated(cts$Accession)), ]
        }
        rownames(cts) <- cts$Accession
        cts$Accession <- NULL
        cts$Description <- NULL

        ftt <- featureall[, c("Accession", "Description")]
        ftt <- ftt[!(duplicated(ftt$Accession)), ]
        rownames(ftt) <- ftt$Accession

        if (analysis == "resfinder"){
            data(resfinderlookup)
            ftt <- left_join(ftt, resfinderlookup)
            ftt[] <- lapply(ftt, as.character)
            #Fill up reasonably with data which might be empty
            ftt[is.na(ftt[])] <- "S"
            for (colm in c("Description", "Class", "ResistancePhenotype", "PMID", "Mechanism")){
                ftt[, colm][which(ftt[, colm] == "S")] <- "none"
            }
            rownames(ftt) <- ftt$Accession
            ftt$Description <- ftt$Class
            #Eliminate resfinder_none (not being a resistance gene) as a category
            cts <- cts[which(rownames(cts) != paste(analysis, "none", sep = "_")), ]
        }

        sampleorder <- rownames(phenoanal)
        ftt <- ftt[rownames(cts), ]
        cts <- cts[, sampleorder]

        ##Create metagenomeSeq MRexperiment
        phenotypeData <- AnnotatedDataFrame(phenoanal)
        ttdata <- AnnotatedDataFrame(ftt)
        mgseq <- newMRexperiment(cts, phenoData=phenotypeData, featureData=ttdata)

        attr(mgseq, "analysis") <- analysis

        expvec[[e]] <- mgseq
        names(expvec)[e] <- analysis
        e <- e + 1

        #Add an extra one if converting resfinder to antibiogram
        if (analysis == "resfinder"){
            print("Converting resfinder to antibiogram")
            expvec[[e]] <- make_antibiogram_experiment(expvec[["resfinder"]])
            names(expvec)[e] <- "antibiogram"
            e <- e + 1
        }
    }

    return(expvec)
}
