#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL,  onlyanalyses = NULL, minnumsampanalysis = NULL, minpropsampanalysis = 0.1, minPctFromCtg = NULL, minProbNumGenomes = NULL, restricttoLKTs = NULL, list.data = NULL)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL,  onlyanalyses = NULL, minnumsampanalysis = NULL, minpropsampanalysis = 0.1, restricttoLKTs = NULL, list.data = NULL){

    require(SummarizedExperiment)

    #Get data for features
    if (!is.null(onlysamples)){
        pheno2 <- pheno[rownames(pheno) %in% onlysamples, ]
    } else {
        pheno2 <- pheno
    }
    Samples <- rownames(pheno2)
    featureobjects <- paste(Samples, "featuredose", sep="_")
    featuredoses <- list.data[featureobjects]
    names(featuredoses) <- Samples

    #Find out which functional analyses can be made
    anallist <- lapply(1:length(featuredoses), function (x) { as.character(unique(featuredoses[[x]][]$Analysis)) } )
    names(anallist) <- Samples

    allanalyses <- Reduce(union, anallist)
    numsampwithanalysis <- NULL
    for (pa in 1:length(allanalyses)){
        possibleanalysis <- allanalyses[pa]
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

    if (!is.null(onlyanalyses)){
        possibleanalyses <- allanalyses[allanalyses %in% onlyanalyses]
    } else {
        possibleanalyses <- allanalyses
    }

    #Stop if unreasonable
    if (length(possibleanalyses) < 1){
        stop("There are no analyses fitting the criteria to make.")
    }

    #Make a vector for holding experiment list
    expvec <- NULL
    expvec <- vector("list",length = (length(possibleanalyses) + 1))
    e <- 1

    if (is.null(restricttoLKTs)){
        #Start by making LKT experiment
        flog.info("Making LKT SummarizedExperiment")
        LKTobjects <- paste(Samples, "LKTdose", sep="_")
        LKTdoses <- list.data[LKTobjects]
        names(LKTdoses) <- Samples

        #Data should be non-redundant. But if for some reason it is, make LKTs unique in a safe manner.
        for (s in 1:length(Samples)){
            #get rid of bogus LKT dupes due to faulty NCBI taxonomy.
            LKTdoses[[s]] <- LKTdoses[[s]][!(duplicated(LKTdoses[[s]][]$LKT)), ]
            LKTdoses[[s]][]$LKT <- as.character(LKTdoses[[s]][]$LKT)
        }
        LKTdosesall <- bind_rows(LKTdoses, .id = "id")
        colnames(LKTdosesall)[which(colnames(LKTdosesall) == "id")] <- "Sample"

        #Make tax table
        taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]

        tt <- LKTdosesall[ , taxlvlspresent]
        tt <- tt[!(duplicated(tt)), ]
        #get rid of LKT dupes due to NCBI taxonomy names at species level including subspecies nomenclature. Sigh. These are usually unclassified species.
        tt <- tt[!(duplicated(tt$LKT)), ]
        rownames(tt) <- tt$LKT
        ttm <- as.matrix(tt)

        #Make counts table
        LKTallcounts <- LKTdosesall[, c("Sample", "LKT", "NumBases")]
        #Be sure that data is not empty or redundant
        LKTallcounts[is.na(LKTallcounts)] <- 0

        cts <- spread(LKTallcounts, Sample, NumBases, fill = 0, drop = FALSE)
        rownames(cts) <- cts$LKT
        cts$LKT <- NULL
        cts <- cts[order(rowSums(cts), decreasing = TRUE), ]
        featureorder <- rownames(cts)

        sampleorder <- rownames(pheno2)
        tt <- tt[featureorder, ]
        cts <- cts[, sampleorder]

        #Register the total number of NAHS bases sequenced for each sample
        TotBasesSamples <- colSums(cts)
        TotalBasesSequenced <- t(as.matrix(TotBasesSamples))
        rownames(TotalBasesSequenced) <- "NumBases"

        #Get percentage from contig matrix
        LKTallPctFromCtg <- LKTdosesall[, c("Sample", "LKT", "PctFromCtg")]
        #Be sure that data is not empty or redundant
        LKTallPctFromCtg[is.na(LKTallPctFromCtg)] <- 0
        LKTPctFromCtgcts <- spread(LKTallPctFromCtg, Sample, PctFromCtg, fill = 0, drop = FALSE)
        rownames(LKTPctFromCtgcts) <- LKTPctFromCtgcts$LKT
        LKTPctFromCtgcts$LKT <- NULL
        LKTPctFromCtgcts <- LKTPctFromCtgcts[featureorder, sampleorder]

        #Get genome completeness matrix
        LKTallGenComp <- LKTdosesall[, c("Sample", "LKT", "ProbNumGenomes")]
        #Be sure that data is not empty or redundant
        LKTallGenComp[is.na(LKTallGenComp)] <- 0
        LKTallGenCompcts <- spread(LKTallGenComp, Sample, ProbNumGenomes, fill = 0, drop = FALSE)
        rownames(LKTallGenCompcts) <- LKTallGenCompcts$LKT
        LKTallGenCompcts$LKT <- NULL
        LKTallGenCompcts <- LKTallGenCompcts[featureorder, sampleorder]

        assays <- list(cts, LKTPctFromCtgcts, LKTallGenCompcts)
        names(assays) <- c("BaseCounts", "PctFromCtgs", "GenomeCompleteness")

        SEobj <- SummarizedExperiment(assays = assays, rowData = tt, colData = pheno2)
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$analysis <- analysis

        expvec[[e]] <- SEobj
        names(expvec)[e] <- "LKT"
        e <- e + 1
    }

    #Now, for the functional analyses
    for (a in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[a]
        flog.info(paste("Making", analysis, "SummarizedExperiment"))

        #subset doses to contain only the analysis wanted
        analysisdoses <- NULL
        analysisdoses <- list()

        #Get names of Samples which have data for the analysis
        SamplesofInterest <- names(anallist)[which(sapply(1:length(anallist), function(x) { (analysis %in% anallist[[x]]) } ) == TRUE)]

        cts <- NULL

        if (is.null(restricttoLKTs)){
            #Dose of each analysis is total number of bases attributed to it
            for (ad in SamplesofInterest){
                tempanaldose <- subset(featuredoses[[ad]], Analysis == analysis)[, c("Accession", "NumBases", "Description")]
                tempanaldose$Accession <- as.character(tempanaldose$Accession)
                tempanaldose$Description <- as.character(tempanaldose$Description)
                tempanaldose$NumBases <- as.numeric(tempanaldose$NumBases)
                analysisdoses[[ad]] <- tempanaldose
            }
        } else {
            flog.info(paste("Making counts table restricted to the", length(restricttoLKTs),"LKTs specified."))
            #Dose of each analysis is sum of bases present in LKTs to restrict to
            for (ad in SamplesofInterest){
                #First, find out if taxa exist. If so, get numbers, else, write down 0.
                tempanaldose <- subset(featuredoses[[ad]], Analysis == analysis)
                taxapresent <- restricttoLKTs[restricttoLKTs %in% colnames(tempanaldose)]
                if (length(taxapresent) > 1){
                    numbavec <- as.numeric(rowSums(tempanaldose[ , taxapresent]))
                } else if (length(taxapresent) == 1) {
                    numbavec <- as.numeric(tempanaldose[ , taxapresent])
                } else {
                    numbavec <- rep(0,  nrow(tempanaldose))
                }
                tempanaldose <- tempanaldose[ , c("Accession", "Description")]
                tempanaldose$Accession <- as.character(tempanaldose$Accession)
                tempanaldose$Description <- as.character(tempanaldose$Description)
                tempanaldose$NumBases <- numbavec

                analysisdoses[[ad]] <- tempanaldose[, c("Accession", "NumBases", "Description")]
            }
        }

        phenoanal <- pheno2

        featureall <- bind_rows(analysisdoses, .id = "id")
        featureall[is.na(featureall)] <- 0
        colnames(featureall)[1] <- "Sample"

        cts <- spread(featureall, Sample, NumBases, fill = 0, drop = TRUE)
        #Just double check that there are no dupes
        if (length(which(duplicated(cts$Accession) == TRUE)) > 0){
            flog.info(paste("Found", length(which(duplicated(cts$Accession) == TRUE)), "duplicated accessions. Keeping the first one in each case."))
            cts <- cts[!(duplicated(cts$Accession)), ]
        }
        rownames(cts) <- cts$Accession
        cts$Accession <- NULL
        cts$Description <- NULL

        #Deal with samples with NO results for this analysis
        emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(cts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(cts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(cts)
            cts <- cbind(cts, complementarycts)
        }
        cts <- as.matrix(cts)

        ftt <- featureall[, c("Accession", "Description")]
        ftt <- ftt[!(duplicated(ftt$Accession)), ]
        rownames(ftt) <- ftt$Accession

        if (analysis == "resfinder"){
            data(resfinderlookup)
            ftt <- suppressMessages(left_join(ftt, resfinderlookup))
            ftt[] <- lapply(ftt, as.character)
            #Fill up reasonably with data which might be empty
            ftt[is.na(ftt[])] <- "S"
            for (colm in c("Description", "Class", "ResistancePhenotype", "PMID", "Mechanism")){
                ftt[, colm][which(ftt[, colm] == "S")] <- "none"
            }
            rownames(ftt) <- ftt$Accession
            ftt$Description <- ftt$Class
        }

        cts <- cts[order(rowSums(cts), decreasing = TRUE), ]
        featureorder <- rownames(cts)
        sampleorder <- rownames(phenoanal)
        ftt <- ftt[featureorder, ]
        cts <- cts[, sampleorder]

        assays <- list(cts)
        names(assays) <- "BaseCounts"

        ##Create SummarizedExperiment
        SEobj <- SummarizedExperiment(assays = assays, rowData = ftt, colData = phenoanal)
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$analysis <- analysis

        expvec[[e]] <- SEobj
        names(expvec)[e] <- analysis
        e <- e + 1

        #Add an extra one if converting resfinder to antibiogram
        if (analysis == "resfinder"){
            flog.info("Converting resfinder to antibiogram")
            expvec[[e]] <- make_antibiogram_experiment(SEobjresfinder = expvec[["resfinder"]])
            names(expvec)[e] <- "antibiogram"
            e <- e + 1
        }
    }

    return(expvec)
}
