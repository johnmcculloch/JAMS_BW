#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL,  onlyanalyses = NULL, minnumsampanalysis = NULL, minpropsampanalysis = 0.1, minPctFromCtg = NULL, minProbNumGenomes = NULL, restricttoLKTs = NULL, list.data = NULL)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL,  onlyanalyses = NULL, minnumsampanalysis = NULL, minpropsampanalysis = 0.1, restricttoLKTs = NULL, stratify_functions_by_taxon = TRUE, list.data = NULL, phenolabels = NULL){

    require(SummarizedExperiment)
    require(Matrix)
    data(MetaCycAccession2Description)
    data(ECdescmap)
    data(GOtermdict)
    data(InterproDict)

    #Get data for features
    if (!is.null(onlysamples)){
        pheno2 <- pheno[rownames(pheno) %in% onlysamples, ]
    } else {
        pheno2 <- pheno
    }

    blastanalyses <- c("abricate", "vfdb", "resfinder", "plasmidfinder")

    Samples <- rownames(pheno2)
    if (any(c(is.null(onlyanalyses), !(all(c((length(onlyanalyses) == 1), ("LKT" %in% onlyanalyses))))))) {

        featuredoses <- list.data[paste(Samples, "featuredose", sep="_")]
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
    }

    if (!is.null(onlyanalyses)){
        possibleanalyses <- possibleanalyses[possibleanalyses %in% onlyanalyses]
    }

    #Make a vector for holding experiment list
    expvec <- list()
    e <- 1

    if (TRUE){
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
        #Add Gram information
        data(Gram)
        phylum2gram <- Gram[,c("Phylum", "Gram")]
        tt <- left_join(as.data.frame(tt), phylum2gram, by = "Phylum")
        tt[which(is.na(tt[ , "Gram"])), "Gram"] <- "na"
        tt <- tt[ , c("Gram", taxlvlspresent)]
        rownames(tt) <- tt$LKT
        tt <- as.matrix(tt)

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

        assays <- list(as.matrix(cts), as.matrix(LKTPctFromCtgcts), as.matrix(LKTallGenCompcts))
        names(assays) <- c("BaseCounts", "PctFromCtgs", "GenomeCompleteness")

        SEobj <- SummarizedExperiment(assays = assays, rowData = as.matrix(tt), colData = as.matrix(pheno2))
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequenced #LKT is the special case in which all bases sequenced are for the analysis
        metadata(SEobj)$analysis <- "LKT"
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }
        if (is.null(restricttoLKTs)){
            expvec[[e]] <- SEobj
            names(expvec)[e] <- "LKT"
            e <- e + 1
            #Clean memory up
            gc()
        }
    }

    if (!is.null(onlyanalyses)){
        if (all(c((length(onlyanalyses) == 1), ("LKT" %in% onlyanalyses)))) {
            #stop here and return LKT SummarizedExperiment
            return(expvec)
        }
    }

    #Now, for the functional analyses
    featuredoses <- list.data[paste(Samples, "featuredose", sep="_")]
    names(featuredoses) <- Samples
    featuredata <- list.data[paste(Samples, "featuredata", sep="_")]
    names(featuredata) <- Samples

    for (a in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[a]
        flog.info(paste("Making", analysis, "SummarizedExperiment"))

        #subset doses to contain only the analysis wanted
        analysisdoses <- NULL
        analysisdoses <- list()
        analysisdata <- NULL
        analysisdata <- list()

        #Get names of Samples which have data for the analysis
        SamplesofInterest <- names(anallist)[which(sapply(1:length(anallist), function(x) { (analysis %in% anallist[[x]]) } ) == TRUE)]

        cts <- NULL

            if (is.null(restricttoLKTs)){
                #Dose of each analysis is total number of bases attributed to it
                for (ad in SamplesofInterest){
                    tempanaldose <- subset(featuredoses[[ad]], Analysis == analysis)[, c("Accession", "NumBases", "Description")]
                    tempanaldose$Accession <- as.character(tempanaldose$Accession)
                    #For ECNumber, GO and MetaCyc, add descriptions safely, based on latest info, later in ftt
                    if (analysis %in% c("ECNumber", "GO", "MetaCyc")){
                        tempanaldose$Description <- "none"
                    } else {
                        tempanaldose$Description <- as.character(tempanaldose$Description)
                    }
                    tempanaldose$NumBases <- as.numeric(tempanaldose$NumBases)
                    analysisdoses[[ad]] <- tempanaldose
                    tempanaldata <- featuredata[[ad]]
                    #Careful to make it backwards-compatible
                    featdatacolsavail <- colnames(tempanaldata)[colnames(tempanaldata) %in% c("Feature", "LKT", "GeneName", analysis)]
                    tempanaldata <- tempanaldata[ , featdatacolsavail]
                    #If analysis contains several accessions per gene, repeat rows for each kind
                    #Also, make it backwards compatible with older jamsfiles
                    if (all(c((analysis %in% colnames(tempanaldata)), (analysis %in% c("GO", "MetaCyc"))))){
                        #flog.info("Analysis contains several accessions per gene, splitting to get number of genes per accession.")
                        tempanaldata <- tidyr::separate_rows(tempanaldata, all_of(analysis), sep = fixed("\\|"))
                    }
                    #If possible, get feature count
                    if (analysis %in% colnames(tempanaldata)){
                        featcount <- as.data.frame(table(tempanaldata[ , analysis]), stringsAsFactors = FALSE)
                        colnames(featcount) <- c("Accession", "Count")
                        rownames(featcount) <- featcount$Accession
                        analysisdata[[ad]] <- featcount
                    } else {
                        analysisdata[[ad]] <- NULL
                    }
                    #Add GeneName to description if Product
                    if (all(c((analysis %in% c("Product", "vfdb")), ("GeneName" %in% featdatacolsavail)) )){
                        tempanaldose$Description <- tempanaldata$GeneName[match(tempanaldose$Accession, tempanaldata$Product)]
                        analysisdoses[[ad]] <- tempanaldose
                    }
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
                    if (analysis %in% c("ECNumber", "GO", "MetaCyc")){
                        tempanaldose$Description <- "none"
                    } else {
                        tempanaldose$Description <- as.character(tempanaldose$Description)
                    }
                    tempanaldose$NumBases <- numbavec
                    analysisdoses[[ad]] <- tempanaldose[, c("Accession", "NumBases", "Description")]
                }
            } #End conditional of restricting to LKTs

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
        #Fix descriptions if ECNumber, GO or MetaCyc
        if (analysis == "ECNumber"){
            ftt$Description <- NULL
            ftt <- left_join(ftt, ECdescmap, by = "Accession")
        } else if (analysis == "MetaCyc"){
            ftt$Description <- NULL
            ftt <- left_join(ftt, MetaCycAccession2Description, by = "Accession")
        } else if (analysis == "GO"){
            ftt$Description <- NULL
            ftt <- left_join(ftt, GOtermdict, by = "Accession")
        }
        rownames(ftt) <- ftt$Accession

        #Make sure there is no missing information
        for (colm in 1:ncol(ftt)){
            ftt[which(is.na(ftt[ , colm])), colm] <- "none"
        }

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
        cts <- as.matrix(cts)

        #Register the total number of bases sequenced for each sample within that analysis
        TotalBasesSequencedinAnalysis <- colSums(cts)
        TotalBasesSequencedinAnalysis <- t(as.matrix(TotalBasesSequencedinAnalysis))
        rownames(TotalBasesSequencedinAnalysis) <- "NumBases"

        #If there is how to count numbers of genes for each accession, then do that.
        if (length(analysisdata) == length(analysisdoses)){
            featuredataall <- bind_rows(analysisdata, .id = "id")
            featuredataall[is.na(featuredataall)] <- 0
            colnames(featuredataall)[1] <- "Sample"

            featcts <- spread(featuredataall, Sample, Count, fill = 0, drop = TRUE)
            #Just double check that there are no dupes
            if (length(which(duplicated(featcts$Accession) == TRUE)) > 0){
                flog.info(paste("Found", length(which(duplicated(featcts$Accession) == TRUE)), "duplicated accessions. Keeping the first one in each case."))
                featcts <- featcts[!(duplicated(featcts$Accession)), ]
            }
            featcts$Accession[which(featcts$Accession == "none")] <- paste(analysis, "none", sep = "_")
            rownames(featcts) <- featcts$Accession
            featcts$Accession <- NULL

            if (length(emptySamples) > 0){
                complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(featcts), data = 0)
                complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
                colnames(complementarycts) <- emptySamples
                rownames(complementarycts) <- rownames(featcts)
                featcts <- cbind(featcts, complementarycts)
            }
            featcts <- featcts[featureorder, sampleorder]
            featcts <- as.matrix(featcts)
        } else {
            featcts <- NULL
        }

        assays <- list()
        assays$BaseCounts <- cts
        if (!is.null(featcts)){
            assays$GeneCounts <- featcts
        }

        ##Create SummarizedExperiment
        SEobj <- SummarizedExperiment(assays = assays, rowData = as.matrix(ftt), colData = as.matrix(phenoanal))
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequencedinAnalysis
        metadata(SEobj)$analysis <- analysis
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }

        #Split functions by taxon, if applicable
        if (stratify_functions_by_taxon){
            flog.info(paste("Splitting", analysis, "features by taxa. This may take a while."))

            featurebytaxonlist <- list()
            SamplesWanted <- colnames(cts)[!(colnames(cts) %in% emptySamples)]
            wantedfeatures <- rownames(cts)
            NumBases1PPMthreshold <- TotalBasesSequencedinAnalysis / 1000000
            for(samp in SamplesWanted){
                currFeatdose <- NULL
                #Get appropriate object with counts
                wantedobj <- paste(samp, "featuredose", sep = "_")
                currFeatdose <- list.data[[wantedobj]]
                currFeatdose <- subset(currFeatdose, Analysis == analysis)
                currFeatdose$Analysis <- NULL
                rownames(currFeatdose) <- currFeatdose$Accession
                currwantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(currFeatdose)]
                currFeatdose <- currFeatdose[currwantedfeatures, ]
                currFeatdose$Accession <- NULL
                currFeatdose$Description <- NULL
                currFeatdose$NumBases <- NULL
                signalthreshold <- NumBases1PPMthreshold["NumBases" , samp]
                NonEmptyTaxa <- names(which(colSums(currFeatdose) > 0))
                currnonemptyFeatdose <- as.data.frame(currFeatdose[, NonEmptyTaxa])
                rownames(currnonemptyFeatdose) <- currwantedfeatures
                colnames(currnonemptyFeatdose) <- NonEmptyTaxa
                currnonemptyFeatdose$Sample <- samp
                currnonemptyFeatdose$Accession <- rownames(currnonemptyFeatdose)
                featurebytaxonlist[[samp]] <- currnonemptyFeatdose
            }

            allfeaturesbytaxa <- plyr::rbind.fill(featurebytaxonlist)
            allfeaturesbytaxa[is.na(allfeaturesbytaxa)] <- 0

            allfeaturesbytaxa <- allfeaturesbytaxa[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa)[which(!colnames(allfeaturesbytaxa) %in% c("Sample", "Accession"))])))]

            SampleAccession2Row <- allfeaturesbytaxa[c("Sample", "Accession")]
            SampleAccession2Row$RowNumber <- 1:nrow(SampleAccession2Row)

            allfeaturesbytaxa$Sample <- NULL
            allfeaturesbytaxa$Accession <- NULL
            allfeaturesbytaxa <- as.matrix(allfeaturesbytaxa)
            rownames(allfeaturesbytaxa) <- 1:nrow(allfeaturesbytaxa)

            allfeaturesbytaxa_matrix <- Matrix::Matrix(data = allfeaturesbytaxa, sparse = TRUE)
            allfeaturesbytaxa <- NULL

            #Add features-by-taxon matrix and index into SummarizedExperiment object
            metadata(SEobj)$allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix
            metadata(SEobj)$allfeaturesbytaxa_index <- SampleAccession2Row

        } #End of splitting features by taxonomy

        expvec[[e]] <- SEobj
        names(expvec)[e] <- analysis
        e <- e + 1
        #clean memory up
        gc()

        #Add an extra one if converting resfinder to antibiogram
        #if (analysis == "resfinder"){
        #    flog.info("Converting resfinder to antibiogram")
        #    expvec[[e]] <- make_antibiogram_experiment(SEobjresfinder = expvec[["resfinder"]])
        #    names(expvec)[e] <- "antibiogram"
        #    e <- e + 1
        #}
    }

    return(expvec)
}
