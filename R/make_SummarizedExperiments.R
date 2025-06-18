#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL,  onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, restricttoLKTs = NULL, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, restricttoLKTs = NULL, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8){

    require(SummarizedExperiment)
    require(Matrix)
    data(MetaCycAccession2Description)
    data(ECdescmap)
    data(GOtermdict)
    data(InterproDict)
    data(blast_lookup)

    #Get data for features
    if (!is.null(onlysamples)){
        pheno2 <- pheno[rownames(pheno) %in% onlysamples, ]
    } else {
        pheno2 <- pheno
    }

    blastanalyses <- c("abricate", "vfdb", "resfinder", "plasmidfinder")

    Samples <- rownames(pheno2)
    possibleanalyses <- NULL

    if (any(c(is.null(onlyanalyses), !(all(c((length(onlyanalyses) == 1), ("LKT" %in% onlyanalyses))))))) {

        #Find out which functional analyses can be made
        anallist <- lapply(Samples, function (x) { as.character(unique(list.data[[paste(x, "abundances", sep = "_")]]$functional$Contig_LKT$Analysis)) } )
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

    if (all(c(!is.null(onlyanalyses), !(all(c((length(onlyanalyses) == 1), ("LKT" %in% onlyanalyses))))))){
        possibleanalyses <- possibleanalyses[possibleanalyses %in% onlyanalyses]
    } else {
        if (!is.null(possibleanalyses)) {
            if (!is.null(onlyanalyses)){
                possibleanalyses <- possibleanalyses[possibleanalyses %in% onlyanalyses]
            }
        } else {
            possibleanalyses <- onlyanalyses
        }

    }

    #Add colour table to expvec, if passed
    if (!is.null(cdict)){
        ctable <- plyr::rbind.fill(cdict)
        ctable <- ctable[!duplicated(ctable), ]
        rownames(ctable) <- ctable$Name
    } else {
        ctable <- NULL
    }

    #Make a vector for holding experiment list
    expvec <- list()
    e <- 1

    #Find out which taxonomic spaces are present
    taxspaceslist <- lapply(Samples, function (x) { as.character(names(list.data[[paste(x, "abundances", sep = "_")]]$taxonomic)) } )
    names(taxspaceslist) <- Samples

    alltaxspaces <- Reduce(union, taxspaceslist)
    numsampwithtaxspace <- NULL
    for (pa in 1:length(alltaxspaces)){
        possibletaxspace <- alltaxspaces[pa]
        numsampwithtaxspace[pa] <- length(which(sapply(1:length(taxspaceslist), function(x) { (possibletaxspace %in% taxspaceslist[[x]]) } ) == TRUE))
    }
    names(numsampwithtaxspace) <- alltaxspaces
    propsampleswithtaxspace <- numsampwithtaxspace / length(taxspaceslist)

    #Consider building an SEobj for a taxonomic space in which there are at least 80% of samples with that taxonomic space available.
    valid_taxonomic_spaces <- c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin")[c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin") %in% names(propsampleswithtaxspace)[propsampleswithtaxspace > 0.8]]

    for (taxonomic_space in valid_taxonomic_spaces){

        #Start by making LKT experiment
        flog.info(paste("Making", taxonomic_space, "SummarizedExperiment"))
        LKTdoses <- lapply(Samples, function (x) { list.data[[paste(x, "abundances", sep = "_")]]$taxonomic[[taxonomic_space]] } )
        names(LKTdoses) <- Samples

        LKTdosesall <- bind_rows(LKTdoses, .id = "id")
        colnames(LKTdosesall)[which(colnames(LKTdosesall) == "id")] <- "Sample"

        #De-duplicate features with the same name within the same sample #WILL REVIEW THIS ALGORITHM LATER
        LKTdosesall$SampFeat <- paste(LKTdosesall$Sample, LKTdosesall[ , taxonomic_space], sep = "§")
        dupes <- sort(LKTdosesall$SampFeat[duplicated(LKTdosesall$SampFeat)])
        for (dupe in dupes){
            LKTdosesall[which(LKTdosesall$SampFeat == dupe), taxonomic_space] <- make.unique(LKTdosesall[which(LKTdosesall$SampFeat == dupe), taxonomic_space])
        }

        #Make counts table
        LKTallcounts <- LKTdosesall[, c("Sample", taxonomic_space, "NumBases")]
        #Be sure that data is not empty or redundant
        LKTallcounts[is.na(LKTallcounts)] <- 0

        cts <- spread(LKTallcounts, Sample, NumBases, fill = 0, drop = FALSE)
        rownames(cts) <- cts[ , taxonomic_space]
        cts[ , taxonomic_space] <- NULL
        cts <- cts[order(rowSums(cts), decreasing = TRUE), ]
        featureorder <- rownames(cts)

        #Make tax table
        taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Taxid", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", taxonomic_space, "RefScore")]

        tt <- LKTdosesall[ , taxlvlspresent]
        tt <- tt[!(duplicated(tt)), ]
        #Fill in dark matter Taxids with "0"
        tt$Taxid[which(is.na(tt$Taxid))] <- "0"
        #Add Gram information
        data(JAMStaxtable)

        Taxid2gram <- JAMStaxtable[ , c("Taxid", "Gram")]
        Taxid2gram <- subset(Taxid2gram, Taxid %in% tt$Taxid)
        Taxid2gram <- Taxid2gram[!(duplicated(Taxid2gram$LKT)), ]
        tt <- left_join(as.data.frame(tt), Taxid2gram, by = "Taxid")
        tt[which(is.na(tt[ , "Gram"])), "Gram"] <- "na"
        tt <- tt[ , c("Gram", taxlvlspresent)]
        rownames(tt) <- tt[ , taxonomic_space]
        #tt <- as.matrix(tt)

        sampleorder <- rownames(pheno2)
        tt <- tt[featureorder, ]
        cts <- cts[, sampleorder]

        #Register the total number of NAHS bases sequenced for each sample
        TotBasesSamples <- colSums(cts)
        TotalBasesSequenced <- t(as.matrix(TotBasesSamples))
        rownames(TotalBasesSequenced) <- "NumBases"

        #Get genome completeness matrix
        LKTallGenComp <- LKTdosesall[, c("Sample", taxonomic_space, "Completeness")]
        #Be sure that data is not empty or redundant
        LKTallGenComp[is.na(LKTallGenComp)] <- 0
        LKTallGenCompcts <- spread(LKTallGenComp, Sample, Completeness, fill = 0, drop = FALSE)
        rownames(LKTallGenCompcts) <- LKTallGenCompcts[ , taxonomic_space]
        LKTallGenCompcts[ , taxonomic_space] <- NULL
        LKTallGenCompcts <- LKTallGenCompcts[featureorder, sampleorder]

        #Get genome Quality matrix
        LKTallGenQual <- LKTdosesall[, c("Sample", taxonomic_space, "Quality")]
        #Be sure that data is not empty or redundant
        LKTallGenQual[is.na(LKTallGenQual)] <- "N_A"
        LKTallGenQualcts <- spread(LKTallGenQual, Sample, Quality, fill = "N_A", drop = FALSE)
        rownames(LKTallGenQualcts) <- LKTallGenQualcts[ , taxonomic_space]
        LKTallGenQualcts[ , taxonomic_space] <- NULL
        LKTallGenQualcts <- LKTallGenQualcts[featureorder, sampleorder]

        assays <- list(as.matrix(cts), as.matrix(LKTallGenCompcts), as.matrix(LKTallGenQualcts))
        names(assays) <- c("BaseCounts", "GenomeCompleteness", "GenomeQuality")

        SEobj <- SummarizedExperiment(assays = assays, rowData = tt, colData = pheno2)
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequenced #LKT is the special case in which all bases sequenced are for the analysis
        metadata(SEobj)$analysis <- taxonomic_space
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }

        if (!is.null(ctable)){
            metadata(SEobj)$ctable <- ctable
        }

        if (is.null(restricttoLKTs)){
            expvec[[e]] <- SEobj
            names(expvec)[e] <- taxonomic_space
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

    #Now, for the functional analyses.
    #Obtain feature doses for every taxonomic space available before stratifying

    if (stratify_functions_by_taxon){
        flog.info(paste("Stratifying all features by taxa, this might take a while."))
        master_sparse_featcounts_matrix_list <- list()
        master_sparse_featcounts_index_list <- list()
        master_sparse_taxon_to_space_df <- NULL

        for (taxonomic_space in valid_taxonomic_spaces){
            featuredoses <- lapply(Samples, function (x) { list.data[[paste(x, "abundances", sep = "_")]]$functional[[taxonomic_space]] } )
            names(featuredoses) <- Samples
            featuredosesall <- dplyr::bind_rows(featuredoses, .id = "Sample")
            featuredosesall[is.na(featuredosesall)] <- 0
            Taxoncols <- colnames(featuredosesall)[6:ncol(featuredosesall)]
            featuredosesall$SampAcc <- paste(featuredosesall$Sample, featuredosesall$Accession, sep = "§")
            featuredosesall$Taxonomic_space <- taxonomic_space
            rownames(featuredosesall) <- featuredosesall$SampAcc
            #Feature index order will come from Contig_LKT space to ensure consistency, as will the overall counts
            if (taxonomic_space == "Contig_LKT"){
                master_index_row_order <- rownames(featuredosesall)
                featurecountsall <- featuredosesall[ , c("Sample", "Analysis", "Accession", "Description", "NumBases")]
            } else {
                featuredosesall <- featuredosesall[master_index_row_order, ]
            }
            #Rearrange for clarity
            featuredosesall <- featuredosesall[ , c("Sample", "Analysis", "Accession", "Description", "NumBases", "SampAcc", "Taxonomic_space", Taxoncols)]
            #Split counts from feature information
            currFeatdose <- featuredosesall[ , Taxoncols]
            #Rename Taxon columns with taxonomic_space because taxon names may be repeated across spaces
            colnames(currFeatdose) <- paste(colnames(currFeatdose), taxonomic_space, sep = "§")
            master_sparse_featcounts_matrix_list[[taxonomic_space]] <- Matrix::Matrix(data = as.matrix(currFeatdose), sparse = TRUE)
            curr_master_index <- featuredosesall[ , c("Sample", "Analysis", "Accession", "Description", "NumBases", "Taxonomic_space")]
            master_sparse_featcounts_index_list[[taxonomic_space]] <- curr_master_index
            #Add taxonomic_space to feature information 

            master_sparse_taxon_to_space_df <- rbind(master_sparse_taxon_to_space_df, data.frame(Taxon = Taxoncols, Taxonomic_space = taxonomic_space, Sparse_matrix_column_name = colnames(currFeatdose)))
            rownames(master_sparse_taxon_to_space_df) <- master_sparse_taxon_to_space_df$Sparse_matrix_column_name

            #Clean up for keeping RAM down, although I'm not sure this'll do anything at all.
            rm(featuredoses)
            rm(featuredosesall)
            gc()
        }

        #bind columns of all taxonomic spaces into a single matrix
        master_sparse_featcounts_matrix <- do.call(cbind, master_sparse_featcounts_matrix_list)
        master_sparse_featcounts_index <- master_sparse_featcounts_index_list[["Contig_LKT"]][ , c("Sample", "Analysis", "Accession", "Description")]
        master_sparse_featcounts_index$RowNumber <- 1:nrow(master_sparse_featcounts_index)

        #Clean up for keeping RAM down, although I'm not sure this'll do anything at all.
        rm(master_sparse_featcounts_matrix_list)
        gc()

    } else {

        featuredoses <- lapply(Samples, function (x) { list.data[[paste(x, "abundances", sep = "_")]]$functional[["Contig_LKT"]] } )
        names(featuredoses) <- Samples
        featuredosesall <- bind_rows(featuredoses, .id = "Sample")
        featurecountsall <- featuredoses[ , c("Sample", "Analysis", "Accession", "Description", "NumBases")]
        rm(featuredoses)
        gc()

    }

    #Gather all gene information
    featuredata <- list.data[paste(Samples, "featuredata", sep = "_")]
    names(featuredata) <- Samples
    featuredataall <- bind_rows(featuredata, .id = "Sample")
    rm(featuredata)
    gc()

    for (analnumb in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[analnumb]
        flog.info(paste("Making", analysis, "SummarizedExperiment"))

        #subset doses to contain only the analysis wanted
        analysisdoses <- NULL
        analysisdata <- NULL

        #Make BaseCounts table
        cts <- NULL
        phenoanal <- pheno2
        analysisdoses <- featurecountsall[which(featurecountsall$Analysis == analysis), ]

        #Fix for not getting duplicate accessions. Delete the Description column to get clean pivot from long to wide.
        cts <- analysisdoses[ , c("Sample", "Accession", "NumBases")] %>% pivot_wider(names_from = Sample, values_from = NumBases, values_fill = 0)
        cts <- as.data.frame(cts)
        cts[is.na(cts)] <- 0

        #Just double check that there are no dupes
        if (length(which(duplicated(cts$Accession) == TRUE)) > 0){
            flog.info(paste("Found", length(which(duplicated(cts$Accession) == TRUE)), "duplicated accessions. Keeping the first one in each case."))
            cts <- cts[!(duplicated(cts$Accession)), ]
        }
        rownames(cts) <- cts$Accession
        cts$Accession <- NULL

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

        #Make GeneCounts counts table
        gennumcts <- NULL
        genlengthcts <- NULL
        analysisdata <- featuredataall[ , c("Sample", "LengthDNA", analysis)]
        analysisdata$LengthDNA <- as.numeric(analysisdata$LengthDNA)
        colnames(analysisdata)[which(colnames(analysisdata) == analysis)] <- "Accession"
        #Fix "none" accession names to match analysisdoses
        analysisdata[which(analysisdata$Accession == "none"), "Accession"] <- paste(analysis, analysisdata[which(analysisdata$Accession == "none"), "Accession"], sep = "_")

        if (analysis %in% c("GO", "Interpro", "MetaCyc")){

            #Separate Accesions delimited by a "|", as a single gene may have multiple accessions within this analysis space 
            gennumlong <- analysisdata[ , c("Sample", "Accession")] %>% tidyr::separate_rows(all_of("Accession"), sep = fixed("\\|")) %>% group_by(across(all_of(c("Sample", "Accession")))) %>% summarize(NumGenes = length(Accession), .groups = "keep")
            genlengthlong <- analysisdata %>% tidyr::separate_rows(all_of("Accession"), sep = fixed("\\|")) %>% group_by(across(all_of(c("Sample", "Accession")))) %>% summarize(GeneLengthSum = sum(LengthDNA), .groups = "keep")

        } else {

            gennumlong <- analysisdata[ , c("Sample", "Accession")] %>% group_by(across(all_of(c("Sample", "Accession")))) %>% summarize(NumGenes = length(Accession), .groups = "keep")
            genlengthlong <- analysisdata %>% group_by(across(all_of(c("Sample", "Accession")))) %>% summarize(GeneLengthSum = sum(LengthDNA), .groups = "keep")

        }

        gennumcts <- gennumlong %>% pivot_wider(names_from = Sample, values_from = NumGenes, values_fill = 0)
        gennumcts <- as.data.frame(gennumcts)
        rownames(gennumcts) <- gennumcts$Accession
        gennumcts$Accession <- NULL
        #Deal with samples with NO results for this analysis
        emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(gennumcts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(gennumcts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(gennumcts)
            genlengthcts <- cbind(gennumcts, complementarycts)
        }
        gennumcts <- as.matrix(gennumcts)

        #Make GeneLengths counts table
        genlengthcts <- genlengthlong %>% pivot_wider(names_from = Sample, values_from = GeneLengthSum, values_fill = 0)
        genlengthcts <- as.data.frame(genlengthcts)
        rownames(genlengthcts) <- genlengthcts$Accession
        genlengthcts$Accession <- NULL
        #Deal with samples with NO results for this analysis
        emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(genlengthcts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(genlengthcts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(genlengthcts)
            genlengthcts <- cbind(genlengthcts, complementarycts)
        }
        genlengthcts <- as.matrix(genlengthcts)

        ## Make feature table
        # Account for the fact that descriptions for the same accession may diverge, due to using jams files annotated by different versions. The function below will choose the most prevalent weighted by sequencing depth.
        find_best_description <- function(featdf = NULL){
            descstats <- as.data.frame(featdf) %>% group_by(Description) %>% summarize(Sum = sum(NumBases))
            best_description <- descstats$Description[which(descstats$Sum == max(descstats$Sum))]

            return(best_description)
        }

        ftt <- data.frame(Accession = unique(analysisdoses$Accession))
        if (!(analysis %in% c("ECNumber", "MetaCyc", "GO"))){
            ftt$Description <- sapply(ftt$Accession, function (x) { find_best_description(featdf = subset(analysisdoses, Accession == x)) } )
        }

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
            ftt <- suppressMessages(left_join(ftt, blast_lookup$resfinder_lookup))
            ftt[] <- lapply(ftt, as.character)
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

        #Order gennumcts and genlengthcts
        gennumcts <- gennumcts[featureorder, sampleorder]
        genlengthcts <- genlengthcts[featureorder, sampleorder]

        assays <- list()
        assays$BaseCounts <- cts
        assays$GeneCounts <- gennumcts
        assays$GeneLengthSum <- genlengthcts

        ##Create SummarizedExperiment
        SEobj <- SummarizedExperiment(assays = assays, rowData = as.matrix(ftt), colData = as.matrix(phenoanal))
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequencedinAnalysis
        metadata(SEobj)$analysis <- analysis
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }
        if (!is.null(ctable)){
            metadata(SEobj)$ctable <- ctable
        }

        #Split functions by taxon, if applicable
        if (stratify_functions_by_taxon){
            flog.info(paste("Splitting", analysis, "features by taxa."))
            flog.info(paste("Splitting", analysis, "feature BaseCounts into their contributing taxa."))

            allfeaturesbytaxa_index <- master_sparse_featcounts_index[which(master_sparse_featcounts_index$Analysis == analysis), ]
            allfeaturesbytaxa_matrix <- master_sparse_featcounts_matrix[rownames(allfeaturesbytaxa_index), ]

            allfeaturesbytaxa_index$RowNumber <- 1:nrow(allfeaturesbytaxa_matrix)

            #Prune empty LKTs
            NonEmptyLKTs <- names(which(Matrix::colSums(allfeaturesbytaxa_matrix) != 0))
            allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix[ , NonEmptyLKTs]

            #Add features-by-taxon matrix and index into SummarizedExperiment object
            metadata(SEobj)$allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix
            metadata(SEobj)$allfeaturesbytaxa_index <- allfeaturesbytaxa_index
            metadata(SEobj)$sparse_taxon_to_space_df <- master_sparse_taxon_to_space_df

            rm(allfeaturesbytaxa_matrix)
            rm(allfeaturesbytaxa_index)


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
