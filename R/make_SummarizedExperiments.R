#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL,  onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8){

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

        #To avoid duplicated names of MB2 bins across samples, paste Sample name to these bins. This will ensure there are no two different entities with the same SGB name.
        #Sorry about the for loop, but it is safer like so.
        for (SN in names(LKTdoses)){
            LKTdoses[[SN]][grep("^MB2__", LKTdoses[[SN]][ , taxonomic_space]), taxonomic_space] <- sapply(grep("^MB2__", LKTdoses[[SN]][ , taxonomic_space]), function (x) { gsub("^MB2__", paste0("MB2__", SN , "__"), LKTdoses[[SN]][x , taxonomic_space]) })
            #Make sure there are no leftover duplicates.
            LKTdoses[[SN]][ , taxonomic_space] <- make.unique(LKTdoses[[SN]][ , taxonomic_space])
        }

        LKTdosesall <- bind_rows(LKTdoses, .id = "Sample")

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
        taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Taxid", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT", taxonomic_space, "NCBI_taxonomic_rank", "RefScore")]

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
        metadata(SEobj)$version <- opt$verstr
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }

        if (!is.null(ctable)){
            metadata(SEobj)$ctable <- ctable
        }

        expvec[[e]] <- SEobj
        names(expvec)[e] <- taxonomic_space
        e <- e + 1
        #Clean memory up
        gc()
    }

    if (!is.null(onlyanalyses)){
        if (all(c((length(onlyanalyses) == 1), ("LKT" %in% onlyanalyses)))) {
            #stop here and return LKT SummarizedExperiment
            return(expvec)
        }
    }

    #Now, for the functional analyses.
    for (analnumb in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[analnumb]
        flog.info(paste("Making", analysis, "SummarizedExperiment"))

        #subset doses to contain only the analysis wanted
        analysisdoses <- NULL
        analysisdata <- NULL

        #Make BaseCounts table
        cts <- NULL
        phenoanal <- pheno2

        #This data frame will hold the unstratified overall counts for all features and all samples.
        analysisdoses <- NULL
        analysisgencts <- NULL

        if (!stratify_functions_by_taxon){
            #If not stratifying functions by taxa, get functional information Contig_LKT abundances
            valid_taxonomic_space <- "Contig_LKT"
        } else {
            #Declare empty lists to hold data
            master_sparse_featcounts_index <- NULL
            master_sparse_taxon_to_space_basecounts_df <- NULL
            master_sparse_taxon_to_space_genecounts_df <- NULL

            #When stratifying, prioritize Consolidated Genome bin if present, if not, default to Contig_LKT
            valid_taxonomic_space <- c("ConsolidatedGenomeBin", "Contig_LKT")[c("ConsolidatedGenomeBin", "Contig_LKT") %in% valid_taxonomic_spaces][1]
        }

        for (SN in Samples){
            curr_FD <- list.data[[paste(SN, "abundances", sep = "_")]]$functional[[valid_taxonomic_space]]
            #Subset to current analysis
            curr_FD <- curr_FD[which(curr_FD$Analysis == analysis), , drop = FALSE]
            #Check that there is any data for that particular sample - there may not be.
            if (nrow(curr_FD) > 0){
                #There shouldn't be any NAs but if there are make them 0
                curr_FD[is.na(curr_FD)] <- 0
                Taxoncols <- colnames(curr_FD)[5:ncol(curr_FD)]
                if (valid_taxonomic_space == "ConsolidatedGenomeBin"){
                    #To avoid duplicated names of MB2 bins across samples, paste Sample name to these bins. This will ensure there are no two different entities with the same SGB name.
                    Taxoncols[grep("^MB2__", Taxoncols)] <- sapply(Taxoncols[grep("^MB2__", Taxoncols)], function (x) { gsub("^MB2__", paste0("MB2__", SN , "__"), x) })
                    colnames(curr_FD)[5:ncol(curr_FD)] <- Taxoncols
                }
                curr_FD$Sample <- SN
                curr_FD$SampAcc <- paste(curr_FD$Sample, curr_FD$Accession, sep = "§")
                rownames(curr_FD) <- curr_FD$SampAcc

                #Deal with gene number tallies
                curr_featdf <- list.data[[paste(SN, "featuredata", sep = "_")]][ , c("Feature", "LengthDNA", analysis, valid_taxonomic_space)]
                #Fix taxon names if ConsolidatedGenomeBin
                if (valid_taxonomic_space == "ConsolidatedGenomeBin"){
                    curr_featdf$ConsolidatedGenomeBin <- sapply(curr_featdf$ConsolidatedGenomeBin, function (x) { gsub("^MB2__", paste0("MB2__", SN , "__"), x) })
                }
                colnames(curr_featdf)[which(colnames(curr_featdf) == analysis)] <- "Accession"
                colnames(curr_featdf)[which(colnames(curr_featdf) == valid_taxonomic_space)] <- "Taxon"
                #Fix "none" to analysis_none
                curr_featdf$Accession[which(curr_featdf$Accession == "none")] <- paste(analysis, "none", sep = "_")

                #Separate rows if GO, Interpro or MataCyc
                if (analysis %in% c("GO", "Interpro", "MetaCyc")){
                    curr_featdf <- curr_featdf %>% tidyr::separate_rows(all_of("Accession"), sep = fixed("\\|"))
                }

                if (stratify_functions_by_taxon){

                    curr_SM <- Matrix::Matrix(data = as.matrix(curr_FD[ , Taxoncols]), sparse = TRUE)
                    colnames(curr_SM) <- Taxoncols

                    #In order to conserve RAM and the size of the taxonomic stratification sparse matrix, do not include any taxon whose MAXIMUM number of bases for a SINGLE feature is less than 1 PPM of the sequencing depth IN THAT SAMPLE, as this is likely a spurious taxonomic hit.
                    curr_Bases_equal_to_1PPM <- round((metadata(expvec[[valid_taxonomic_space]])$TotalBasesSequenced["NumBases", SN] / 1E6), 0)
                    curr_taxaToKeep <- colnames(curr_SM)[colSums(curr_SM[rownames(curr_SM)[!(rownames(curr_SM) %in% paste(SN, paste(analysis, "none", sep = "_"), sep = "§"))], , drop = FALSE]) >=  curr_Bases_equal_to_1PPM]

                    #Proceed only if there is anything worthwile to merge
                    if (length(curr_taxaToKeep) > 0) {
                        #Ok, we're going to bank stuff so commit to index
                        master_sparse_featcounts_index <- rbind(master_sparse_featcounts_index, curr_FD[ , c("Sample", "Accession", "Description", "NumBases")])

                        curr_SM <- curr_SM[ , curr_taxaToKeep, drop = FALSE]
                        master_sparse_taxon_to_space_basecounts_df <- merge_sparse_matrix(matlist = list(master_sparse_taxon_to_space_basecounts_df, curr_SM))
                        #Ensure correct row order
                        master_sparse_taxon_to_space_basecounts_df <- master_sparse_taxon_to_space_basecounts_df[rownames(master_sparse_featcounts_index), , drop = FALSE]

                        #Deal with gene number stratifications
                        #Pivot to wide form
                        curr_split_numgenes <- as.data.frame(curr_featdf %>% group_by(across(all_of(c("Taxon", "Accession")))) %>% summarise(NumGenes = length(Accession), .groups = "keep") %>% pivot_wider(names_from = Taxon, values_from = NumGenes, values_fill = 0))
                        rownames(curr_split_numgenes) <- paste(SN, curr_split_numgenes$Accession, sep = "§")
                        curr_split_numgenes$Accession <- NULL

                        #Curtail to the same taxa as for basecounts
                        #All taxa from curr_SM should be in curr_split_numgenes, but will account for any missing ones. Features are non-negotiable as to maintain the same row index.
                        curr_split_numgenes <- Matrix::Matrix(data = as.matrix(curr_split_numgenes[rownames(curr_SM), colnames(curr_SM)[colnames(curr_SM) %in% colnames(curr_split_numgenes)], drop = FALSE]), sparse = TRUE)

                        #Merge sparse version in to master gene count matrix
                        master_sparse_taxon_to_space_genecounts_df <- merge_sparse_matrix(matlist = list(master_sparse_taxon_to_space_genecounts_df, curr_split_numgenes))

                        #Ensure correct row order
                        master_sparse_taxon_to_space_genecounts_df <- master_sparse_taxon_to_space_genecounts_df[rownames(master_sparse_featcounts_index), ]

                    }# end conditional that there are any taxa > 1 PPM

                    #Clean up
                    curr_SM <- NULL
                    curr_split_numgenes <- NULL
                    curr_Bases_equal_to_1PPM <- NULL
                    curr_taxaToKeep <- NULL
                    gc()
                }

                #Eliminate Taxoncols to save RAM
                analysisdoses <- rbind(analysisdoses, curr_FD[ , c("Sample", "Accession", "Description", "NumBases")])

                curr_featdf$Sample <- SN
                curr_featdf$Taxon <- NULL
                curr_analysisgencts <- as.data.frame(curr_featdf %>% group_by(Accession) %>% summarise(NumGenes = length(Accession), .groups = "keep"))
                curr_analysisgencts$Sample <- SN
                rownames(curr_analysisgencts) <- paste(SN, curr_analysisgencts$Accession, sep = "§")
                curr_analysisgencts <- curr_analysisgencts[rownames(curr_FD), c("Sample", "Accession", "NumGenes")]
                analysisgencts <- rbind(analysisgencts, curr_analysisgencts)

                #Clean up
                curr_FD <- NULL
                curr_featdf <- NULL
                curr_analysisgencts <- NULL

                gc()
            } #End conditional for there being data for that analysis for that sample
        } #End loop for obtaining data for each sample

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


       #Fix for not getting duplicate accessions. Delete the Description column to get clean pivot from long to wide.
        gennumcts <- analysisgencts[ , c("Sample", "Accession", "NumGenes")] %>% pivot_wider(names_from = Sample, values_from = NumGenes, values_fill = 0)
        gennumcts <- as.data.frame(gennumcts)
        gennumcts[is.na(gennumcts)] <- 0
        rownames(gennumcts) <- gennumcts$Accession
        gennumcts$Accession <- NULL

        #Deal with samples with NO results for this analysis
        emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(gennumcts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(gennumcts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(gennumcts)
            gennumcts <- cbind(gennumcts, complementarycts)
        }
        gennumcts <- as.matrix(gennumcts)

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
            ftt <- left_join(ftt, blast_lookup$resfinder_lookup, by = "Accession")
            ftt <- replace_NAs_with_character(df = ftt, replacement = "N_A")
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

        assays <- list()
        assays$BaseCounts <- cts
        assays$GeneCounts <- gennumcts
 
        ##Create SummarizedExperiment
        SEobj <- SummarizedExperiment(assays = assays, rowData = as.matrix(ftt), colData = as.matrix(phenoanal))
        metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
        metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequencedinAnalysis
        metadata(SEobj)$analysis <- analysis
        metadata(SEobj)$version <- opt$verstr
        if (!is.null(phenolabels)){
            metadata(SEobj)$phenolabels <- phenolabels
        }
        if (!is.null(ctable)){
            metadata(SEobj)$ctable <- ctable
        }

        #Split functions by taxon, if applicable
        if (stratify_functions_by_taxon){
 
            master_sparse_featcounts_index$RowNumber <- 1:nrow(master_sparse_featcounts_index)
            metadata(SEobj)$allfeaturesbytaxa_index <- master_sparse_featcounts_index

            #Add features-by-taxon matrix and index into SummarizedExperiment object
            #Prune empty LKTs, there shoud be none, but be safe
            NonEmptyLKTs <- names(which(Matrix::colSums(master_sparse_taxon_to_space_basecounts_df) != 0))
            metadata(SEobj)$allfeaturesbytaxa_matrix <- master_sparse_taxon_to_space_basecounts_df[rownames(master_sparse_featcounts_index), NonEmptyLKTs]

            metadata(SEobj)$allfeaturesbytaxa_GeneCounts_matrix <- master_sparse_taxon_to_space_genecounts_df[rownames(master_sparse_featcounts_index), NonEmptyLKTs]

            #Clean up
            master_sparse_featcounts_index <- NULL
            master_sparse_taxon_to_space_basecounts_df <- NULL
            master_sparse_taxon_to_space_genecounts_df <- NULL
            gc()

        } #End of splitting features by taxonomy

        expvec[[e]] <- SEobj
        names(expvec)[e] <- analysis
        e <- e + 1
        #clean memory up
        gc()
    }

    return(expvec)
}
