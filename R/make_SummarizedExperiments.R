#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8, functional_contextualization_dissimilarity_cutoff = 0.15)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8, functional_contextualization_dissimilarity_cutoff = 0.15){

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
    #MB2 space is being phased out after JAMS ver 2.0.3.
    valid_taxonomic_spaces <- c("Contig_LKT", "ConsolidatedGenomeBin")[c("Contig_LKT", "ConsolidatedGenomeBin") %in% names(propsampleswithtaxspace)[propsampleswithtaxspace > 0.8]]

    retrieve_abundance_table <- function(Sample = NULL, taxonomic_space = NULL, colsToIgnore = NULL){
        abundance_df <- list.data[[paste(Sample, "abundances", sep = "_")]]$taxonomic[[taxonomic_space]][ , ]
        abundance_df <- abundance_df[ , which(!colnames(abundance_df) %in% colsToIgnore)]

        return(abundance_df)
    }

    for (taxonomic_space in valid_taxonomic_spaces){

        #Start by making taxonomic space SummarizedExperiment objects
        flog.info(paste("Making", taxonomic_space, "SummarizedExperiment"))
        LKTdoses <- lapply(Samples, function (x) { retrieve_abundance_table(Sample = x, taxonomic_space = taxonomic_space, colsToIgnore = c("Genome_Size", "PPM", "Genome_Size", "Total_Contigs", "Contig_N50", "Max_Contig_Length", "Quality", "GC_Content", "Total_Coding_Sequences", "Coding_Density", "Num_assemblies_in_taxid", "Num_isolate", "Proportion_of_MAG_in_taxid")) } )
        names(LKTdoses) <- Samples

        #Sorry about the for loop, but it is safer like so.
        for (SN in names(LKTdoses)){
            #Add PPM estimate
            LKTdoses[[SN]]$PPM <- round(((LKTdoses[[SN]]$NumBases / sum(LKTdoses[[SN]]$NumBases)) * 1E6), 0)
            rownames(LKTdoses[[SN]]) <- 1:nrow(LKTdoses[[SN]])
        }

        LKTdosesall <- dplyr::bind_rows(LKTdoses, .id = "Sample")

        #If taxonomic space is ConsolidatedGenomeBin, evaluate, contextualize and rename redundant taxonomies using functional annotation.
        if (taxonomic_space == "ConsolidatedGenomeBin"){
            LKTdosesall <- contextualize_taxonomy(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = functional_contextualization_dissimilarity_cutoff)

            #Bequeath to opt for using later when building functional experiments
            #Keep this here for the time being. make_SummarizedExperiments does not return opt, so at a later date I might add this to the SEobj itself, depending on the object size.
            opt$CGB2LKTdict <- LKTdosesall[ , c("Sample", "ConsolidatedGenomeBin", "LKT"), drop = FALSE] %>% dplyr::mutate(across(where(is.character), as.factor))

            #Fix list.data objects with contextualized taxonomy

            flog.info("Writing contextualized taxonomy to functional data counts")
            for (S2C in as.character(unique(opt$CGB2LKTdict$Sample))){
                curr_CGB2LKT_df <- opt$CGB2LKTdict[which(opt$CGB2LKTdict$Sample == S2C), , drop = FALSE]
                curr_CGB2LKT_df <- curr_CGB2LKT_df %>% dplyr::mutate(across(where(is.factor), as.character))
                #Fix featuredata first
                list.data[[paste(S2C, "featuredata", sep = "_")]]$LKT <- NA
                list.data[[paste(S2C, "featuredata", sep = "_")]]$LKT <- curr_CGB2LKT_df$LKT[match(list.data[[paste(S2C, "featuredata", sep = "_")]]$ConsolidatedGenomeBin, curr_CGB2LKT_df$ConsolidatedGenomeBin)]

                #Now, fix functional abundance
                #Subset curr_CGB2LKT_df to only CGBs present in the functional abundance dataframe. Some may be missing because JAMSalpha only registers functional stratification above a certain small (albeit non-zero) relative abundance, in order to not inflate data frame size with noise.
                curr_CGB2LKT_df <- curr_CGB2LKT_df[curr_CGB2LKT_df$ConsolidatedGenomeBin %in% colnames(list.data[[paste(S2C, "abundances", sep = "_")]][["functional"]][["ConsolidatedGenomeBin"]]), , drop = FALSE]
                colpos <- match(curr_CGB2LKT_df$ConsolidatedGenomeBin, colnames(list.data[[paste(S2C, "abundances", sep = "_")]][["functional"]][["ConsolidatedGenomeBin"]]))
                colnames(list.data[[paste(S2C, "abundances", sep = "_")]][["functional"]][["ConsolidatedGenomeBin"]])[colpos] <- curr_CGB2LKT_df$LKT
                curr_CGB2LKT_df <- NULL
                colpos <- NULL
            }

        } else {
            LKTdosesall$LKT <- LKTdosesall[ , taxonomic_space]
        }

        #Make counts table
        #n.b. In the case of ConsolidatedGenomeBin, there may be still two different MAGs (from MetaBAT2 bins) bearing the same LKT name even after functional clusterization. I have found that these are the exceedingly rare, occurring at a rate of about 0.005% of all Sample LKT pairs. Why this is so, it is not known, but MetaBAT is binning near-identical entities within the same sample into two different bins. It is thus necessary to sum the NumBases and genome completenesses of these. This will have minimal effect on an analysis and genome completeness overshoot can still be ascertained via GenomeCompleteness matrices in SEobj.

        if (taxonomic_space == "ConsolidatedGenomeBin"){
            LKTcountsall <- LKTdosesall[ , c("Sample", "LKT", "NumBases", "Completeness", "Contamination")] %>% group_by(Sample, LKT) %>% summarize(NumBases = sum(NumBases), Completeness = sum(Completeness), Contamination = sum(Contamination), .groups = "drop") 
        } else {
            LKTcountsall <- LKTdosesall[ , c("Sample", "LKT", "NumBases", "Completeness", "Contamination")]
        }

        cts <- LKTcountsall[ , c("Sample", "LKT", "NumBases")] %>% group_by(Sample, LKT) %>% tidyr::pivot_wider(names_from = Sample, values_from = NumBases, values_fill = 0)
        cts <- as.data.frame(cts)
        rownames(cts) <- cts$LKT
        cts$LKT <- NULL
        cts <- cts[order(rowSums(cts), decreasing = TRUE), ]
        featureorder <- rownames(cts)

        #Make tax table
        if (taxonomic_space == "ConsolidatedGenomeBin"){
            taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT")]
            tt <- LKTdosesall[ , taxlvlspresent]
            tt <- tt[!(duplicated(tt)), ]
            #Extract taxid from LKT. Some of the contextualized infraspecies may have moved up to species.
            tt$Taxid <- sapply(tt$LKT, function(x) { extract_NCBI_taxid_from_featname(Taxon = x) } )
        } else {
            taxlvlspresent <- colnames(LKTdosesall)[colnames(LKTdosesall) %in% c("Taxid", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "LKT", "NCBI_taxonomic_rank")]
            tt <- LKTdosesall[ , taxlvlspresent]
            tt <- tt[!(duplicated(tt)), ] 
        }

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
        rownames(tt) <- tt[ , "LKT"]

        sampleorder <- rownames(pheno2)
        tt <- tt[featureorder, ]
        cts <- cts[, sampleorder]

        #Register the total number of NAHS bases sequenced for each sample
        TotBasesSamples <- colSums(cts)
        TotalBasesSequenced <- t(as.matrix(TotBasesSamples))
        rownames(TotalBasesSequenced) <- "NumBases"

        #Get genome completeness matrix
        LKTallGenComp <- LKTcountsall[, c("Sample", "LKT", "Completeness")]
        #Be sure that data is not empty or redundant
        LKTallGenComp[is.na(LKTallGenComp)] <- 0
        LKTallGenCompcts <- LKTallGenComp %>% tidyr::pivot_wider(names_from = Sample, values_from = Completeness, values_fill = 0)
        LKTallGenCompcts <- as.data.frame(LKTallGenCompcts)
        rownames(LKTallGenCompcts) <- LKTallGenCompcts$LKT
        LKTallGenCompcts$LKT <- NULL
        LKTallGenCompcts <- LKTallGenCompcts[featureorder, sampleorder]

        #Get genome contamination matrix
        LKTallGenCont <- LKTcountsall[, c("Sample", "LKT", "Contamination")]
        #Be sure that data is not empty or redundant
        LKTallGenCont[is.na(LKTallGenCont)] <- 0
        LKTallGenContcts <- LKTallGenCont %>% group_by(Sample, LKT) %>% tidyr::pivot_wider(names_from = Sample, values_from = Contamination, values_fill = 0)
        LKTallGenContcts <- as.data.frame(LKTallGenContcts)
        rownames(LKTallGenContcts) <- LKTallGenContcts$LKT
        LKTallGenContcts$LKT <- NULL 
        LKTallGenContcts <- LKTallGenContcts[featureorder, sampleorder]

        assays <- list(as.matrix(cts), as.matrix(LKTallGenCompcts), as.matrix(LKTallGenContcts))
        names(assays) <- c("BaseCounts", "GenomeCompleteness", "GenomeContamination")

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
        analysislencts <- NULL

        if (!stratify_functions_by_taxon){
            #If not stratifying functions by taxa, get functional information Contig_LKT abundances
            valid_taxonomic_space <- "Contig_LKT"
        } else {
            #Declare empty lists to hold data
            master_sparse_featcounts_index_longform <- NULL
            master_sparse_taxon_to_space_basecounts_df <- NULL
            master_sparse_taxon_to_space_genecounts_df <- NULL

            #When stratifying, prioritize Consolidated Genome bin if present, if not, default to Contig_LKT
            valid_taxonomic_space <- c("ConsolidatedGenomeBin", "Contig_LKT")[c("ConsolidatedGenomeBin", "Contig_LKT") %in% valid_taxonomic_spaces][1]
        }

        if (valid_taxonomic_space != "ConsolidatedGenomeBin"){
            appropriate_featurdata_tax_colm <- valid_taxonomic_space
        } else {
            appropriate_featurdata_tax_colm <- "LKT"
        }

        batch_start_time <- Sys.time()
        for (sampnum in 1:length(Samples)){
            SN <- Samples[sampnum]
            curr_FD <- list.data[[paste(SN, "abundances", sep = "_")]]$functional[[valid_taxonomic_space]]
            #Subset to current analysis
            curr_FD <- curr_FD[which(curr_FD$Analysis == analysis), , drop = FALSE]
            #If not stratifying, eliminate taxon columns at this point
            if (!stratify_functions_by_taxon){
                curr_FD <- curr_FD[ , c("Accession", "Description", "NumBases")]
            }

            #Check that there is any data for that particular sample - there may not be.
            if (nrow(curr_FD) > 0){
                #There shouldn't be any NAs but if there are make them 0
                curr_FD[is.na(curr_FD)] <- 0
                Taxoncols <- colnames(curr_FD)[5:ncol(curr_FD)]
                curr_FD$Sample <- SN

                #Deal with gene number tallies
                curr_featdf <- list.data[[paste(SN, "featuredata", sep = "_")]][ , c("Feature", "LengthDNA", analysis, appropriate_featurdata_tax_colm)]
                colnames(curr_featdf)[which(colnames(curr_featdf) == analysis)] <- "Accession"
                colnames(curr_featdf)[which(colnames(curr_featdf) == appropriate_featurdata_tax_colm)] <- "Taxon"

                #Fix "none" to analysis_none
                curr_featdf$Accession[which(curr_featdf$Accession == "none")] <- paste(analysis, "none", sep = "_")
                #Ensure numeric
                curr_featdf$LengthDNA <- as.numeric(curr_featdf$LengthDNA)

                #Separate rows if not Product or ECNumber, as there may be more than a single accession for a single gene separated by a "|".
                if (!(analysis %in% c("Product", "ECNumber"))){
                    curr_featdf <- curr_featdf %>% tidyr::separate_rows(all_of("Accession"), sep = fixed("\\|"))
                }

                if (stratify_functions_by_taxon){
                    curr_SM <- Matrix::Matrix(data = as.matrix(curr_FD[ , Taxoncols]), sparse = TRUE)
                    colnames(curr_SM) <- Taxoncols
                    #In order to conserve RAM and the size of the taxonomic stratification sparse matrix, do not include any taxon whose total sum of bases across the entire metagenomic proteome is smaller than 10000 bp (10 kbp). 10 kbp / ~3 Mbp is a completeness of 0.3%, with a sequencing depth of 1X if the entire sample consisted of a single species (which it doesn't). So this filtration is very conseravative and will shrink memory usage considerably.
                    curr_taxaToMerge <- names(which(Matrix::colSums(curr_SM) < 10000))
                    if (length(curr_taxaToMerge) > 0){
                        #Sum all low abundance values into a single column.
                        ULA_sums <- Matrix::rowSums(curr_SM[ , curr_taxaToMerge, drop = FALSE])
                        compl_SM <- Matrix::Matrix(data = as.matrix(ULA_sums), sparse = TRUE)
                        colnames(compl_SM) <- "Ultra_low_abundance_LKTs"
                        #Eliminate ULAs from matrix
                        curr_SM <- curr_SM[ , !(colnames(curr_SM) %in% curr_taxaToMerge), drop = FALSE]
                        #Add aggregated column
                        curr_SM <- cbind(curr_SM, compl_SM)
                    }
                    #Make row pairs temporarily named as sample-feature pairs to avoid double entries when merging sparse matrix with the next sample. For instance Urease§Sample1 should be a different row to Urease§Sample17, etc. Sparse matrix merging will preserve only column names for LKTs which are shared across samples.
                    rownames(curr_SM) <- paste(SN, rownames(curr_SM), sep = "§")

                    #Create a temporary long-form index
                    curr_sparse_featcounts_index <- curr_FD[ , c("Sample", "Accession")]
                    rownames(curr_sparse_featcounts_index) <- paste(SN, curr_sparse_featcounts_index$Accession, sep = "§")
                    #Ensure row order is identical
                    curr_sparse_featcounts_index <- curr_sparse_featcounts_index[rownames(curr_SM), ]
                    #Add matrix row number to temporary index, counting from the last row of the previous matrix
                    if (!is.null(master_sparse_taxon_to_space_basecounts_df)){
                        #A previous sample was accrued so, increment the row number from the current size of the master matrix
                        curr_sparse_featcounts_index$RowNumber <- nrow(master_sparse_taxon_to_space_basecounts_df) + (1:nrow(curr_SM))
                    } else {
                        #This is the first sample
                        curr_sparse_featcounts_index$RowNumber <- 1:nrow(curr_SM)
                    }
                    #Commit index to long form dataframe
                    master_sparse_featcounts_index_longform <- rbind(master_sparse_featcounts_index_longform, curr_sparse_featcounts_index)
                    #Commit rownames as index number
                    rownames(master_sparse_featcounts_index_longform) <- master_sparse_featcounts_index_longform$RowNumber

                    #Merge curr_SM to the master matrix
                    master_sparse_taxon_to_space_basecounts_df <- merge_sparse_matrix(matlist = list(master_sparse_taxon_to_space_basecounts_df, curr_SM))
                    #Reset row names to numerical to conserve RAM
                    rownames(master_sparse_taxon_to_space_basecounts_df) <- 1:nrow(master_sparse_taxon_to_space_basecounts_df)

                    #Deal with gene number stratifications
                    #Pivot to wide form
                    curr_split_numgenes <- as.data.frame(curr_featdf %>% group_by(across(all_of(c("Taxon", "Accession")))) %>% summarise(NumGenes = length(Accession), .groups = "keep") %>% pivot_wider(names_from = Taxon, values_from = NumGenes, values_fill = 0))
                    rownames(curr_split_numgenes) <- paste(SN, curr_split_numgenes$Accession, sep = "§")
                    curr_split_numgenes$Accession <- NULL

                    #Curtail to the same features as for basecounts and make sparse
                    curr_split_numgenes <- Matrix::Matrix(data = as.matrix(curr_split_numgenes[rownames(curr_sparse_featcounts_index), ]), sparse = TRUE)

                    #Aggregate low abundance taxa, if applicable.
                    if (length(curr_taxaToMerge) > 0){
                        #Sum all low abundance values into a single column.
                        ULA_gene_sums <- Matrix::rowSums(curr_split_numgenes[ , curr_taxaToMerge, drop = FALSE])
                        compl_split_numgenes <- Matrix::Matrix(data = as.matrix(ULA_gene_sums), sparse = TRUE)
                        colnames(compl_split_numgenes) <- "Ultra_low_abundance_LKTs"
                        #Eliminate ULAs from matrix
                        curr_split_numgenes <- curr_split_numgenes[ , !(colnames(curr_split_numgenes) %in% curr_taxaToMerge), drop = FALSE]
                        #Add aggregated column
                        curr_split_numgenes <- cbind(curr_split_numgenes, compl_split_numgenes)
                    }
                    #Ensure exact row and column order
                    curr_split_numgenes <- curr_split_numgenes[rownames(curr_sparse_featcounts_index), colnames(curr_SM)]

                    #Merge sparse version in to master gene count matrix
                    master_sparse_taxon_to_space_genecounts_df <- merge_sparse_matrix(matlist = list(master_sparse_taxon_to_space_genecounts_df, curr_split_numgenes))

                    #Clean up
                    curr_SM <- NULL
                    curr_split_numgenes <- NULL
                    curr_taxaToMerge <- NULL
                    curr_sparse_featcounts_index <- NULL
                    gc()
                }

                #Eliminate Taxoncols to save RAM
                curr_FD$Sample <- SN
                rownames(curr_FD) <- paste(SN, rownames(curr_FD), sep = "§")
                analysisdoses <- rbind(analysisdoses, curr_FD[ , c("Sample", "Accession", "Description", "NumBases")])
                curr_featdf$Sample <- SN
                curr_featdf$Taxon <- NULL
                curr_analysisgencts <- as.data.frame(curr_featdf %>% group_by(Accession) %>% summarise(NumGenes = length(Accession), .groups = "keep"))
                curr_analysislencts <- as.data.frame(curr_featdf %>% group_by(Accession) %>% summarise(LenGenes = sum(LengthDNA), .groups = "keep"))
                curr_analysisgencts$Sample <- SN
                curr_analysislencts$Sample <- SN
                rownames(curr_analysisgencts) <- paste(SN, curr_analysisgencts$Accession, sep = "§")
                rownames(curr_analysislencts) <- paste(SN, curr_analysislencts$Accession, sep = "§")
                curr_analysisgencts <- curr_analysisgencts[rownames(curr_FD), c("Sample", "Accession", "NumGenes"), drop = FALSE]
                curr_analysislencts <- curr_analysislencts[rownames(curr_FD), c("Sample", "Accession", "LenGenes"), drop = FALSE]
                analysisgencts <- rbind(analysisgencts, curr_analysisgencts)
                analysislencts <- rbind(analysislencts, curr_analysislencts)

                #Clean up
                curr_FD <- NULL
                curr_featdf <- NULL
                curr_analysisgencts <- NULL
                curr_analysislencts <- NULL
                gc()
            } #End conditional for there being data for that analysis for that sample

            #Give message on progress
            if (sampnum == length(Samples)) {
                elapsed <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
                avg_per_sample <- elapsed / sampnum
                flog.info(sprintf("Completed %d/%d samples | elapsed time: %s", sampnum, length(Samples), pretty_time(elapsed)))
            } else if (length(Samples) > 20 && sampnum %% 20 == 0) {
                elapsed <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
                avg_per_sample <- elapsed / sampnum
                eta_secs <- avg_per_sample * (length(Samples) - sampnum)
                flog.info(sprintf("Completed %d/%d samples | elapsed time: %s | ETA: %s", sampnum, length(Samples), pretty_time(elapsed), pretty_time(eta_secs)))
            }

        } #End loop for obtaining data for each sample

        #Obtain unstratified counts matrix
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

        ################################
        ## Make GeneCounts counts table
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

        #################################
        ## Make GeneLengths counts table
        genlencts <- NULL

        #Fix for not getting duplicate accessions. Delete the Description column to get clean pivot from long to wide.
        genlencts <- analysislencts[ , c("Sample", "Accession", "LenGenes")] %>% pivot_wider(names_from = Sample, values_from = LenGenes, values_fill = 0)
        genlencts <- as.data.frame(genlencts)
        genlencts[is.na(genlencts)] <- 0
        rownames(genlencts) <- genlencts$Accession
        genlencts$Accession <- NULL

        #Deal with samples with NO results for this analysis
        emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(genlencts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(genlencts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(genlencts)
            genlencts <- cbind(genlencts, complementarycts)
        }
        genlencts <- as.matrix(genlencts)


        ######################
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

        #Order gennumcts and genlencts
        gennumcts <- gennumcts[featureorder, sampleorder]
        genlencts <- genlencts[featureorder, sampleorder]

        assays <- list()
        assays$BaseCounts <- cts
        assays$GeneCounts <- gennumcts
        assays$GeneLengths <- genlencts
        #Add stratification index as assay if applicable
        if (stratify_functions_by_taxon){
            master_sparse_featcounts_index <- master_sparse_featcounts_index_longform %>% group_by(Sample, Accession) %>% tidyr::pivot_wider(names_from = Sample, values_from = RowNumber, values_fill = NA)
            master_sparse_featcounts_index <- as.data.frame(master_sparse_featcounts_index)
            rownames(master_sparse_featcounts_index) <- master_sparse_featcounts_index$Accession
            master_sparse_featcounts_index$Accession <- NULL
            #Add empty samples if necessary
            emptySamples <- rownames(phenoanal)[!(rownames(phenoanal) %in% colnames(master_sparse_featcounts_index))]

            if (length(emptySamples) > 0){
                complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(master_sparse_featcounts_index), data = 0)
                complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
                colnames(complementarycts) <- emptySamples
                rownames(complementarycts) <- rownames(master_sparse_featcounts_index)
                master_sparse_featcounts_index <- cbind(master_sparse_featcounts_index, complementarycts)
            }

            master_sparse_featcounts_index <- master_sparse_featcounts_index[rownames(cts), colnames(cts)]
            assays$SparseIndex <- as.matrix(master_sparse_featcounts_index)
        }

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
 
            #Add features-by-taxon matrix and index into SummarizedExperiment object
            metadata(SEobj)$allfeaturesbytaxa_matrix <- master_sparse_taxon_to_space_basecounts_df
            metadata(SEobj)$allfeaturesbytaxa_GeneCounts_matrix <- master_sparse_taxon_to_space_genecounts_df

            #Clean up
            master_sparse_featcounts_index <- NULL
            master_sparse_taxon_to_space_basecounts_df <- NULL
            master_sparse_taxon_to_space_genecounts_df <- NULL
            master_sparse_featcounts_index_longform <- NULL
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
