#' make_SummarizedExperiments(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8, totmembytes = NULL, stratification_threads = NULL, functional_contextualization_dissimilarity_cutoff = 0.15)
#'
#' Makes a SummarizedExperiment object for every analysis that is possible to make given loaded jams files in list.data.
#'
#' @param stratification_threads Number of forked workers used for the per-sample pass that
#'   builds the functional (non-taxonomic) SummarizedExperiment objects. Defaults to
#'   min(threads, 16), which is deliberately conservative: each worker inherits list.data
#'   copy-on-write for free, but also holds its own transient per-sample objects. Raise it if
#'   there is RAM headroom; set it to 1 to force the old single-threaded behaviour. Ignored
#'   (forced to 1) on platforms without fork.
#' @export

make_SummarizedExperiments <- function(pheno = NULL, onlysamples = NULL, onlyanalyses = c("LKT", "Product", "ECNumber", "Pfam", "Interpro", "resfinder", "PRINTS", "GO"), minnumsampanalysis = NULL, minpropsampanalysis = 0.1, stratify_functions_by_taxon = TRUE, add_TNF_data = FALSE, list.data = NULL, phenolabels = NULL, cdict = cdict, threads = 8, totmembytes = NULL, stratification_threads = NULL, functional_contextualization_dissimilarity_cutoff = 0.15){

    require(SummarizedExperiment)
    require(Matrix)
    data(MetaCycAccession2Description)
    data(ECdescmap)
    data(GOtermdict)
    data(InterproDict)
    data(blast_lookup)

    if (threads > 1){
        #Assure necessary packages are loaded
        library(foreach)
        library(doParallel)
    }

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
    #Contig_LKT space is being phased out after JAMS ver 2.2.6
    valid_taxonomic_spaces <- c("ConsolidatedGenomeBin")[c("ConsolidatedGenomeBin") %in% names(propsampleswithtaxspace)[propsampleswithtaxspace > 0.8]]

    #############################################################################
    ## Stratification consistency audit.
    ## After contextualize_taxonomy has rewritten LKT names and CGB2LKTdict has been
    ## built, this checks that every ConsolidatedGenomeBin the dict intends to rewrite
    ## can actually be found (a) as a column in the sample's functional
    ## ConsolidatedGenomeBin abundance table and (b) in the sample's featuredata
    ## ConsolidatedGenomeBin column. Mismatches here would silently desynchronise the
    ## -u taxon-stratified functional matrices from the taxonomic SEobj, so we surface
    ## them precisely rather than letting a downstream surprise occur. A clean run logs
    ## a single reassuring line. This is cheap: it is O(number of CGBs), not O(features).
    #############################################################################
    audit_stratification_consistency <- function(CGB2LKTdict = NULL, list.data = NULL){
        n_checked <- 0
        missing_from_featuredata <- character(0)
        missing_from_functional <- character(0)

        for (S2C in as.character(unique(CGB2LKTdict$Sample))){
            curr_dict <- CGB2LKTdict[which(CGB2LKTdict$Sample == S2C), , drop = FALSE]
            curr_dict <- curr_dict %>% dplyr::mutate(across(where(is.factor), as.character))

            curr_featuredata <- list.data[[paste(S2C, "featuredata", sep = "_")]]
            curr_func <- list.data[[paste(S2C, "abundances", sep = "_")]][["functional"]][["ConsolidatedGenomeBin"]]

            featuredata_CGBs <- if (!is.null(curr_featuredata)) unique(curr_featuredata$ConsolidatedGenomeBin) else character(0)
            functional_CGBs <- if (!is.null(curr_func)) colnames(curr_func) else character(0)

            for (cgb in unique(curr_dict$ConsolidatedGenomeBin)){
                n_checked <- n_checked + 1
                if (!(cgb %in% featuredata_CGBs)){
                    missing_from_featuredata <- c(missing_from_featuredata, paste0("Sample=", S2C, " | CGB=", cgb))
                }
                #Note: it is legitimate and expected for a CGB to be absent from the functional
                #ConsolidatedGenomeBin abundance table, because JAMSalpha only records functional
                #stratification above a small (non-zero) relative abundance. So absence here is NOT
                #flagged as an error; it is tallied separately for information only.
                if (!(cgb %in% functional_CGBs)){
                    missing_from_functional <- c(missing_from_functional, paste0("Sample=", S2C, " | CGB=", cgb))
                }
            }
        }

        if (length(missing_from_featuredata) == 0){
            flog.info(paste0("make_SummarizedExperiments stratification audit: all ", n_checked, " ConsolidatedGenomeBins in the contextualized dictionary reconcile with their sample featuredata. Taxonomic and functional (-u) domains are consistent."))
        } else {
            flog.warn(paste0("make_SummarizedExperiments stratification audit: ", length(missing_from_featuredata), " of ", n_checked, " ConsolidatedGenomeBins in the contextualized dictionary are NOT present in their sample featuredata ConsolidatedGenomeBin column. This may desynchronise the -u taxon-stratified functional matrices. Offending bins follow:"))
            max_report <- 50
            to_report <- missing_from_featuredata[1:min(length(missing_from_featuredata), max_report)]
            for (msg in to_report){
                flog.warn(paste0("   ", msg))
            }
            if (length(missing_from_featuredata) > max_report){
                flog.warn(paste0("   ... and ", (length(missing_from_featuredata) - max_report), " more not shown."))
            }
        }

        #Informational only: how many dict CGBs are below the functional-stratification recording
        #threshold in JAMSalpha. This is expected and benign, reported at info level only if present.
        if (length(missing_from_functional) > 0){
            flog.info(paste0("make_SummarizedExperiments stratification audit: ", length(missing_from_functional), " of ", n_checked, " ConsolidatedGenomeBins are absent from the functional ConsolidatedGenomeBin abundance table (expected for bins below the JAMSalpha functional-stratification abundance floor; not an error)."))
        }

        return(invisible(NULL))
    }

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

        #Check whether any is NULL
        null_doses <- sapply(LKTdoses, is.null)
        if (any(null_doses)){
            null_tax_samples <- names(which(null_doses))
            flog.warn(paste("WARNING: Could not find counts data for taxonomic space", taxonomic_space, "for samples:", paste0(null_tax_samples, collapse = ", ")))
            flog.warn(paste("Eliminating", length(null_tax_samples),"sample(s) from", taxonomic_space, "analysis"))
            LKTdoses <- LKTdoses[names(LKTdoses)[!names(LKTdoses) %in% null_tax_samples]]
        }

        #Sorry about the for loop, but it is safer like so.
        for (SN in names(LKTdoses)){
            #Add PPM estimate
            LKTdoses[[SN]]$PPM <- round(((LKTdoses[[SN]]$NumBases / sum(LKTdoses[[SN]]$NumBases)) * 1E6), 0)
            rownames(LKTdoses[[SN]]) <- 1:nrow(LKTdoses[[SN]])
        }

        LKTdosesall <- dplyr::bind_rows(LKTdoses, .id = "Sample")

        #If taxonomic space is ConsolidatedGenomeBin, evaluate, contextualize and rename redundant taxonomies using functional annotation.
        if (taxonomic_space == "ConsolidatedGenomeBin"){
            #If isolates, the "ConsolidatedGenomeBin" column may be missing.
            if (("isolate" %in% opt$phenotable$JAMS_Run_type) && (!("ConsolidatedGenomeBin" %in% colnames(LKTdosesall)))){
                LKTdosesall$ConsolidatedGenomeBin <- LKTdosesall$LKT
                #Switch LKT tag to CGB. For isolates, the CGB IS the LKT.
                LKTdosesall$ConsolidatedGenomeBin <- gsub("^LKT__", "CGB__", LKTdosesall$ConsolidatedGenomeBin)
            }

            LKTdosesall <- contextualize_taxonomy(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = functional_contextualization_dissimilarity_cutoff, threads = threads, totmembytes = totmembytes)

            #Bequeath to opt for using later when building functional experiments
            #Keep this here for the time being. make_SummarizedExperiments does not return opt, so at a later date I might add this to the SEobj itself, depending on the object size.
            opt$CGB2LKTdict <- LKTdosesall[ , c("Sample", "ConsolidatedGenomeBin", "LKT"), drop = FALSE] %>% dplyr::mutate(across(where(is.character), as.factor))

            #Audit consistency between the contextualized dictionary, featuredata, and functional
            #abundance tables before we rely on them for -u stratification.
            audit_stratification_consistency(CGB2LKTdict = opt$CGB2LKTdict %>% dplyr::mutate(across(where(is.factor), as.character)), list.data = list.data)

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
        Taxid2gram <- Taxid2gram[!(duplicated(Taxid2gram$Taxid)), ]
        tt <- left_join(as.data.frame(tt), Taxid2gram, by = "Taxid")

        #Preserve all four canonical Gram categories from JAMStaxtable (positive, negative,
        #N_A = does not apply, e.g. eukaryotes/viruses; undetermined = unknown). Join misses
        #(dark matter with no taxid match) are "undetermined", since we genuinely do not know
        #their Gram status, whereas N_A is a positive assertion that Gram does not apply.
        canonical_gram <- c("positive", "negative", "N_A", "undetermined")
        tt$Gram[which(is.na(tt$Gram))] <- "undetermined"
        tt$Gram[!(tt$Gram %in% canonical_gram)] <- "undetermined"
        tt <- tt[ , c("Gram", taxlvlspresent)]

        #The taxonomic table must have exactly ONE row per LKT, because LKT is the row key of
        #the SummarizedExperiment (counts were already collapsed by Sample x LKT above). Two
        #entities can legitimately share an LKT name while carrying different upstream lineage
        #strings; the typical case is dark matter: a cLKT__Unclassified bin (taxid 0) and an
        #unanchored MetaBAT2 bin (classified only to root/cellular organisms) both collapse to
        #LKT__Unclassified but may differ at Domain/Kingdom. Previously this crashed JAMSbeta with
        #"duplicate 'row.names' are not allowed". Keep the first lineage seen and say so in the log.
        if (any(duplicated(tt$LKT))){
            dupe_LKTs <- unique(tt$LKT[duplicated(tt$LKT)])
            flog.warn(paste0("make_SummarizedExperiments: ", length(dupe_LKTs), " LKT name(s) in the ", taxonomic_space, " space map to more than one taxonomic lineage. Keeping the first lineage for each. LKTs: ", paste0(dupe_LKTs, collapse = ", ")))
            tt <- tt[!(duplicated(tt$LKT)), , drop = FALSE]
        }
        rownames(tt) <- tt[ , "LKT"]

        sampleorder <- rownames(pheno2)
        tt <- tt[featureorder, ]
        cts <- cts[, sampleorder]

        #Save CGB taxonomic table for including in metadata of taxon-stratified functional SEobj.
        if (taxonomic_space == "ConsolidatedGenomeBin"){
            CGB_tt <- tt
        }

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

    #########################################################################################
    ## HOW THE FUNCTIONAL PHASE WORKS (refactored for large, cross-cohort sample numbers).
    ##
    ## It used to be a single threaded loop over samples which, on EVERY iteration, rbind()ed
    ## onto four growing data frames and called merge_sparse_matrix() twice to widen and extend
    ## two growing sparse matrices. Each of those accumulations copies the whole accumulated
    ## object, so the total cost grew with the SQUARE of the number of samples. A 1409 sample
    ## cohort projected to roughly 26 hours for one analysis. Two full gc() calls per sample made
    ## it worse still, because the cost of a full garbage collection is set by the size of the
    ## live heap (list.data, hundreds of GB for a cohort this size), not by the garbage that the
    ## iteration actually produced.
    ##
    ## It now works in two stages:
    ##
    ##   STAGE 1  Every sample is reduced INDEPENDENTLY by build_sample_stratification_block()
    ##            into a compact per-sample contribution. Independent means it can be run with
    ##            parallel::mclapply. Fork is used deliberately: workers inherit list.data
    ##            copy-on-write, so not one byte of it is serialised to them, which is the same
    ##            reasoning behind the contextualize_taxonomy refactor. Workers never call gc()
    ##            and never log, so they neither disturb the inherited heap nor interleave output.
    ##
    ##   STAGE 2  The master stacks all per-sample contributions ONCE. The sparse blocks go
    ##            through rbind_sparse_blocks() (a single linear pass, no repeated copying) and
    ##            the tabular parts are concatenated and written into their matrices by direct
    ##            index assignment instead of pivot_wider().
    ##
    ## Row numbering is unchanged. Samples are visited in the order of `Samples` and each
    ## contributes a contiguous band of row numbers, exactly as the old serial loop did, so the
    ## result is deterministic and does not depend on how many threads were used.
    #########################################################################################

    #Resolve how many forked workers Stage 1 may use. Each worker holds only its own transient
    #per-sample objects on top of the shared, inherited list.data, but with very wide featuredata
    #that still adds up, so the default is deliberately conservative and can be raised.
    fork_is_available <- (.Platform$OS.type == "unix")
    if (is.null(stratification_threads)){
        stratification_threads <- min(threads, 16)
    }
    stratification_threads <- suppressWarnings(as.integer(stratification_threads))
    if (is.na(stratification_threads) || (stratification_threads < 1)){
        stratification_threads <- 1
    }
    stratification_threads <- min(stratification_threads, max(1, as.integer(threads)))
    if (!fork_is_available && (stratification_threads > 1)){
        flog.warn("make_SummarizedExperiments: parallel::mclapply needs a unix-alike (fork) in order to share list.data with its workers without copying it. This does not look like one, so the per-sample functional pass will run single threaded.")
        stratification_threads <- 1
    }
    flog.info(paste0("make_SummarizedExperiments: the per-sample functional pass will use ", stratification_threads, " of the ", threads, " CPU(s) available. Raise this with the stratification_threads argument (JAMSbeta: --stratification_threads) if there is RAM headroom, or set it to 1 to force the old serial behaviour."))

    #How many samples are handed out per mclapply call. Larger means fewer fork() calls (each
    #fork of a master holding hundreds of GB has to copy its page tables, which is not free);
    #smaller means more frequent progress lines and a tighter bound on peak memory.
    stratification_chunk_size <- max(1L, (stratification_threads * 4L))

    for (analnumb in 1:length(possibleanalyses)) {
        analysis <- possibleanalyses[analnumb]
        flog.info(paste("Making", analysis, "SummarizedExperiment"))

        phenoanal <- pheno2
        sampleorder <- rownames(phenoanal)

        if (!stratify_functions_by_taxon){
            #If not stratifying functions by taxa, get functional information Contig_LKT abundances
            valid_taxonomic_space <- "Contig_LKT"
        } else {
            #When stratifying, prioritize Consolidated Genome bin if present, if not, default to Contig_LKT
            valid_taxonomic_space <- c("ConsolidatedGenomeBin", "Contig_LKT")[c("ConsolidatedGenomeBin", "Contig_LKT") %in% valid_taxonomic_spaces][1]
        }

        if (valid_taxonomic_space != "ConsolidatedGenomeBin"){
            appropriate_featurdata_tax_colm <- valid_taxonomic_space
        } else {
            appropriate_featurdata_tax_colm <- "LKT"
        }

        #####################################################################################
        ## STAGE 1: reduce every sample, in parallel where possible.
        #####################################################################################
        reduce_one_sample <- function(SN){
            #A single odd sample must never take the whole JAMSbeta run down with it, so the
            #error is trapped here and reported by the master afterwards.
            tryCatch({
                build_sample_stratification_block(Sample = SN, list.data = list.data, analysis = analysis, valid_taxonomic_space = valid_taxonomic_space, appropriate_featurdata_tax_colm = appropriate_featurdata_tax_colm, stratify_functions_by_taxon = stratify_functions_by_taxon)
            }, error = function(e){
                list(Sample = SN, has_data = FALSE, failed = TRUE, error = conditionMessage(e), messages = character(0), accessions = character(0), descriptions = character(0), numbases = numeric(0), numgenes = numeric(0), lengenes = numeric(0), base_block = NULL, gene_block = NULL)
            })
        }

        per_sample_results <- vector("list", length(Samples))
        chunk_starts <- seq(1, length(Samples), by = stratification_chunk_size)
        batch_start_time <- Sys.time()

        for (chunk_start in chunk_starts){
            chunk_end <- min((chunk_start + stratification_chunk_size - 1), length(Samples))
            chunk_idx <- chunk_start:chunk_end

            if (stratification_threads > 1){
                chunk_res <- parallel::mclapply(Samples[chunk_idx], reduce_one_sample, mc.cores = stratification_threads, mc.preschedule = TRUE, mc.cleanup = TRUE)
            } else {
                chunk_res <- lapply(Samples[chunk_idx], reduce_one_sample)
            }

            for (k in seq_along(chunk_idx)){
                res <- chunk_res[[k]]
                #A forked child that is killed outright (almost always the OOM killer) hands
                #back a try-error or a NULL rather than our list, so say so plainly.
                if (is.null(res) || inherits(res, "try-error") || !is.list(res) || is.null(res$Sample)){
                    res <- list(Sample = Samples[chunk_idx[k]], has_data = FALSE, failed = TRUE, error = "the parallel worker died without returning a result, which almost always means it ran out of memory. Lower stratification_threads.", messages = character(0), accessions = character(0), descriptions = character(0), numbases = numeric(0), numgenes = numeric(0), lengenes = numeric(0), base_block = NULL, gene_block = NULL)
                }
                per_sample_results[[chunk_idx[k]]] <- res
            }
            chunk_res <- NULL

            #Give message on progress
            elapsed <- as.numeric(difftime(Sys.time(), batch_start_time, units = "secs"))
            if (chunk_end == length(Samples)) {
                flog.info(sprintf("Completed %d/%d samples | elapsed time: %s", chunk_end, length(Samples), pretty_time(elapsed)))
            } else {
                avg_per_sample <- elapsed / chunk_end
                eta_secs <- avg_per_sample * (length(Samples) - chunk_end)
                flog.info(sprintf("Completed %d/%d samples | elapsed time: %s | ETA: %s", chunk_end, length(Samples), pretty_time(elapsed), pretty_time(eta_secs)))
            }
        } #End of the per-sample reduction

        #Replay anything the workers wanted said, in sample order, capped so a systemic problem
        #cannot bury the log.
        worker_messages <- unlist(lapply(per_sample_results, function(x) { x$messages } ), use.names = FALSE)
        if (length(worker_messages) > 0){
            max_report <- 50
            for (msg in worker_messages[1:min(length(worker_messages), max_report)]){
                flog.warn(msg)
            }
            if (length(worker_messages) > max_report){
                flog.warn(paste0("   ... and ", (length(worker_messages) - max_report), " further message(s) of the same kind not shown."))
            }
        }
        worker_messages <- NULL

        failed_flags <- vapply(per_sample_results, function(x) { isTRUE(x$failed) }, logical(1))
        if (any(failed_flags)){
            failed_samples <- Samples[which(failed_flags)]
            flog.warn(paste0("make_SummarizedExperiments: ", length(failed_samples), " sample(s) failed while being reduced for analysis ", analysis, " and will contribute NO data to it. This is not fatal, but the affected samples will read as all zero for this analysis. Samples and errors follow:"))
            max_report <- 20
            for (fs in failed_samples[1:min(length(failed_samples), max_report)]){
                flog.warn(paste0("   Sample=", fs, " | error: ", per_sample_results[[which(Samples == fs)[1]]]$error))
            }
            if (length(failed_samples) > max_report){
                flog.warn(paste0("   ... and ", (length(failed_samples) - max_report), " more not shown."))
            }
        }

        #Silent tally of taxon-stratification coverage across this analysis, so a systemic
        #failure of the -u path is surfaced here rather than passing unnoticed.
        has_data_flags <- vapply(per_sample_results, function(x) { isTRUE(x$has_data) }, logical(1))
        strat_samples_with_data <- sum(has_data_flags)
        strat_samples_without_data <- sum(!has_data_flags & !failed_flags)

        if (stratify_functions_by_taxon){
            if (strat_samples_with_data == 0){
                flog.warn(paste0("make_SummarizedExperiments: for analysis ", analysis, ", NO samples yielded usable taxon-stratified data (", strat_samples_without_data, " samples had no data for this analysis, ", sum(failed_flags), " failed). The -u stratification metadata for this analysis will be empty. Check your inputs if this is unexpected."))
            } else {
                flog.info(paste0("make_SummarizedExperiments: taxon-stratification tally for analysis ", analysis, ": ", strat_samples_with_data, " sample(s) contributed stratified data, ", strat_samples_without_data, " sample(s) had no data for this analysis, ", sum(failed_flags), " sample(s) failed. Consistent."))
            }
        }

        if (strat_samples_with_data == 0){
            flog.warn(paste0("make_SummarizedExperiments: not a single sample has data for analysis ", analysis, ", so no SummarizedExperiment object can be built for it. Skipping this analysis and carrying on with the rest."))
            per_sample_results <- NULL
            gc()
            next
        }

        #####################################################################################
        ## STAGE 2: stack every sample's contribution, once.
        #####################################################################################
        has_idx <- which(has_data_flags)
        nrows_each <- vapply(per_sample_results[has_idx], function(x) { length(x$accessions) }, integer(1))
        cumulative_rows <- cumsum(as.numeric(nrows_each))
        total_rows <- cumulative_rows[length(cumulative_rows)]
        if (total_rows > .Machine$integer.max){
            stop(paste0("make_SummarizedExperiments: analysis ", analysis, " would need ", total_rows, " stratification rows, beyond what a single sparse matrix can address. Split the cohort, or build fewer analyses in one run."))
        }
        total_rows <- as.integer(total_rows)
        #Number of rows preceding each contributing sample. This reproduces, exactly, the
        #RowNumber bookkeeping the old serial loop did as it went.
        row_offsets <- as.integer(c(0, cumulative_rows[-length(cumulative_rows)]))

        analysis_long <- data.frame(
            Sample = rep(Samples[has_idx], times = nrows_each),
            Accession = unlist(lapply(per_sample_results[has_idx], function(x) { x$accessions } ), use.names = FALSE),
            Description = unlist(lapply(per_sample_results[has_idx], function(x) { x$descriptions } ), use.names = FALSE),
            NumBases = unlist(lapply(per_sample_results[has_idx], function(x) { x$numbases } ), use.names = FALSE),
            NumGenes = unlist(lapply(per_sample_results[has_idx], function(x) { x$numgenes } ), use.names = FALSE),
            LenGenes = unlist(lapply(per_sample_results[has_idx], function(x) { x$lengenes } ), use.names = FALSE),
            RowNumber = seq_len(total_rows),
            stringsAsFactors = FALSE)

        #The tabular per-sample payload has been copied into analysis_long, so let it go. The
        #sparse blocks are deliberately kept until rbind_sparse_blocks has consumed them.
        for (k in seq_along(per_sample_results)){
            per_sample_results[[k]][["accessions"]] <- NULL
            per_sample_results[[k]][["descriptions"]] <- NULL
            per_sample_results[[k]][["numbases"]] <- NULL
            per_sample_results[[k]][["numgenes"]] <- NULL
            per_sample_results[[k]][["lengenes"]] <- NULL
        }

        #####################################################################################
        ## Unstratified assay matrices. Built by direct index assignment into a pre-made matrix,
        ## which is both far faster and far lighter than pivot_wider() on tens of millions of
        ## rows, and which makes the old "deal with samples with NO results" blocks unnecessary:
        ## a sample that contributed nothing simply keeps its pre-filled column.
        #####################################################################################
        feature_accessions <- unique(analysis_long$Accession)
        acc_pos <- match(analysis_long$Accession, feature_accessions)
        samp_pos <- match(analysis_long$Sample, sampleorder)
        placeable <- which(!is.na(acc_pos) & !is.na(samp_pos))
        cell_index <- cbind(acc_pos[placeable], samp_pos[placeable])

        #A given accession must appear at most once per sample within an analysis. If that is
        #violated the abundance table is malformed upstream, so say so loudly rather than let
        #one of the values win silently.
        if (length(placeable) > 1){
            cell_id <- (as.numeric(samp_pos[placeable]) - 1) * as.numeric(length(feature_accessions)) + as.numeric(acc_pos[placeable])
            n_dupe_cells <- sum(duplicated(cell_id))
            if (n_dupe_cells > 0){
                flog.warn(paste0("make_SummarizedExperiments: found ", n_dupe_cells, " duplicated Sample/Accession pair(s) in analysis ", analysis, ". Keeping the last occurrence of each. This points at a malformed functional abundance table upstream and is worth investigating."))
            }
            cell_id <- NULL
        }

        make_assay_matrix <- function(values = NULL, fill = 0){
            mat <- matrix(fill, nrow = length(feature_accessions), ncol = length(sampleorder), dimnames = list(feature_accessions, sampleorder))
            mat[cell_index] <- values[placeable]

            return(mat)
        }

        cts <- make_assay_matrix(values = analysis_long$NumBases, fill = 0)
        gennumcts <- make_assay_matrix(values = analysis_long$NumGenes, fill = 0)
        genlencts <- make_assay_matrix(values = analysis_long$LenGenes, fill = 0)

        if (stratify_functions_by_taxon){
            #n.b. the fill is NA, NOT 0. retrieve_features_by_taxa drops NA row numbers, but a 0
            #would survive that filter and then fail to index the stratification matrix. The old
            #code filled columns for empty samples with 0, which was a latent crash.
            master_sparse_featcounts_index <- matrix(NA_integer_, nrow = length(feature_accessions), ncol = length(sampleorder), dimnames = list(feature_accessions, sampleorder))
            master_sparse_featcounts_index[cell_index] <- analysis_long$RowNumber[placeable]
        }

        ######################
        ## Make feature table
        # Account for the fact that descriptions for the same accession may diverge, due to using jams files annotated by different versions. The function below will choose the most prevalent weighted by sequencing depth.
        # n.b. this is now done for every accession in ONE grouped pass. It used to be a sapply
        # over accessions, each call subsetting the whole long form table, which is quadratic and
        # on a cohort this size would have taken longer than everything else put together.
        ftt <- data.frame(Accession = feature_accessions, stringsAsFactors = FALSE)
        if (!(analysis %in% c("ECNumber", "MetaCyc", "GO"))){
            descstats <- as.data.frame(analysis_long[ , c("Accession", "Description", "NumBases")] %>% group_by(Accession, Description) %>% summarize(Sum = sum(NumBases), .groups = "drop"))
            #n.b. the old code returned EVERY description tied on the maximum, because it had no
            #[1] on the which(). sapply() then degraded to a list, ftt$Description became a list
            #column, and rowData ended up holding the literal text c("desc A", "desc B"). Here
            #exactly one description is always returned, and Description is carried as an
            #explicit third sort key so the tie-break is the alphabetically first one rather
            #than whatever order() happened to leave behind. method = "radix" is named because
            #order()'s default "shell" is NOT stable, and it is what gets chosen as soon as one
            #of the keys is character, which Accession is.
            descstats <- descstats[order(descstats$Accession, -(descstats$Sum), descstats$Description, method = "radix"), , drop = FALSE]
            descstats <- descstats[!(duplicated(descstats$Accession)), , drop = FALSE]
            ftt$Description <- descstats$Description[match(ftt$Accession, descstats$Accession)]
            descstats <- NULL
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

        cts <- cts[order(rowSums(cts), decreasing = TRUE), , drop = FALSE]
        featureorder <- rownames(cts)
        ftt <- ftt[featureorder, , drop = FALSE]
        cts <- cts[ , sampleorder, drop = FALSE]
        cts <- as.matrix(cts)

        #Register the total number of bases sequenced for each sample within that analysis
        TotalBasesSequencedinAnalysis <- colSums(cts)
        TotalBasesSequencedinAnalysis <- t(as.matrix(TotalBasesSequencedinAnalysis))
        rownames(TotalBasesSequencedinAnalysis) <- "NumBases"

        #Order gennumcts and genlencts
        gennumcts <- gennumcts[featureorder, sampleorder, drop = FALSE]
        genlencts <- genlencts[featureorder, sampleorder, drop = FALSE]

        assays <- list()
        assays$BaseCounts <- cts
        assays$GeneCounts <- gennumcts
        assays$GeneLengths <- genlencts
        #Add stratification index as assay if applicable
        if (stratify_functions_by_taxon){
            assays$SparseIndex <- master_sparse_featcounts_index[featureorder, sampleorder, drop = FALSE]
            master_sparse_featcounts_index <- NULL
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

            #Global row names of the stratification matrices. These are the RowNumbers, as
            #characters, which is what retrieve_features_by_taxa indexes them with.
            stratification_rownames <- as.character(seq_len(total_rows))

            base_blocks <- lapply(per_sample_results[has_idx], function(x) { x$base_block } )
            for (k in seq_along(per_sample_results)){
                per_sample_results[[k]][["base_block"]] <- NULL
            }
            master_sparse_taxon_to_space_basecounts_df <- rbind_sparse_blocks(blocks = base_blocks, row_offsets = row_offsets, total_rows = total_rows, rownames_out = stratification_rownames)
            base_blocks <- NULL

            gene_blocks <- lapply(per_sample_results[has_idx], function(x) { x$gene_block } )
            for (k in seq_along(per_sample_results)){
                per_sample_results[[k]][["gene_block"]] <- NULL
            }
            master_sparse_taxon_to_space_genecounts_df <- rbind_sparse_blocks(blocks = gene_blocks, row_offsets = row_offsets, total_rows = total_rows, rownames_out = stratification_rownames)
            gene_blocks <- NULL

            if (!is.null(master_sparse_taxon_to_space_basecounts_df)){
                flog.info(sprintf("make_SummarizedExperiments: %s stratification matrix is %d features-in-samples x %d taxa holding %d non-zero values (%.1f MB).", analysis, nrow(master_sparse_taxon_to_space_basecounts_df), ncol(master_sparse_taxon_to_space_basecounts_df), length(master_sparse_taxon_to_space_basecounts_df@x), (as.numeric(object.size(master_sparse_taxon_to_space_basecounts_df)) / 1e6)))
            }

            #Add features-by-taxon matrix and index into SummarizedExperiment object
            metadata(SEobj)$allfeaturesbytaxa_matrix <- master_sparse_taxon_to_space_basecounts_df
            metadata(SEobj)$allfeaturesbytaxa_GeneCounts_matrix <- master_sparse_taxon_to_space_genecounts_df
            #Add taxonomic level information for downstream aggregation of taxonomy stratification
            metadata(SEobj)$tt <- CGB_tt

            #Clean up
            master_sparse_taxon_to_space_basecounts_df <- NULL
            master_sparse_taxon_to_space_genecounts_df <- NULL
            stratification_rownames <- NULL

        } #End of splitting features by taxonomy

        #Clean up before the next analysis
        per_sample_results <- NULL
        analysis_long <- NULL
        cell_index <- NULL
        acc_pos <- NULL
        samp_pos <- NULL
        placeable <- NULL

        expvec[[e]] <- SEobj
        names(expvec)[e] <- analysis
        e <- e + 1
        #clean memory up
        gc()
    }

    return(expvec)
}