#' build_sample_stratification_block(Sample = NULL, list.data = NULL, analysis = NULL, valid_taxonomic_space = NULL, appropriate_featurdata_tax_colm = "LKT", stratify_functions_by_taxon = TRUE, ULA_basecount_threshold = 10000)
#'
#' PER-SAMPLE WORKER for the functional phase of make_SummarizedExperiments (see). Given ONE
#' sample, it returns everything that sample contributes to a functional SummarizedExperiment
#' for ONE analysis, and nothing else. It is a pure function of its arguments: it reads
#' list.data but never writes to it, it emits no log lines (messages are returned so the master
#' can log them in order), and it never calls gc(). This makes it safe to run inside a forked
#' parallel::mclapply worker.
#'
#' WHY THIS EXISTS. The functional loop used to do all of this inline, on the master, one sample
#' at a time, accumulating into growing objects with rbind() and merge_sparse_matrix(). With
#' ~1500 samples that accumulation is quadratic and single threaded, and it dominated the run
#' (ETA ~26 h for 1409 samples). Splitting the per-sample work out into this function allows the
#' master to (a) run it on many cores at once and (b) combine all the results in a single pass at
#' the end (see rbind_sparse_blocks).
#'
#' WHAT IS RETURNED. A list with:
#'   Sample       - the sample name
#'   has_data     - TRUE if this sample has any feature for this analysis
#'   failed       - TRUE if the worker trapped an error (the master decides what to do)
#'   error        - the error message, if any
#'   messages     - character vector of things the master should flog.warn
#'   accessions   - character vector, ONE ENTRY PER ROW the sample contributes, in row order
#'   descriptions - character vector, aligned to accessions
#'   numbases     - numeric, aligned to accessions (unstratified BaseCounts)
#'   numgenes     - numeric, aligned to accessions (unstratified GeneCounts)
#'   lengenes     - numeric, aligned to accessions (unstratified GeneLengths)
#'   base_block   - dgCMatrix, nrow = length(accessions), cols named by taxon. NULL if not stratifying.
#'   gene_block   - dgCMatrix, same dims and SAME COLUMN ORDER as base_block. NULL if not stratifying.
#'
#' The blocks deliberately carry NO rownames. Global row numbering can only be decided by the
#' master, once it knows how many rows every sample contributed, so it is applied at assembly time.
#'
#' @param Sample Character. Name of the sample to process.
#' @param list.data The JAMSbeta list.data object. Passed by reference; never modified here.
#' @param analysis Character. The functional analysis being built, e.g. "ECNumber", "Product".
#' @param valid_taxonomic_space Character. Taxonomic space of the functional abundance table.
#' @param appropriate_featurdata_tax_colm Character. Column of featuredata holding the taxon.
#' @param stratify_functions_by_taxon Logical. If FALSE, no sparse blocks are built.
#' @param ULA_basecount_threshold Numeric. Taxa whose total bases in this sample fall below this
#'   are collapsed into a single "Ultra_low_abundance_LKTs" column. Default 10000 (10 kbp), which
#'   is the long-standing JAMS value; see the rationale in make_SummarizedExperiments.
#' @export

build_sample_stratification_block <- function(Sample = NULL, list.data = NULL, analysis = NULL, valid_taxonomic_space = NULL, appropriate_featurdata_tax_colm = "LKT", stratify_functions_by_taxon = TRUE, ULA_basecount_threshold = 10000){

    SN <- Sample
    msgs <- character(0)

    empty_result <- function(messages = character(0)){
        list(Sample = SN, has_data = FALSE, failed = FALSE, error = NULL, messages = messages,
            accessions = character(0), descriptions = character(0), numbases = numeric(0),
            numgenes = numeric(0), lengenes = numeric(0), base_block = NULL, gene_block = NULL)
    }

    curr_FD <- list.data[[paste(SN, "abundances", sep = "_")]]$functional[[valid_taxonomic_space]]
    if (is.null(curr_FD)){
        return(empty_result(msgs))
    }

    wantedrows <- which(curr_FD$Analysis == analysis)
    if (length(wantedrows) < 1){
        return(empty_result(msgs))
    }
    curr_FD <- curr_FD[wantedrows, , drop = FALSE]

    accessions <- as.character(curr_FD$Accession)
    nfeat <- length(accessions)

    if (!is.null(curr_FD$Description)){
        descriptions <- as.character(curr_FD$Description)
    } else {
        descriptions <- rep("none", nfeat)
    }
    descriptions[is.na(descriptions)] <- "none"

    numbases <- suppressWarnings(as.numeric(curr_FD$NumBases))
    if (length(numbases) != nfeat){
        numbases <- rep(0, nfeat)
    }
    numbases[is.na(numbases)] <- 0

    #####################################################################################
    ## Featuredata side. This supplies gene tallies (stratified and unstratified) and gene
    ## lengths. It is read ONCE here, rather than being sliced repeatedly by dplyr.
    #####################################################################################
    curr_featdf_full <- list.data[[paste(SN, "featuredata", sep = "_")]]

    acc_vec <- NULL
    tax_vec <- NULL
    len_vec <- NULL

    if (is.null(curr_featdf_full)){
        msgs <- c(msgs, paste0("make_SummarizedExperiments: featuredata is missing for sample ", SN, ". Its BaseCounts for analysis ", analysis, " are retained, but GeneCounts and GeneLengths will be zero and it will contribute no taxon-stratified gene counts."))
    } else {
        nrows_fd <- nrow(curr_featdf_full)

        #Defensive: the column named by `analysis` (e.g. resfinder) may be ABSENT from a
        #sample's featuredata if that blast analysis never ran at JAMSalpha time (as opposed
        #to running and finding nothing, which yields a column full of "none"). Synthesize it
        #as "none" so the sample stays in the analysis, mirroring how ECNumber is handled in
        #make_featuredata_from_bedfile. n.b. we build a standalone vector rather than writing a
        #column into curr_featdf_full, because the latter would deep-copy the whole featuredata.
        if (!(analysis %in% colnames(curr_featdf_full))){
            msgs <- c(msgs, paste0("make_SummarizedExperiments: analysis column '", analysis, "' is absent from featuredata for sample ", SN, " (this analysis likely did not run for this sample at JAMSalpha time). Synthesizing it as 'none' so the sample is retained."))
            acc_vec <- rep("none", nrows_fd)
        } else {
            acc_vec <- as.character(curr_featdf_full[[analysis]])
        }

        if (!("LengthDNA" %in% colnames(curr_featdf_full))){
            msgs <- c(msgs, paste0("make_SummarizedExperiments: featuredata for sample ", SN, " has no LengthDNA column. GeneLengths for this sample will be zero."))
            len_vec <- rep(0, nrows_fd)
        } else {
            len_vec <- suppressWarnings(as.numeric(curr_featdf_full[["LengthDNA"]]))
        }

        if (!(appropriate_featurdata_tax_colm %in% colnames(curr_featdf_full))){
            msgs <- c(msgs, paste0("make_SummarizedExperiments: featuredata for sample ", SN, " has no ", appropriate_featurdata_tax_colm, " column, so this sample cannot contribute taxon-stratified gene counts."))
            tax_vec <- rep(NA_character_, nrows_fd)
        } else {
            tax_vec <- as.character(curr_featdf_full[[appropriate_featurdata_tax_colm]])
        }

        #Fix "none" to analysis_none
        acc_vec[which(acc_vec == "none")] <- paste(analysis, "none", sep = "_")
        len_vec[is.na(len_vec)] <- 0

        #Separate rows if not Product or ECNumber, as there may be more than a single accession
        #for a single gene separated by a "|". Kept as the original tidyr call, deliberately, so
        #that whatever this currently does for Pfam/Interpro/resfinder/PRINTS/GO keeps doing it.
        if (!(analysis %in% c("Product", "ECNumber"))){
            sep_df <- data.frame(Accession = acc_vec, LengthDNA = len_vec, Taxon = tax_vec, stringsAsFactors = FALSE)
            sep_df <- as.data.frame(sep_df %>% tidyr::separate_rows(all_of("Accession"), sep = fixed("\\|")))
            acc_vec <- as.character(sep_df$Accession)
            len_vec <- as.numeric(sep_df$LengthDNA)
            tax_vec <- as.character(sep_df$Taxon)
            sep_df <- NULL
        }
    }

    #Row position in this sample's block for every featuredata row, NA when the accession is not
    #among the accessions this analysis reports for this sample.
    if (!is.null(acc_vec)){
        feat_rowpos <- match(acc_vec, accessions)
    } else {
        feat_rowpos <- integer(0)
    }
    usable_feat <- which(!is.na(feat_rowpos))

    #####################################################################################
    ## Unstratified gene tallies and gene lengths, aligned to the accession order above.
    #####################################################################################
    numgenes <- numeric(nfeat)
    lengenes <- numeric(nfeat)
    if (length(usable_feat) > 0){
        numgenes <- as.numeric(tabulate(feat_rowpos[usable_feat], nbins = nfeat))
        #rowsum() labels its result rows with as.character(group). group is integer here (it
        #comes from match()), so the labels are plain digits and round-trip exactly; this would
        #NOT be safe with a double group, where as.character(1e+05) is "1e+05".
        len_groups <- as.integer(feat_rowpos[usable_feat])
        len_by_row <- rowsum(len_vec[usable_feat], group = len_groups, reorder = FALSE)
        lengenes[as.integer(rownames(len_by_row))] <- as.numeric(len_by_row[ , 1])
        len_by_row <- NULL
        len_groups <- NULL
    }

    if (!stratify_functions_by_taxon){
        return(list(Sample = SN, has_data = TRUE, failed = FALSE, error = NULL, messages = msgs,
            accessions = accessions, descriptions = descriptions, numbases = numbases,
            numgenes = numgenes, lengenes = lengenes, base_block = NULL, gene_block = NULL))
    }

    #####################################################################################
    ## Taxon-stratified BaseCounts block.
    ##
    ## The abundance table is Analysis, Accession, Description, NumBases, then one column per
    ## taxon. We walk those taxon columns one at a time and keep only the non-zero entries,
    ## so we never materialise the dense nfeat x ntaxa matrix that as.matrix() used to build.
    #####################################################################################
    if (ncol(curr_FD) >= 5){
        Taxoncols <- colnames(curr_FD)[5:ncol(curr_FD)]
    } else {
        Taxoncols <- character(0)
    }

    if (length(Taxoncols) < 1){
        msgs <- c(msgs, paste0("make_SummarizedExperiments: the ", valid_taxonomic_space, " functional abundance table for sample ", SN, " carries no taxon columns, so this sample contributes no taxon-stratified data for analysis ", analysis, "."))
        return(list(Sample = SN, has_data = TRUE, failed = FALSE, error = NULL, messages = msgs,
            accessions = accessions, descriptions = descriptions, numbases = numbases,
            numgenes = numgenes, lengenes = lengenes, base_block = NULL, gene_block = NULL))
    }

    idx_list <- vector("list", length(Taxoncols))
    val_list <- vector("list", length(Taxoncols))
    taxon_colsums <- numeric(length(Taxoncols))

    for (j in seq_along(Taxoncols)){
        #n.b. the column is fetched POSITIONALLY, not by name. Taxoncols can legitimately contain
        #repeated names (see the collapsing block below), and curr_FD[[ "someName" ]] returns the
        #FIRST column with that name, so a by-name fetch would read the first bin's column twice
        #and never read the second one at all. Taxon column j is physical column j + 4, because
        #the table is Analysis, Accession, Description, NumBases, then the taxa (see
        #harvest_functions, which pins that layout).
        v <- curr_FD[[j + 4L]]
        #as.numeric unconditionally: these columns are often INTEGER, and sum() on an integer
        #vector overflows to NA past 2^31, which one dominant taxon in a deeply sequenced
        #metagenome will reach. Factors go via as.character first, because as.numeric() on a
        #factor returns its integer codes rather than its values.
        if (is.factor(v)){
            v <- as.character(v)
        }
        v <- suppressWarnings(as.numeric(v))
        v[is.na(v)] <- 0
        nz <- which(v != 0)
        idx_list[[j]] <- nz
        val_list[[j]] <- v[nz]
        taxon_colsums[j] <- sum(v)
    }

    #########################################################################################
    ## Collapse taxon columns that share a name.
    ##
    ## This is NOT hypothetical. make_SummarizedExperiments renames the functional abundance
    ## table's ConsolidatedGenomeBin columns to their contextualized LKT, and two different
    ## MetaBAT2 bins in the SAME sample can end up with the same LKT name (rare, around 0.005%
    ## of Sample/LKT pairs, which over a cohort of a few hundred samples is a near certainty).
    ## The taxonomic SummarizedExperiment already sums such pairs with
    ## group_by(Sample, LKT) %>% summarize(NumBases = sum(NumBases)); we do exactly the same
    ## here so the stratified matrix stays consistent with it, and so no sparse matrix is ever
    ## handed two columns of the same name.
    ##
    ## Summing happens BEFORE the ultra-low-abundance threshold is applied, so the threshold
    ## sees the taxon's real total rather than one bin's share of it.
    #########################################################################################
    if (anyDuplicated(Taxoncols) > 0){
        unique_tax <- unique(Taxoncols)
        tax_group <- match(Taxoncols, unique_tax)
        merged_idx <- vector("list", length(unique_tax))
        merged_val <- vector("list", length(unique_tax))
        merged_colsums <- numeric(length(unique_tax))

        for (g in seq_along(unique_tax)){
            members <- which(tax_group == g)
            if (length(members) == 1){
                merged_idx[[g]] <- idx_list[[members]]
                merged_val[[g]] <- val_list[[members]]
            } else {
                acc <- numeric(nfeat)
                for (mm in members){
                    if (length(idx_list[[mm]]) > 0){
                        acc[idx_list[[mm]]] <- acc[idx_list[[mm]]] + val_list[[mm]]
                    }
                }
                acc_nz <- which(acc != 0)
                merged_idx[[g]] <- acc_nz
                merged_val[[g]] <- acc[acc_nz]
                acc <- NULL
            }
            merged_colsums[g] <- sum(taxon_colsums[members])
        }

        msgs <- c(msgs, paste0("make_SummarizedExperiments: sample ", SN, " has ", (length(Taxoncols) - length(unique_tax)), " taxon column(s) in its ", valid_taxonomic_space, " functional abundance table that share an LKT name with another column (two genome bins contextualized to the same taxon). Their bases and genes have been summed into one column, which is what the taxonomic SummarizedExperiment does with the same pairs."))

        Taxoncols <- unique_tax
        idx_list <- merged_idx
        val_list <- merged_val
        taxon_colsums <- merged_colsums
        merged_idx <- NULL
        merged_val <- NULL
    }

    #In order to conserve RAM and the size of the taxonomic stratification sparse matrix, do not
    #include any taxon whose total sum of bases across the entire metagenomic proteome is smaller
    #than 10000 bp (10 kbp). 10 kbp / ~3 Mbp is a completeness of 0.3%, with a sequencing depth of
    #1X if the entire sample consisted of a single species (which it doesn't). So this filtration
    #is very conservative and will shrink memory usage considerably.
    lowcols <- which(taxon_colsums < ULA_basecount_threshold)
    keptcols <- setdiff(seq_along(Taxoncols), lowcols)

    taxon_names <- Taxoncols[keptcols]
    final_idx <- idx_list[keptcols]
    final_val <- val_list[keptcols]

    if (length(lowcols) > 0){
        #Sum all low abundance values into a single column, appended last (exactly where the old
        #cbind() put it, so column order is unchanged).
        ULA_sums <- numeric(nfeat)
        for (j in lowcols){
            if (length(idx_list[[j]]) > 0){
                ULA_sums[idx_list[[j]]] <- ULA_sums[idx_list[[j]]] + val_list[[j]]
            }
        }
        ULA_nz <- which(ULA_sums != 0)
        taxon_names <- c(taxon_names, "Ultra_low_abundance_LKTs")
        final_idx <- c(final_idx, list(ULA_nz))
        final_val <- c(final_val, list(ULA_sums[ULA_nz]))
        ULA_sums <- NULL
    }

    idx_list <- NULL
    val_list <- NULL

    base_i <- unlist(final_idx, use.names = FALSE)
    base_x <- unlist(final_val, use.names = FALSE)
    if (is.null(base_i)){
        base_i <- integer(0)
    }
    if (is.null(base_x)){
        base_x <- numeric(0)
    }

    #Assembled directly as compressed-column storage. which() returns ascending indices, so the
    #row indices are already sorted within each column, which is the one invariant CSC needs.
    #Going through the dgCMatrix constructor rather than Matrix::sparseMatrix() keeps the
    #all-zero case (a sample where every taxon column is empty) unambiguous. n.b. the i slot is
    #0 based, hence the -1L.
    base_block <- new("dgCMatrix",
        i = as.integer(base_i) - 1L,
        p = as.integer(c(0, cumsum(lengths(final_idx)))),
        x = as.numeric(base_x),
        Dim = c(as.integer(nfeat), as.integer(length(taxon_names))),
        Dimnames = list(NULL, taxon_names))

    final_idx <- NULL
    final_val <- NULL
    base_i <- NULL
    base_x <- NULL

    #####################################################################################
    ## Taxon-stratified GeneCounts block. Must have EXACTLY the same dimensions, column
    ## order and ULA collapsing as base_block, because downstream code indexes both with
    ## the same row numbers and expects matched columns.
    ##
    ## Built directly as a sparse matrix from the (accession, taxon) pairs in featuredata.
    ## This replaces a group_by %>% summarise %>% pivot_wider that produced a DENSE
    ## accessions x taxa data frame per sample before being made sparse again.
    #####################################################################################
    #Map every original taxon column to its destination column in base_block: itself if kept,
    #the Ultra_low_abundance_LKTs column if it was collapsed.
    tax_target <- rep(NA_integer_, length(Taxoncols))
    if (length(keptcols) > 0){
        tax_target[keptcols] <- seq_along(keptcols)
    }
    if (length(lowcols) > 0){
        tax_target[lowcols] <- length(taxon_names)
    }

    if (length(usable_feat) > 0){
        gene_colpos <- tax_target[match(tax_vec[usable_feat], Taxoncols)]
        gene_rowpos <- feat_rowpos[usable_feat]
        gene_ok <- which(!is.na(gene_colpos))
    } else {
        gene_colpos <- integer(0)
        gene_rowpos <- integer(0)
        gene_ok <- integer(0)
    }

    #n.b. use.last.ij = FALSE is the sparseMatrix() default and is what we want: repeated
    #(accession, taxon) pairs are SUMMED, which is precisely the gene tally. It is passed
    #explicitly so the intent is not mistaken for an oversight.
    if (length(gene_ok) > 0){
        gene_block <- Matrix::sparseMatrix(i = as.integer(gene_rowpos[gene_ok]),
            j = as.integer(gene_colpos[gene_ok]),
            x = rep(1, length(gene_ok)),
            dims = c(as.integer(nfeat), as.integer(length(taxon_names))),
            dimnames = list(NULL, taxon_names),
            index1 = TRUE, use.last.ij = FALSE)
    } else {
        #No gene at all could be placed. Hand back an explicitly empty matrix of the right shape
        #rather than relying on sparseMatrix() to cope with zero-length inputs.
        gene_block <- new("dgCMatrix",
            i = integer(0),
            p = integer(length(taxon_names) + 1),
            x = numeric(0),
            Dim = c(as.integer(nfeat), as.integer(length(taxon_names))),
            Dimnames = list(NULL, taxon_names))
    }

    return(list(Sample = SN, has_data = TRUE, failed = FALSE, error = NULL, messages = msgs,
        accessions = accessions, descriptions = descriptions, numbases = numbases,
        numgenes = numgenes, lengenes = lengenes,
        base_block = base_block, gene_block = gene_block))
}
