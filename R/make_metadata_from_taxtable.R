#' make_metadata_from_taxtable(taxlvlinterest = NULL, taxainterest = NULL, nobs = TRUE,  metadatafn = NULL)
#'
#' Generates JAMS-style metadata in xlsx format for comparing bacterial isolates which have been deposited in NCBI GenBank
#' @export

make_metadata_from_taxtable <- function(taxlvlinterest = NULL, taxainterest = NULL, nobs = TRUE, assembly_upto = NULL, metadatafn = NULL){

    #Define useful functions
    is.redundant <- function(vec){
        propunique <- length(unique(vec)) / length(vec)
        if (propunique < 1){
            redundant <- TRUE
        } else {
            redundant <- FALSE
        }

        return(redundant)
    }

    is.useful <- function(vec){
        numunique <- length(unique(vec))
        if (numunique < 2){
            useful = FALSE
        } else {
            useful = TRUE
        }

        return(useful)
    }

    tt <- get_genomes_NCBI(assembly_upto = assembly_upto, nobs = nobs)
    data(JAMStaxtable)
    cinterest <- c("assembly_accession", "species_taxid", "organism_name", "infraspecific_name", "seq_rel_date", "refseq_category", "excluded_from_refseq", "assembly_level", "paired_asm_comp", "relation_to_type_material", "bioproject", "biosample")
    if ((taxlvlinterest %in% colnames(JAMStaxtable)) == FALSE){
        stop("The taxon level of interest could not be found. Try something such as Genus.")
	  } else {
        if (all(taxainterest %in% JAMStaxtable[, taxlvlinterest]) == FALSE){
            stop("Unable to find your taxon of interest. Try something such as g__Blautia.")
        }
    }

    #Build Metadata
    flog.info("Found taxa within taxlevel you are interested in. Making Metadata.")
    taxidsinterest <- JAMStaxtable[ , "Taxid"][which(JAMStaxtable[ , taxlvlinterest] %in% taxainterest)]
    ttinterest <- subset(tt, taxid %in% taxidsinterest)
    phenotable <- ttinterest[ , cinterest]

    #Clean up metadata
    phenotable$Year <- sapply(1:nrow(phenotable), function (x) { as.character(unlist(strsplit(phenotable$seq_rel_date[x], split = "\\/"))[1]) } )
    phenotable$Year[which(phenotable$Year %in% c(NA, ""))] <- "1982"
    phenotable$Year <- as.numeric(phenotable$Year)
    phenotable$seq_rel_date <- as.Date(phenotable$seq_rel_date, optional = TRUE)
    phenotable$Numdate <- as.numeric(phenotable$seq_rel_date)
    phenotable$Numdate[which(is.na(phenotable$Numdate))] <- 4383 #If date is absent consider GenBanks creation date 1982-01-01
    phenotable$NCBI_LKT <- sapply(1:nrow(phenotable), function (x) { paste0("LKT__s__", paste(unlist(strsplit(phenotable$organism_name[x], split = "_"))[1:2], collapse="_")) })
    phenotable$excluded_from_refseq[which(phenotable$excluded_from_refseq == "")] <- "OK"
    phenotable$refseq_category[which(phenotable$refseq_category == "na")] <- "none"
    phenotable$infraspecific_name[which(phenotable$infraspecific_name == "")] <- "none"
    phenotable$infraspecific_name[which(is.na(phenotable$infraspecific_name))] <- "none"
    phenotable$paired_asm_comp[which(phenotable$paired_asm_comp == "")] <- "none"
    phenotable$relation_to_type_material[which(phenotable$relation_to_type_material == "")] <- "none"
    phenotable$Sample <- gsub("_", "", phenotable$assembly_accession)
    phenotable$Sample <- sapply(1:nrow(phenotable), function (x) { as.character(unlist(strsplit(phenotable$Sample[x], split = "\\."))[1]) } )
    phenotable$species_taxid <- paste("Species_taxid", phenotable$species_taxid, sep = "_")

    colsofinterest <- c("Sample", "assembly_accession", "assembly_level", "NCBI_LKT", "infraspecific_name", "species_taxid", "Year", "seq_rel_date", "Numdate", "refseq_category", "excluded_from_refseq", "paired_asm_comp", "relation_to_type_material", "bioproject", "biosample")
    colsofinterestcats <- c("Sample", "ID", "discrete", "discrete", "discrete", "discrete", "continuous", "date", "continuous", "discrete", "discrete", "discrete", "discrete", "ID", "ID")
    phenotable <- phenotable[, colsofinterest]

    #vet columns in phenotable
    #Are prefixes redundant
    redundancy <- sapply(1:length(c("Sample", "assembly_accession")), function (x) { is.redundant(phenotable[, c("Sample", "assembly_accession")[x]]) })
    names(redundancy) <- c("Sample", "assembly_accession")
    if (!all(!(redundancy))){
        #Warn that sample prexifes are not OK
        notOK <- names(redundancy)[which(redundancy == TRUE)]
        print(paste("WARNING! Column(s)", paste0(notOK, collapse = ", "), "are not unique. Curate phenotable manually before using."))
     }
    #Are other columns useful
    nonsampcols <- colnames(phenotable)[!(colnames(phenotable) %in% c("NCBI_LKT", "Sample", "assembly_accession"))]
    usefulness <- sapply(1:length(nonsampcols), function (x) { is.useful(phenotable[, nonsampcols[x]]) })
    names(usefulness) <- nonsampcols
    if (any(!usefulness)){
        #Warn that sample prexifes are not OK
        notOK <- names(usefulness)[which(usefulness == FALSE)]
        print(paste("Column(s)", paste0(notOK, collapse = ", "), "are not informative. Eliminating these variables from the metadata"))
        phenotable[ , notOK] <- NULL
     }

    phenotablecolspresent <- colnames(phenotable)[colnames(phenotable) %in% colsofinterest]
    phenolabelspresent <- colsofinterestcats[colsofinterest %in% colnames(phenotable)]

    phenolabels <- data.frame(Var_label = phenotablecolspresent, Var_type = phenolabelspresent, stringsAsFactors = FALSE)

    metadata <- list()
    metadata[[1]] <- phenotable
    metadata[[2]] <- phenolabels
    names(metadata) <- c("phenotable", "phenolabels")

    if (!(is.null(metadatafn))){
		    write.xlsx(metadata, file = metadatafn, asTable = TRUE, rowNames = FALSE, colNames = TRUE, borders = "surrounding", colWidths = "auto")
    }

    return(metadata)
}
