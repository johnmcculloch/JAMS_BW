#' get_genomes_NCBI(organisms = "bacteria", assembly_summary = NULL, outputdir = NULL, assembly_accession = NULL, organism_name = NULL, infraspecific_name = NULL, taxid = NULL, species_taxid = NULL, asm_name = NULL, as_general_expression = FALSE, assembly_level = NULL, assembly_upto = NULL, bioproject = NULL, biosample = NULL, refseq_category_upto = NULL, ntop = NULL, nobs = TRUE, allow_assemblies_from_metagenomes = TRUE, fileformat = "fasta", simulate = TRUE)
#'
#' Downloads Genomes from NCBI GenBank based on matching criteria
#' @export

get_genomes_NCBI <- function(organisms = "bacteria", assembly_summary = NULL, outputdir = NULL, assembly_accession = NULL, organism_name = NULL, infraspecific_name = NULL, taxid = NULL, species_taxid = NULL, asm_name = NULL, as_general_expression = FALSE, assembly_level = NULL, assembly_upto = NULL, bioproject = NULL, biosample = NULL, refseq_category_upto = NULL, ntop = NULL, nobs = TRUE, allow_assemblies_from_metagenomes = TRUE, fileformat = "fasta", simulate = TRUE){

    require(data.table)

    # outputdir is only fundamental when we are actually going to download.
    # When simulating (exploring available genomes), outputdir may remain NULL.
    if (simulate == FALSE){
        if (is.null(outputdir)){
            outputdir <- getwd()
        }
    }

    if (all(c(is.null(organisms), is.null(taxid), is.null(species_taxid), is.null(assembly_accession), is.null(bioproject), is.null(biosample), is.null(organism_name), is.null(infraspecific_name), is.null(refseq_category_upto), is.null(asm_name), is.null(assembly_level), is.null(assembly_upto)))){
        stop("You must provide at least one criterion as input.")
    }

    if (is.null(assembly_summary)){
        flog.info("Downloading assembly_summary.txt from GenBank")
        GBURL <- "https://ftp.ncbi.nlm.nih.gov/genomes/genbank"
        options(timeout = 600)
        if (!is.null(organisms)){
            ASURL <- paste(GBURL, organisms, "assembly_summary.txt", sep = "/")
        } else {
            ASURL <- paste(GBURL, "assembly_summary_genbank.txt", sep = "/")
        }

        # NCBI assembly_summary.txt has TWO header lines:
        #   Line 1: a "# See ftp://..." README comment  <-- must be skipped
        #   Line 2: the actual column names
        # Some rows have more tab-separated fields than the header declares
        # (e.g. embedded tabs in free-text fields like asm_submitter).
        # We pre-read the true column names, then use fill=Inf so fread never
        # stops early, and drop any overflow columns by name afterwards.
        con <- url(ASURL, open = "r")
        readLines(con, n = 1L)                 # discard the README comment line
        header_line <- readLines(con, n = 1L)  # the real column-name line
        close(con)
        true_colnames <- trimws(strsplit(header_line, "\t")[[1]])
        true_colnames <- trimws(gsub("#", "", true_colnames))
        flog.info(paste("Expected number of columns in assembly_summary:", length(true_colnames)))

        assembly_summary <- fread(ASURL, stringsAsFactors = FALSE, sep = "\t", header = TRUE,
                                  quote = "", data.table = FALSE, fill = Inf, skip = 1)
        assembly_summary <- as.data.frame(assembly_summary)

        # Drop any overflow columns (they will be named V39, V40, etc. by fread)
        assembly_summary <- assembly_summary[ , seq_len(length(true_colnames)), drop = FALSE]

        # ── Defensive column renaming ──────────────────────────────────────────
        # NCBI occasionally changes the header line (e.g. adds/removes the '#'
        # prefix or tweaks spacing). Strip '#' and whitespace from every name
        # and warn if the result contains duplicates.
        raw_names <- colnames(assembly_summary)
        clean_names <- trimws(gsub("#", "", raw_names))   # remove # and surrounding whitespace
        clean_names <- make.names(clean_names, unique = TRUE)  # ensure valid + unique R names
        clean_names <- gsub("\\.", "_", clean_names)           # replace dots introduced by make.names

        if (any(duplicated(clean_names))) {
            warning("Duplicate column names detected after cleaning. Check the assembly_summary header.")
        }
        colnames(assembly_summary) <- clean_names

        flog.info(paste("Columns found in assembly_summary:", paste(colnames(assembly_summary), collapse = ", ")))

        # ── Validate that the columns we depend on actually exist ──────────────
        required_cols <- c("taxid", "species_taxid", "infraspecific_name",
                           "asm_name", "organism_name", "assembly_accession",
                           "assembly_level", "refseq_category",
                           "excluded_from_refseq", "ftp_path",
                           "bioproject", "biosample")
        missing_cols <- setdiff(required_cols, colnames(assembly_summary))
        if (length(missing_cols) > 0) {
            stop(paste(
                "The following expected columns are missing from assembly_summary.txt.",
                "NCBI may have changed their file format. Missing columns:",
                paste(missing_cols, collapse = ", ")
            ))
        }

        # ── Type coercions ─────────────────────────────────────────────────────
        assembly_summary$taxid          <- as.character(assembly_summary$taxid)
        assembly_summary$species_taxid  <- as.character(assembly_summary$species_taxid)
        assembly_summary$infraspecific_name <- gsub("strain=", "", assembly_summary$infraspecific_name)
        assembly_summary$asm_name       <- gsub(" ", "_", assembly_summary$asm_name)

        fix_names <- function(asdf, colm){
            asdf[, colm] <- gsub("[^[:alnum:] ]", "_", asdf[, colm])
            asdf[, colm] <- gsub("^root",         "Unclassified", asdf[, colm])
            asdf[, colm] <- gsub("unclassified",  "Unclassified", asdf[, colm])
            asdf[, colm] <- gsub(" ",  "_", asdf[, colm])
            asdf[, colm] <- gsub("__",  "_", asdf[, colm])
            asdf[, colm] <- gsub("___", "_", asdf[, colm])
            asdf[, colm] <- gsub("^_", "",  asdf[, colm])
            asdf[, colm] <- gsub("_sp_", "_Unclassified_", asdf[, colm])
            asdf[, colm] <- gsub("_$", "",  asdf[, colm])
            return(asdf)
        }
        assembly_summary <- fix_names(asdf = assembly_summary, colm = "organism_name")
        assembly_summary <- fix_names(asdf = assembly_summary, colm = "infraspecific_name")
    }

    flog.info(paste("There are", nrow(assembly_summary), "entries for", organisms, "genomes in NCBI GenBank."))
    wanted_assembly_summary <- assembly_summary

    if (!is.null(assembly_accession)){
        .acc <- assembly_accession
        flog.info(paste("Subsetting to entries only with assembly accession(s)", paste0(.acc, collapse=", ")))
        wanted_assembly_summary <- subset(wanted_assembly_summary, (assembly_accession %in% .acc))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the assembly accession criteria."))
    }

    if (!is.null(bioproject)){
        .bp <- bioproject
        flog.info(paste("Subsetting to entries only within BioProject(s)", paste0(.bp, collapse=", ")))
        wanted_assembly_summary <- subset(wanted_assembly_summary, (bioproject %in% .bp))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the BioProject criteria."))
    }

    if (!is.null(biosample)){
        .bs <- biosample
        flog.info(paste("Subsetting to entries only within BioSample(s)", paste0(.bs, collapse=", ")))
        wanted_assembly_summary <- subset(wanted_assembly_summary, (biosample %in% .bs))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the BioSample criteria."))
    }

    # ── Species / taxid filtering ──────────────────────────────────────────────
    if ((!is.null(species_taxid)) | (!is.null(taxid)) | (!is.null(organism_name))){

        if (is.null(species_taxid)){
            flog.info("Subsetting by Organism Name")
            if (as_general_expression != TRUE){
                wanted_species_taxid <- unique(wanted_assembly_summary$species_taxid[which(wanted_assembly_summary[ , "organism_name"] %in% organism_name)])
            } else {
                flog.info("Considering organism_name supplied to be a general expression.")
                wantedrowslist <- unlist(lapply(1:length(organism_name), function(x){ grep(organism_name[x], wanted_assembly_summary[ , "organism_name"]) }))
                wanted_species_taxid <- unique(wanted_assembly_summary$species_taxid[wantedrowslist])
            }
        } else {
            flog.info("Subsetting by Species Tax ID (species_taxid)")
            wanted_species_taxid <- species_taxid
        }

        if (!is.null(taxid)){
            flog.info("Subsetting by Tax ID (taxid)")
            wanted_species_taxid <- unique(append(wanted_species_taxid, unique(wanted_assembly_summary$species_taxid[which(wanted_assembly_summary[ , "taxid"] %in% taxid)])))
        }

        .stid <- wanted_species_taxid
        wanted_assembly_summary <- subset(wanted_assembly_summary, (species_taxid %in% .stid))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Species criteria."))
    }

    # ── Infraspecies / strain / asm_name filtering ─────────────────────────────
    if (!is.null(taxid)){
        flog.info("Filtering by Tax ID")
        .tid <- taxid
        wanted_assembly_summary <- subset(wanted_assembly_summary, (taxid %in% .tid))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Tax ID criteria."))
    } else {
        if (!is.null(infraspecific_name)){
            flog.info("Filtering by strain name")
            .inf <- infraspecific_name
            wanted_assembly_summary <- subset(wanted_assembly_summary, (infraspecific_name %in% .inf))
            flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Strain Name criteria."))
        }

        if (!is.null(asm_name)){
            flog.info("Filtering by asm_name")
            if (as_general_expression != TRUE){
                .asm <- asm_name
                wanted_assembly_summary <- subset(wanted_assembly_summary, (asm_name %in% .asm))
            } else {
                flog.info("Considering asm_name supplied to be a general expression.")
                wantedrowslist <- unlist(lapply(1:length(asm_name), function(x){ grep(asm_name[x], wanted_assembly_summary[ , "asm_name"]) }))
                wanted_assembly_summary <- wanted_assembly_summary[unique(wantedrowslist), ]
            }
            flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the asm name criteria."))
        }
    }

    # ── Assembly level filtering ───────────────────────────────────────────────
    wanted_assembly_level <- c("Contig", "Scaffold", "Chromosome", "Complete Genome")
    if (is.null(assembly_level)){
        if (!is.null(assembly_upto)){
            wanted_assembly_level <- switch(assembly_upto,
                "Contig"      = c("Contig", "Scaffold", "Chromosome", "Complete Genome"),
                "Scaffold"    = c("Scaffold", "Chromosome", "Complete Genome"),
                "Chromosome"  = c("Chromosome", "Complete Genome"),
                "Complete"    = "Complete Genome")
        }
    } else {
        wanted_assembly_level <- assembly_level
    }

    .lvl <- wanted_assembly_level
    wanted_assembly_summary <- subset(wanted_assembly_summary, (assembly_level %in% .lvl))
    flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Assembly Level criteria."))

    if (!is.null(refseq_category_upto)){
        wanted_refseq_category <- switch(refseq_category_upto,
            "reference"      = "reference genome",
            "representative" = c("reference genome", "representative genome"),
            "all"            = c("reference genome", "representative genome", "na"))
        .rsc <- wanted_refseq_category
        flog.info(paste("Subsetting by RefSeq category to include only:", paste0(.rsc, collapse=", ")))
        wanted_assembly_summary <- subset(wanted_assembly_summary, (refseq_category %in% .rsc))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the RefSeq criteria."))
    }

    if (nobs == TRUE){
        if (allow_assemblies_from_metagenomes){
            flog.info("Assemblies derived from metagenomes will be included.")
            permissible_flags <- c("", "na", "derived from metagenome", "from large multi-isolate project")
        } else {
            flog.info("Assemblies derived from metagenomes will be EXCLUDED.")
            permissible_flags <- c("", "na", "from large multi-isolate project")
        }
        flog.info("Filtering out entries flagged as problematic or low quality.")
        noproblem <- which(wanted_assembly_summary$excluded_from_refseq %in% permissible_flags)
        wanted_assembly_summary <- wanted_assembly_summary[noproblem, ]
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries NOT flagged as problematic or low quality."))
    }

    if (nrow(wanted_assembly_summary) > 0){

        filesuffix <- switch(fileformat,
            "genbank" = "genomic.gbff.gz",
            "fasta"   = "genomic.fna.gz",
            "gff"     = "genomic.gff.gz",
            "faa"     = "protein.faa.gz",
            "rna"     = "rna_from_genomic.fna.gz")

        # ── Build the download URL ─────────────────────────────────────────────
        # If a fresh summary was downloaded, ftp_path is present and we build the
        # url from it. If a previously-processed summary was passed back in,
        # ftp_path is gone but a url column already exists, so we reuse it.
        if ("ftp_path" %in% colnames(wanted_assembly_summary)){

            # Drop rows with missing/placeholder ftp_path before URL construction
            valid_ftp <- !is.na(wanted_assembly_summary$ftp_path) &
                         wanted_assembly_summary$ftp_path != "na" &
                         wanted_assembly_summary$ftp_path != ""
            n_dropped <- sum(!valid_ftp)
            if (n_dropped > 0){
                flog.info(paste("Dropping", n_dropped, "entries with missing ftp_path."))
                wanted_assembly_summary <- wanted_assembly_summary[valid_ftp, ]
            }

            sufreps <- sapply(seq_len(nrow(wanted_assembly_summary)), function(x){
                parts <- unlist(strsplit(wanted_assembly_summary$ftp_path[x], split = "/"))
                parts[length(parts)]
            })

            wanted_assembly_summary$url    <- paste(wanted_assembly_summary$ftp_path, paste(sufreps, filesuffix, sep = "_"), sep = "/")
            wanted_assembly_summary$ftp_path <- NULL

        } else if ("url" %in% colnames(wanted_assembly_summary)){
            flog.info("Input summary has no ftp_path; reusing existing url column.")
        } else {
            stop("assembly_summary has neither an ftp_path nor a url column; cannot construct download URLs.")
        }

        # ── Always (re)compute destfn from the CURRENT outputdir ────────────────
        # This makes the outputdir of the current call authoritative, even if a
        # summary carrying a stale destfn from an earlier call is passed in.
        # destfn only matters when actually downloading; when simulating with a
        # NULL outputdir we skip it (there is nowhere to write to yet).
        if (!is.null(outputdir)){
            wanted_assembly_summary$destfn <- file.path(
                outputdir,
                paste(wanted_assembly_summary$organism_name,
                      wanted_assembly_summary$infraspecific_name,
                      wanted_assembly_summary$assembly_accession,
                      filesuffix, sep = "_"))
        } else {
            wanted_assembly_summary$destfn <- NULL
        }

        # ── ntop guard: don't request more rows than exist ─────────────────────
        if (!is.null(ntop)){
            ntop_actual <- min(ntop, nrow(wanted_assembly_summary))
            if (ntop_actual < ntop){
                flog.info(paste("Requested ntop =", ntop, "but only", ntop_actual, "entries available. Keeping all."))
            } else {
                flog.info(paste("Keeping only the top", ntop_actual, "genomes."))
            }
            wanted_assembly_summary <- wanted_assembly_summary[seq_len(ntop_actual), ]
        }

        if (simulate == FALSE){
            if (!dir.exists(outputdir)){
                flog.info(paste("Creating output directory:", outputdir))
                dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
            }
            flog.info("Downloading genomes")
            dnl <- sapply(seq_len(nrow(wanted_assembly_summary)), function(x){
                tryCatch(
                    download.file(wanted_assembly_summary$url[x], destfile = wanted_assembly_summary$destfn[x]),
                    error = function(e) flog.warn(paste("Unable to download", wanted_assembly_summary$url[x]))
                )
            })
        }

    } else {
        flog.info("No genomes fit the criteria imposed.")
        wanted_assembly_summary <- NULL
    }

    return(wanted_assembly_summary)
}