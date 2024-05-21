#' get_genomes_NCBI(organisms="bacteria", assembly_summary=NULL, outputdir=NULL, assembly_accession=NULL, organism_name=NULL, infraspecific_name=NULL, taxid=NULL, species_taxid=NULL, asm_name=NULL, as_general_expression=FALSE, assembly_level=NULL, assembly_upto=NULL, bioproject=NULL, biosample=NULL, refseq_category_upto=NULL, ntop=NULL, nobs=TRUE, fileformat="fasta", simulate=TRUE)
#'
#' Downloads Genomes from NCBI GenBank based on matching criteria
#' @export

get_genomes_NCBI <- function(organisms = "bacteria", assembly_summary = NULL, outputdir = NULL, assembly_accession = NULL, organism_name = NULL, infraspecific_name = NULL, taxid = NULL, species_taxid = NULL, asm_name = NULL, as_general_expression = FALSE, assembly_level = NULL, assembly_upto = NULL, bioproject = NULL, biosample = NULL, refseq_category_upto = NULL, ntop = NULL, nobs = TRUE, fileformat = "fasta", simulate = TRUE){

    require(data.table)
    if (is.null(outputdir)){
        outputdir <- getwd()
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

        assembly_summary <- fread(ASURL, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "", data.table = FALSE)

        assembly_summary <- as.data.frame(assembly_summary)
        #Fix column names
        colnames(assembly_summary) <- gsub(" ", "", gsub("#", "", colnames(assembly_summary)))
        #Correct table to make consistent with JAMS taxtable
        assembly_summary$taxid <- as.character(assembly_summary$taxid)
        assembly_summary$species_taxid <- as.character(assembly_summary$species_taxid)
        assembly_summary$infraspecific_name <- gsub("strain=", "", assembly_summary$infraspecific_name)
        assembly_summary$asm_name <- gsub(" ", "_", assembly_summary$asm_name)
        fix_names<-function(asdf, colm){
            asdf[, colm] <- gsub("[^[:alnum:] ]", "_", asdf[, colm])
            asdf[, colm] <- gsub("^root", "Unclassified", asdf[, colm])
            asdf[, colm] <- gsub("unclassified", "Unclassified", asdf[, colm])
            asdf[, colm] <- gsub(" ", "_", asdf[, colm])
            asdf[, colm] <- gsub("__", "_", asdf[, colm])
            asdf[, colm] <- gsub("___", "_", asdf[, colm])
            asdf[, colm] <- gsub("^_", "", asdf[, colm])
            asdf[, colm] <- gsub("_sp_", "_Unclassified_", asdf[, colm])
            asdf[, colm] <- gsub("_$", "", asdf[, colm])

            return(asdf)
        }
        assembly_summary <- fix_names(asdf=assembly_summary, colm="organism_name")
        assembly_summary <- fix_names(asdf=assembly_summary, colm="infraspecific_name")
    }

    flog.info(paste("There are", nrow(assembly_summary), "entries for", organisms, "genomes in NCBI GenBank."))
    wanted_assembly_summary <- assembly_summary

    if (!is.null(assembly_accession)){
        flog.info(paste("Subsetting to entries only with assembly accession(s)", paste0(assembly_accession, collapse=", ")))
        wanted_assembly_accession <- assembly_accession
        wanted_assembly_summary <- subset(wanted_assembly_summary, (assembly_accession %in% wanted_assembly_accession))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the assembly accession criteria."))
    }

    if (!is.null(bioproject)){
        flog.info(paste("Subsetting to entries only within BioProject(s)", paste0(bioproject, collapse=", ")))
        wanted_bioproject<-bioproject
        wanted_assembly_summary<-subset(wanted_assembly_summary, (bioproject %in% wanted_bioproject))
        print(paste("There are", nrow(wanted_assembly_summary), "entries matching the BioProject criteria."))
    }

    if(!is.null(biosample)){
        flog.info(paste("Subsetting to entries only within BioSample(s)", paste0(biosample, collapse=", ")))
        wanted_biosample<-biosample
        wanted_assembly_summary<-subset(wanted_assembly_summary, (biosample %in% wanted_biosample))
        print(paste("There are", nrow(wanted_assembly_summary), "entries matching the BioSample criteria."))
    }

    #Define species input criterion
    if((!is.null(species_taxid)) | (!is.null(taxid)) | (!is.null(organism_name))){

        if(is.null(species_taxid)){
            #Get wanted species taxid from organism_name
            flog.info("Subsetting by Organism Name")
            if(as_general_expression != TRUE){
                 wanted_species_taxid <- unique(wanted_assembly_summary$species_taxid[which(wanted_assembly_summary[ , "organism_name"] %in% organism_name)])
            } else {
                flog.info("Considering organism_name supplied to be a general expression. Will select all Organism Names matching the expression(s) supplied.")
                wantedrowslist <- unlist(lapply(1:length(organism_name), function (x) { grep(organism_name[x], wanted_assembly_summary[ , "organism_name"]) }))
                wanted_species_taxid <- unique(wanted_assembly_summary$species_taxid[wantedrowslist])
            }
        } else {
            flog.info("Subsetting by Species Tax ID (species_taxid)")
            wanted_species_taxid <- species_taxid
        }

        if(!is.null(taxid)){
            print("Subsetting by Tax ID (taxid)")
            wanted_species_taxid<-unique(append(wanted_species_taxid, unique(wanted_assembly_summary$species_taxid[which(wanted_assembly_summary[ , "taxid"] %in% taxid)])))
        }

        wanted_assembly_summary <- subset(wanted_assembly_summary, (species_taxid %in% wanted_species_taxid))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Species criteria."))
    }

    #Subset by infraspecies criteria
    if(!is.null(taxid)){
        flog.info("Filtering by Tax ID")
        wantedtaxid<-taxid
        wanted_assembly_summary<-subset(wanted_assembly_summary, (taxid %in% wantedtaxid))
        print(paste("There are", nrow(wanted_assembly_summary), "entries matching the Tax ID criteria."))
    } else {
        if (!is.null(infraspecific_name)){
            flog.info("Filtering by strain name")
            strainnames<-infraspecific_name
            wanted_assembly_summary<-subset(wanted_assembly_summary, (infraspecific_name %in% strainnames))
            flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Strain Name criteria."))
        }

        if (!is.null(asm_name)){
            flog.info("Filtering by asm_name")
            asmnames <- asm_name
            if (as_general_expression != TRUE){
                wanted_assembly_summary <- subset(wanted_assembly_summary, (asm_name %in% asmnames))
            } else {
                flog.info("Considering asm_name supplied to be a general expression. Will select all asm names matching the expression(s) supplied.")
                wantedrowslist <- unlist(lapply(1:length(asmnames), function (x) { grep(asmnames[x], wanted_assembly_summary[ , "asm_name"]) }))
                wanted_assembly_summary <- wanted_assembly_summary[unique(wantedrowslist) , ]
            }
            flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the asm name criteria."))
        }
    }

    ##Subset by quality criteria
    #By assembly level
    wanted_assembly_level <- c("Contig", "Scaffold", "Chromosome", "Complete Genome")
    if (is.null(assembly_level)){
        if(!is.null(assembly_upto)){
            wanted_assembly_level <- switch(assembly_upto, "Contig"=c("Contig", "Scaffold", "Chromosome", "Complete Genome"), "Scaffold"=c("Scaffold", "Chromosome", "Complete Genome"), "Chromosome"=c("Chromosome", "Complete genome"), "Complete"="Complete Genome")
        }
    } else {
        wanted_assembly_level <- assembly_level
    }

    wanted_assembly_summary <- subset(wanted_assembly_summary, (assembly_level %in% wanted_assembly_level))
    flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the Assembly Level criteria."))

    if (!is.null(refseq_category_upto)){
        wanted_refseq_category <- switch(refseq_category_upto, "reference" = "reference genome", "representative" = c("reference genome", "representative genome"), "all" = c("reference genome", "representative genome", "na"))
        flog.info(paste("Subsetting by RefSeq category to include only those within", paste0(wanted_refseq_category, collapse=", "), "categories."))
        wanted_assembly_summary <- subset(wanted_assembly_summary, (refseq_category %in% wanted_refseq_category))
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries matching the RefSeq criteria."))
    }

    if (nobs == TRUE){
        flog.info("Filtering out entries which have been tagged as being problematic or low quality.")
        noproblem <- which(wanted_assembly_summary$excluded_from_refseq %in% c("", "na"))
        wanted_assembly_summary <- wanted_assembly_summary[noproblem, ]
        flog.info(paste("There are", nrow(wanted_assembly_summary), "entries which have NOT been flagged as being problematic or low quality."))
    }

    if (nrow(wanted_assembly_summary) > 0){
        #Fix URLs to match filetype wanted
        filesuffix <- switch(fileformat, "genbank" = "genomic.gbff.gz", "fasta" = "genomic.fna.gz", "gff" = "genomic.gff.gz", "faa" = "protein.faa.gz", "rna" = "rna_from_genomic.fna.gz")
        #Get the repeat bit needed for each ftp_path. Sigh.
        sufreps <- sapply(1:nrow(wanted_assembly_summary), function (x) {unlist(strsplit(wanted_assembly_summary$ftp_path[x], split="/"))[length(unlist(strsplit(wanted_assembly_summary$ftp_path[x], split="/")))]})
        wanted_assembly_summary$url <- paste(wanted_assembly_summary$ftp_path, paste(sufreps, filesuffix, sep="_")  , sep="/")
        wanted_assembly_summary$destfn <- file.path(outputdir, paste(wanted_assembly_summary$organism_name, wanted_assembly_summary$infraspecific_name, wanted_assembly_summary$assembly_accession, filesuffix, sep="_"))
        wanted_assembly_summary$ftp_path <- NULL

        #Curtail to ntop if requested
        if (!is.null(ntop)){
            flog.info(paste("Keeping only the top", ntop, "genomes, as requested with the ntop argument."))
            wanted_assembly_summary <- wanted_assembly_summary[1:ntop, ]
        }

        #Now download stuff wanted to the output dir
        if(simulate == FALSE){
            flog.info("Downloading genomes")
            dnl <- sapply(1:nrow(wanted_assembly_summary), function (x) { tryCatch(download.file(wanted_assembly_summary$url[x], destfile = wanted_assembly_summary$destfn[x] ), error = function(e) print(paste("Unable to download", wanted_assembly_summary$url[x]))) } )
        }

    } else {
        flog.info("No genomes fit the criteria imposed.")
        wanted_assembly_summary <- NULL
    }

    return(wanted_assembly_summary)
}
