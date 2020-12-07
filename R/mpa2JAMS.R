#' mpa2JAMS(mpa_folder = NULL, only_prefixes = NULL, output_folder = NULL, mpa_file_suffix = "_profile.txt", return_mpa_list.data = TRUE, verbose = TRUE)
#'
#' Makes .jams files from Metaphlann output files. This function is experimental, use at your own risk and peril.
#' @export

mpa2JAMS <- function(mpa_folder = NULL, only_prefixes = NULL, output_folder = NULL, mpa_file_suffix = "_profile.txt", return_mpa_list.data = TRUE, verbose = TRUE, asRDS = TRUE){

    flog.warn("This function is experimental, use at your own risk and peril.")
    mpa_list.data <- list()
    list_elem <- 1
    mpa_taxlevel_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    mpa_folder <- fixrelpath(mpa_folder)
    mpa_files_absolute <- list.files(path = mpa_folder, pattern = mpa_file_suffix, full.names = TRUE, recursive = TRUE)
    mpa_files <- sapply(1:length(mpa_files_absolute), function(x) { tail(unlist(strsplit(mpa_files_absolute[x], split = "/")), n = 1) })

    if (length(mpa_files) < 1){
        stop(paste("No MetaPhlAnn output files were found"))
    }

    mpa_files_df <- data.frame(MPA_path = mpa_files_absolute, MPA_file = mpa_files, stringsAsFactors = FALSE)

    mpa_files_df$Prefix <- gsub(mpa_file_suffix, "", mpa_files_df$MPA_file)

    #Load only prefixes requested
    if (!is.null(only_prefixes)){
        mpa_files_df <- subset(mpa_files_df, Prefix %in% only_prefixes)
        if (nrow(mpa_files_df) < 1){
            flog.warn("Metaphlann files for none of the requested prefixes were found in the path supplied.")

            return(NULL)
        }
    }

    #Loop over files and load data into mpa_list.data
    for(mpan in 1:nrow(mpa_files_df)){
        #Load Metaphlann output into R
        curr_prefix <- mpa_files_df$Prefix[mpan]
        if (verbose){
            flog.info(paste("Processing sample", curr_prefix))
        }
        curr_mpa_output <- read.table(file = mpa_files_df$MPA_path[mpan], header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
        curr_mpa_output <- curr_mpa_output[ , 1:3]
        colnames(curr_mpa_output)<- c("Taxinfo", "Taxids", "PctRelabund")

        #Find how many taxlevels are present
        curr_mpa_output$NumTaxLevels <- sapply(1:nrow(curr_mpa_output), function (x) { length(unlist(strsplit(curr_mpa_output[x, "Taxinfo"], split = "\\|"))) } )

        #Consider LKT to be the level which has the largest number of levels
        maxlvlnum <- max(curr_mpa_output$NumTaxLevels)

        #Subset the relabund table to that
        curr_mpa_output <- subset(curr_mpa_output, NumTaxLevels == maxlvlnum)

        for(curr_mpa_taxlevel in 1:maxlvlnum){
            curr_mpa_output[ , mpa_taxlevel_names[curr_mpa_taxlevel]] <- sapply(1:nrow(curr_mpa_output), function(x) { unlist(strsplit( curr_mpa_output[x, "Taxinfo"], split = "\\|"))[curr_mpa_taxlevel] } )
        }

        #Add final Taxid to a column
        curr_mpa_output$Taxid <- sapply(1:nrow(curr_mpa_output), function (x) { tail(unlist(strsplit(curr_mpa_output$Taxids[x], split = "\\|")), n = 1) } )

        curr_mpa_output$LKT <- paste("LKT", curr_mpa_output[ , mpa_taxlevel_names[maxlvlnum]], sep = "__")
        curr_mpa_output$Domain <- curr_mpa_output$Kingdom

        #Make a JAMS-style LKTdose output
        taxlevelspresent <- c("Domain", mpa_taxlevel_names, "LKT") [c("Domain", mpa_taxlevel_names, "LKT") %in% colnames(curr_mpa_output)]
        LKTdose <- curr_mpa_output[ , c(taxlevelspresent, "PctRelabund")]
        colnames(LKTdose)[which(colnames(LKTdose) == "PctRelabund")] <- "NumBases"
        LKTdose$PctFromCtg <- 100.00
        LKTdose$ProbNumGenomes <- 1.00
        LKTdose$NumBases <- (LKTdose$NumBases * 10000)

        mpa_list.data[[list_elem]] <- LKTdose
        names(mpa_list.data)[list_elem] <- paste(curr_prefix, "LKTdose", sep = "_")
        list_elem <- list_elem + 1

        #Make sampleinfo df
        mpa_dbase <- system2("head", args = c("-1", mpa_files_df$MPA_path[mpan]), stdout = TRUE, stderr = FALSE)
        mpa_dbase <- gsub("#", "", mpa_dbase)
        projinfo <- data.frame(Run_info = c("Sample_name", "Run_type", "Process", "Metaphlann_Version", "JAMS_Kdb_Version"), Run_value = c(curr_prefix, "metagenome", "Metaphlann", mpa_dbase, mpa_dbase), stringsAsFactors = FALSE)
        rownames(projinfo) <- projinfo$Run_info

        mpa_list.data[[list_elem]] <- projinfo
        names(mpa_list.data)[list_elem] <- paste(curr_prefix, "projinfo", sep = "_")
        list_elem <- list_elem + 1
    }

    if (!is.null(output_folder)){
        output_folder <- fixrelpath(output_folder)
        flog.info(paste("Writing .jams files to", output_folder))
        dir.create(output_folder, recursive = TRUE)

        current_dir <- getwd()
        setwd(output_folder)
        dir.create("jamstempfiles")
        setwd("jamstempfiles")

        if (asRDS){
            filesuffix <- "rds"
        } else {
            filesuffix <- "tsv"
        }

        #loop through prefixes and write tsvs to temp folder
        for (prefix in mpa_files_df$Prefix){
            for (objname in c("LKTdose", "projinfo")){
                wantedobj <- paste(prefix, objname, sep="_")
                expfn <- file.path(getwd(), paste(wantedobj, filesuffix, sep = "."))
                if (asRDS){
                    saveRDS(mpa_list.data[[wantedobj]], file = expfn)
                } else {
                    write.table(mpa_list.data[[wantedobj]], file = expfn, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
                }
            }
            #Tar them up
            JAMStsvs <- paste(paste(prefix, c("LKTdose", "projinfo"), sep = "_"), filesuffix, sep = ".")
            jamsfile <- paste(prefix, "jams", sep=".")
            jamsargs <- c("-zcvf", jamsfile, JAMStsvs)
            system2('tar', args = jamsargs, stdout = FALSE, stderr = FALSE)
            #move jamsfiles to destination
            system2('mv', args = c(jamsfile, output_folder), stdout = FALSE, stderr = FALSE)
            #delete current tsvs
            file.remove(JAMStsvs)
        }
        #delete jamstempfolder
        setwd(output_folder)
        unlink("jamstempfiles", recursive = TRUE)
        #Back to where we were before
        setwd(current_dir)
    }

    if (return_mpa_list.data == TRUE){
        return(mpa_list.data)
    }

}
