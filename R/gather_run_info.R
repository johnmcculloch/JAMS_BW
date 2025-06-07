#' gather_run_info
#'
#' JAMSalpha function
#' @export

gather_run_info <- function(opt = NULL){

    flog.info("Gathering information on JAMS analysis run.")
    useful_parameters <- c("analysis", "verstr", "JAMS_Kdb_Version", "krakenconfidencescore", "prefix", "threads", "contigs", "assemblyaccession", "readorigin", "seqtype", "libstructure", "R1", "R2", "SE", "sraaccession", "assembler", "assemblerversion", "host", "outdir", "maxbases_input", "useipro")

    projinfo <- data.frame(Opt_name = names(unlist(opt[useful_parameters])), Run_value = unlist(opt[useful_parameters]))

    projinfo$Run_info <- sapply(projinfo$Opt_name, function (x) { switch(x, "analysis" = "Run_type", "verstr" = "JAMS_version", "JAMS_Kdb_Version" = "JAMS_Kdb_Version", "krakenconfidencescore" = "JAMS_Kdb_Confidence", "prefix" = "Sample_name", "starttime" = "Run_start", "endtime" = "Run_end", "threads" = "CPUs_used", "contigs" = "Contigs_supplied", "assemblyaccession" = "GenBank_assembly_accession_number", "readorigin" = "Reads_source", "seqtype" = "Read_chemistry", "libstructure" = "Library_strategy", "R1" = "R1_path", "R2" = "R2_path", "SE" = "Single_read_path", "sraaccession" = "NCBI_SRA_accession_number", "assembler" = "Contig_assembler", "assemblerversion" = "Contig_assembler_version", "host" = "Host_species", "outdir" = "Output_path", "maxbases_input" = "Input_number_bp_cap", "useipro" = "InterProScan") } )
    rownames(projinfo) <- projinfo$Run_info

    if ("readsdir" %in% names(opt)){
        process <- "Assemble_from_reads"
    } else {
        process <- "Contigs_supplied"
    }

    projinfosupp <- data.frame(Run_info = c("Run_start", "Run_end"), Opt_name = c("starttime", "endtime"), Run_value = c(as.character(opt$starttime), as.character(opt$endtime)))
    rownames(projinfosupp) <- projinfosupp$Run_info
    stunde <- round(as.numeric((opt$endtime - opt$starttime), units = "hours"), 2)
    projinfosupp["Duration_hours", c("Run_info", "Opt_name", "Run_value")] <- c("Duration_hours", NA, stunde)
    projinfosupp["Process", c("Run_info", "Opt_name", "Run_value")] <- c("Process", NA, process)

    projinfo <- rbind(projinfo, projinfosupp)

    opt$projinfo <- projinfo[ , c("Run_info", "Opt_name", "Run_value")]

    return(opt)
}
