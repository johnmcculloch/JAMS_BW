#' gather_run_info
#'
#' JAMSalpha function
#' @export

gather_run_info<-function(opt=NULL){

    flog.info("Gathering information on JAMS analysis run.")
    stunde<-round(as.numeric((opt$endtime - opt$starttime), units="hours"), 2)

    Run_info=c("JAMS_version", "JAMS_Kdb_Version", "Sample_name", "Run_start", "Run_end", "Duration_hours", "CPUs_used")
    Run_value=c(opt$verstr, opt$JAMS_Kdb_Version, opt$prefix, as.character(opt$starttime), as.character(opt$endtime), stunde, opt$threads)
    projinfo<-data.frame(Run_info=Run_info, Run_value=Run_value)
    
    if("readsdir" %in% names(opt)){
        process<-"Assemble_from_reads"
        Run_info=c("Run_type", "Process", "Sequencing_chemistry", "Library_strategy", "Contig_assembler", "Contig_assembler_version", "Host_species")
        Run_value=c(opt$analysis, process, opt$seqtype, opt$libstructure, opt$assembler, opt$assemblerversion, opt$host)

    } else {
        process<-"Contigs_supplied"
        Run_info=c("Run_type", "Process")
        Run_value=c(opt$analysis, process)
    }

    projinfosup<-data.frame(Run_info=Run_info, Run_value=Run_value)

    opt$projinfo<-rbind(projinfo, projinfosup)

    return(opt)
}
