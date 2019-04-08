#' get_feature_to_accession_table_parallel
#'
#' JAMSalpha function
#' @export

get_feature_to_accession_table_parallel<-function(opt=NULL){
    require(parallel)

    iproanalyses<-sort(unique(opt$interproscanoutput$Analysis))
    iproanalyses<-c(iproanalyses, "Interpro", "GO")
    cores<-((opt$threads) - 1)
    ipro_cl <- makeCluster(cores)
    flog.info("Building parallel cluster")
    clusterExport(ipro_cl, list("opt", "iproanalyses", "left_join", "flog.info"))

    flog.info("Adding Interpro signatures to featuredata. Please be patient.")
    #Aggregate accession
    iproanalysislist<-parLapply(ipro_cl, iproanalyses,  function(x) get_feature_to_accession_table(opt=opt, iproanalysis=x))
    names(iproanalysislist)<-iproanalyses
    stopCluster(ipro_cl)

    flog.info("Interpro signatures have been added to featuredata")

    return(iproanalysislist)
}
