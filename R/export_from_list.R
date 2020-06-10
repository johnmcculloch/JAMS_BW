#' export_from_list
#'
#' JAMSalpha function
#' @export

export_from_list <- function(objlist = opt, objnames = NULL, destination = NULL){
    for (objname in objnames){
        wantedobjnum <- which(names(objlist) == objname)
        #Test if object really exists
        if (length(wantedobjnum) > 0){
            wantedobj <- objlist[[wantedobjnum]]
            #Ok, it exists, but what is it
            if (class(wantedobj)[1] == "data.frame"){
                flog.info(paste("Exporting", objname))
                expfn = file.path(destination, paste(paste(opt$prefix, objname, sep = "_"), "tsv", sep = "."))
                write.table(wantedobj, file = expfn, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
            } else {
                flog.info(paste("Unable to export objects like", objname, "using this function for the time being. Stay tuned for newer versions."))
            }
        } else {
            flog.info(paste("There is no such object", objname, "in the list to export."))
        }
    }
}
