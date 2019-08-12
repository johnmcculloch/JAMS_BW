#' fix_interproscanoutput
#'
#' JAMSalpha function
#' @export

fix_interproscanoutput<-function(opt=NULL){

    if("iprodir" %in% names(opt)){
        #Get interprojob
        opt$iprojob <- system2('cat', args=file.path(opt$iprodir, "ipro.job"), stdout = TRUE, stderr = FALSE)
        #See if job finished
        iprojobstatus <- system2('sacct', args=c("-j", opt$iprojob), stdout = TRUE)
        #eliminate header
        iprojobstatus <- iprojobstatus[3:length(iprojobstatus)]
        iprojobstatus <- iprojobstatus[grep("quick", iprojobstatus)]
        totaljobs <- length(iprojobstatus)
        completedjobs <- length(grep("COMPLETED", iprojobstatus))
        ratiojobscomplete <- completedjobs/totaljobs
        runningjobs <- length(grep("RUNNING", iprojobstatus))

        #Delay if there still are jobs to complete
        nattempt=1
        #Wait until jobs are close to complete
        while ((ratiojobscomplete < 0.9) && nattempt < 50) {
            flog.info("Interproscan analysis of proteome is still incomplete.")
            flog.info(paste("There are", runningjobs, " Interpro jobs running for this sample."))
            flog.info(paste("The proportion of Interpro jobs complete is currently", round(ratiojobscomplete, 2)))
            flog.info(paste("Will check again in 5 minutes time. This is attempt", nattempt, "of 50 before giving up."))
            Sys.sleep(300)
            nattempt <- nattempt + 1
            #See if job finished
            iprojobstatus <- system2('sacct', args=c("-j", opt$iprojob), stdout = TRUE)
            #eliminate header
            iprojobstatus <- iprojobstatus[3:length(iprojobstatus)]
            iprojobstatus <- iprojobstatus[grep("quick", iprojobstatus)]
            totaljobs <- length(iprojobstatus)
            completedjobs <- length(grep("COMPLETED", iprojobstatus))
            ratiojobscomplete <- completedjobs/totaljobs
            runningjobs <- length(grep("RUNNING", iprojobstatus))
        }

        flog.info("Harvesting and integrating Interproscan data.")

        interprotsvs<-file.path(opt$iprodir, list.files(path=opt$iprodir, pattern=".tsv"))
        #load tsvs into a single object in memory
        readipro<-function(x){
            read.table(file=x, sep="\t", header=FALSE, quote="", skipNul=FALSE, fill=TRUE, colClasses = "character", col.names=c("Feature","MD5","AALength","Analysis","Accession","Description","Start","Stop","Score","Status","Date","IproAcc","IproDesc","GOterms","Pathways"))
        }
        alliprotsvs <- lapply(interprotsvs, readipro)
        ipro<-plyr::ldply(alliprotsvs, rbind)

        #Clean-up datafrane
        ipro<-subset(ipro, Analysis !="MobiDBLite")
        ipro$MD5<-NULL
        ipro$Start<-NULL
        ipro$Stop<-NULL
        ipro$Start<-NULL
        ipro$Score<-as.numeric(ipro$Score)
        ipro$Score[is.na(ipro$Score)]<- 0
        ipro<-subset(ipro, Score < 0.001)
        ipro<-subset(ipro, Status == "T")
        ipro$Score<-NULL
        ipro$Status<-NULL
        ipro$Date<-NULL
        ipro$AALength<-as.numeric(ipro$AALength)

        #Fill in the blanks
        ipro <- ipro %>% mutate(Description = ifelse(Description == "", "none", Description))
        ipro <- ipro %>% mutate(IproAcc = ifelse(IproAcc == "", "none", IproAcc))
        ipro <- ipro %>% mutate(IproDesc = ifelse(IproDesc == "", "none", IproDesc))
        ipro <- ipro %>% mutate(GOterms = ifelse(GOterms == "", "none", GOterms))
        ipro <- ipro %>% mutate(Pathways = ifelse(Pathways == "", "none", Pathways))

        #remove eventual duplicates
        ipro<-ipro[!duplicated(ipro), ]
        opt$interproscanoutput<-ipro
    }

    return(opt)
}
