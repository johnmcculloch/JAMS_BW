#' harvest_functions
#'
#' JAMSalpha function
#' @export

harvest_functions<-function(opt=opt){

    data(ECdescmap)
    data(GOtermdict)

    flog.info("Harvesting functional data.")
    #Harvest common functions
    basicanalyses<-c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb")[c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb") %in% colnames(opt$featuredata)]
    
    featurenumbaseslist<-lapply(basicanalyses, function(x) { compute_signature_numbases(featuredata=opt$featuredata, columnname=x) })
    featurenumbaseslist<-plyr::ldply(featurenumbaseslist, rbind)
    featurenumbaseslist$`.id`<-NULL
    rownames(featurenumbaseslist)<-featurenumbaseslist$Accession

    #Harvest interpro functions, if applicable
    if(opt$skipipro != TRUE){
        opt<-fix_interproscanoutput(opt=opt)
    }
    interpronumbaseslist<-NULL
    if("interproscanoutput" %in% names(opt)){
        opt<-add_interpro_to_featuredata(opt=opt)
        iproanalyses<-sort(unique(opt$interproscanoutput$Analysis))
        iproanalyses<-c(iproanalyses, "Interpro", "GO")
        flog.info("Creating counts table for Interpro signatures.")
        interpronumbaseslist<-lapply(iproanalyses, function(x) { compute_signature_numbases(featuredata=opt$featuredata, columnname=x) })
        names(interpronumbaseslist)<-iproanalyses
        interpronumbaseslist<-plyr::ldply(interpronumbaseslist, rbind)
        interpronumbaseslist$`.id`<-NULL
        rownames(interpronumbaseslist)<-interpronumbaseslist$Accession
    }

    opt$featuredose<-rbind(featurenumbaseslist, interpronumbaseslist)
 
    #Add descriptions to featuredose
    LKTcols<-colnames(opt$featuredose)[4:ncol(opt$featuredose)]
    opt$featuredose$Description<-rep("none", nrow(opt$featuredose))
    #rearrange
    opt$featuredose<-opt$featuredose[, c("Analysis", "Accession", "Description", "NumBases", LKTcols)]
    #Add EC numbers
    descriptions<-ECdescmap$Description[match(opt$featuredose$Accession, ECdescmap$Accession)]
    opt$featuredose$Description[which(!(is.na(descriptions)))]<-descriptions[which(!(is.na(descriptions)))]

    #Add interpro descriptions, if applicable
    if("interproscanoutput" %in% names(opt)){
        #Add GO descriptions
        descriptions<-GOtermdict$TERM[match(opt$featuredose$Accession, GOtermdict$GOID)]
        opt$featuredose$Description[which(!(is.na(descriptions)))]<-descriptions[which(!(is.na(descriptions)))]
        #Make description dictionary
        dictaccessions<-opt$interproscanoutput$Accession
        dictdescriptions<-opt$interproscanoutput$Description
        dictaccessions<-append(dictaccessions, opt$interproscanoutput$IproAcc, after=length(dictaccessions))
        dictdescriptions<-append(dictdescriptions, opt$interproscanoutput$IproDesc, after=length(dictdescriptions))
        acc2desc<-data.frame(Accession=dictaccessions, Description=dictdescriptions, stringsAsFactors=FALSE)
        acc2desc<-acc2desc[!(duplicated(acc2desc)), ]
        #Add interpro descriptions to featuredose
        descriptions<-acc2desc$Description[match(opt$featuredose$Accession, acc2desc$Accession)]
        opt$featuredose$Description[which(!(is.na(descriptions)))]<-descriptions[which(!(is.na(descriptions)))]
    }

    return(opt)
}
