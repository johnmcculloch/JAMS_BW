#' compute_readcounts(opt=opt)
#'
#' JAMSalpha function
#' @export

compute_readcounts<-function(opt=opt){

    #If not re-computing fastqstats, get them from actually counting reads
    if(is.null(opt$fastqstats)){
        #Check if reads are in tmpdir
        if(opt$workdir != opt$sampledir){
            readsworkdir<-file.path(opt$workdir, "reads")
        } else {
            readsworkdir<-opt$readsdir
        }
        setwd(readsworkdir)

        #Copy NAss reads to readsdir for computation, if they exist
        NAssreadsfiles<-list.files(opt$covdir, pattern="_NAss_")
        NAssreadsfilespaths<-file.path(opt$covdir, NAssreadsfiles)
        if(length(NAssreadsfilespaths) > 0){
            file.copy(NAssreadsfilespaths, readsworkdir)
        }

        #Get read stats if applicable
        flog.info("Obtaining read statistics.")
        #Flush out any previous attempts at counting just in case.
        if(file.exists("reads.stats")){
            file.remove("reads.stats")
        }
        countcmd<-file.path(opt$bindir, "countfastqJAMS.sh")
        countargs<-c("-d", readsworkdir, "-t", opt$threads)
        system2(countcmd, countargs, stderr=FALSE, stdout=FALSE)
        opt$fastqstats<-read.table("reads.stats", header=FALSE, sep="\t", stringsAsFactors=FALSE, skipNul=FALSE, fill=TRUE, colClasses=c("character","numeric","numeric"))
        colnames(opt$fastqstats)<-c("Reads", "Count", "Bases")
    } else {
        NAssreadsfiles<-opt$fastqstats[grep("_NAss_",opt$fastqstats$Reads ),"Reads"]
    }

    #Flush out any lines which are NA
    opt$fastqstats<-opt$fastqstats[!is.na(opt$fastqstats[,2]), ]
    opt$fastqstats<-opt$fastqstats[!is.na(opt$fastqstats[,3]), ]

    rawcount<-sum(subset(opt$fastqstats, Reads %in% opt$rawreads)[]$Count)
    rawbases<-sum(subset(opt$fastqstats, Reads %in% opt$rawreads)[]$Bases)
    trimcount<-sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Count)
    trimbases<-sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Bases)
    NAsscount<-sum(subset(opt$fastqstats, Reads %in% NAssreadsfiles)[]$Count)
    NAssbases<-sum(subset(opt$fastqstats, Reads %in% NAssreadsfiles)[]$Bases)

    #assbases<-(sum(opt$contigsdata$NumBases))

    trimbasestats<-round(100*(c((trimbases/rawbases), 1-(trimbases/rawbases))), 1)

    if (opt$host == "none" ){
        #Compute number of bases which were assembled. Assembly was donne from trimreads
        assbases<-trimbases-NAssbases
        assbasestats<-round(100*(c((assbases/trimbases), 1-(assbases/trimbases))), 1)
        readsdf<-data.frame(Condition=c("Surviving", "Not-surviving"), QC=trimbasestats, Ass=assbasestats)
        opt$readstats<-data.frame(Read_type=c("Raw", "Trimmed", "Assembled"), Base_counts=c(rawbases, trimbases, assbases))
    } else {
        NAHScount<-sum(subset(opt$fastqstats, Reads %in% opt$nahsreads)[]$Count)
        NAHSbases<-sum(subset(opt$fastqstats, Reads %in% opt$nahsreads)[]$Bases)
        NAHSbasestats<-round(100*(c((NAHSbases/trimbases), 1-(NAHSbases/trimbases))), 1)
        #Compute number of bases which were assembled. Assembly was donne from NAHSreads
        assbases<-NAHSbases-NAssbases
        assbasestats<-round(100*(c((assbases/NAHSbases), 1-(assbases/NAHSbases))), 1)
        readsdf<-data.frame(Condition=c("Surviving", "Not-surviving"), QC=trimbasestats, NAHS=NAHSbasestats, Ass=assbasestats)
        opt$readstats<-data.frame(Read_type=c("Raw", "Trimmed", "NonHost", "Assembled"), Base_counts=c(rawbases, trimbases, NAHSbases, assbases))
    }

    rownames(opt$readstats)<-opt$readstats$Read_type

    readsdf$survivingQC<-paste(paste(readsdf$Condition, readsdf$QC, sep="="), "%", sep="")
    p <- ggplot(readsdf, aes(x = "", y = QC, fill = survivingQC)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("red", "green"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Passing QC", sep=" - "))
    p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
    pieQC<-p

   if ( opt$host != "none" ){
        readsdf$Alignment_to_Host<-paste(paste(c("Not Aligned", "Aligned"), readsdf$NAHS, sep="="), "%", sep="")
        p <- ggplot(readsdf, aes(x = "", y = NAHS, fill = Alignment_to_Host)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("red", "green"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Aligned to Host Species", sep=" - "))
        p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
        pieNAHS<-p
        flog.info(paste("Host contamination rate is", subset(readsdf, Condition=="Not-surviving")[]$NAHS, "%"))
    } else {
        pieNAHS<-NULL
    }

    readsdf$Assembled_bases<-paste(paste(c("Assembled", "Not Assembled"), readsdf$Ass, sep="="), "%", sep="")
    p <- ggplot(readsdf, aes(x = "", y = Ass, fill = Assembled_bases)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("green", "red"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Assembled into Contigs", sep=" - "))
    p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
    pieAss<-p
    flog.info(paste("Contig assembly rate is", subset(readsdf, Condition=="Surviving")[]$Ass, "%"))

    if (opt$host != "none" ){
        GBases<-round((c(rawbases, trimbases, NAHSbases, assbases))/1000000000, 2)
        Read_type<-c("Raw Reads", "Trimmed Reads", "Non-Aligned to Host", "Assembled")
        Read_type<-factor(Read_type, levels = Read_type)
        readsdf2<-data.frame(Read_type=Read_type, Gigabases=GBases)
    }else{
        GBases<-round((c(rawbases, trimbases, assbases))/1000000000, 2)
        RLs<-round(c((rawbases/rawcount), (trimbases/trimcount), (assbases/NAsscount)), 1)
        Read_type<-c("Raw Reads", "Trimmed Reads", "Assembled")
        Read_type<-factor(Read_type, levels = Read_type)
        readsdf2<-data.frame(Read_type=Read_type, Gigabases=GBases)
    }

    Bar.gigabases <- ggplot(readsdf2, aes(Read_type, GBases)) + geom_bar(stat = "identity", fill=(rep(4, length(readsdf2[,1])))) + labs(title = paste(opt$prefix, "Total Basepairs per Read Type (Gbp)", sep=" - "))

    opt$readplots<-list(pieQC, pieNAHS, pieAss, Bar.gigabases)

    return(opt)
}
