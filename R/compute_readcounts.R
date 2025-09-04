#' compute_readcounts(opt=opt)
#'
#' JAMSalpha function
#' @export

compute_readcounts <- function(opt = opt){

    ##Aggregate counts
    #Start by aggregating counts which will always be present
    rawbases <- sum(subset(opt$fastqstats, Reads %in% opt$rawreads)[]$Bases)
    NAssbases <- sum(subset(opt$fastqstats, Reads %in% opt$NAssfastqs)[]$Bases)

    #If reads were subsampled, then take that into consideration
    #If you are reading this code: I know this is a verbose way of accounting, but between ambiguity and prolixity, the latter is the lesser of two evils and I am judging clarity to be most important, so please bear with me and be forgiving.
    if (!is.null(opt$totbasesbeforesubsampling)){
        flog.info("Reads used for assembly in this sample have been subsampled. Taking this into account for computing read statistics.")
        if (opt$host == "none" ){
            #Trimmed reads were the ones subsampled
            trimcount <- sum(opt$totbasesbeforesubsampling$Count)
            trimbases <- sum(opt$totbasesbeforesubsampling$Bases)
            subsampledcount <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Count)
            subsampledbases <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Bases)

            #Compute number of bases which were assembled. Assembly was done from subsampledtrimreads
            trimbasestats <- round(100 * (c((trimbases / rawbases), 1 - (trimbases / rawbases))), 1)
            assbases <- subsampledbases - NAssbases
            assbasestats <- round(100 * (c((assbases / subsampledbases), 1 - (assbases / subsampledbases))), 1)
            subsampledstats <- round(100 * (c((subsampledbases / trimbases), 1 - (subsampledbases / trimbases))), 1)
            readsdf <- data.frame(Condition = c("Surviving", "Not-surviving"), QC = trimbasestats, Subsampled = subsampledstats, Ass = assbasestats)
            opt$readstats <- data.frame(Read_type = c("Raw", "Trimmed", "Subsampled", "Assembled"), Base_counts = c(rawbases, trimbases, subsampledbases, assbases), stringsAsFactors = FALSE)

        } else {
            #NAHS reads were the ones subsampled
            trimcount <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Count)
            trimbases <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Bases)
            NAHSbases <- sum(opt$totbasesbeforesubsampling$Bases)
            NAHSbasestats <- round(100 * (c((NAHSbases / trimbases), 1 - (NAHSbases / trimbases))), 1)
            subsampledcount <- sum(subset(opt$fastqstats, Reads %in% opt$nahsreads)[]$Count)
            subsampledbases <- sum(subset(opt$fastqstats, Reads %in% opt$nahsreads)[]$Bases)

            #Compute number of bases which were assembled. Assembly was done from subsampledNAHSreads
            trimbasestats <- round(100 * (c((trimbases / rawbases), 1 - (trimbases / rawbases))), 1)
            assbases <- subsampledbases - NAssbases
            assbasestats <- round(100 * (c((assbases / subsampledbases), 1 - (assbases / subsampledbases))), 1)
            subsampledstats <- round(100 * (c((subsampledbases / NAHSbases), 1 - (subsampledbases / NAHSbases))), 1)
            readsdf <- data.frame(Condition = c("Surviving", "Not-surviving"), QC = trimbasestats, NAHS = NAHSbasestats, Subsampled = subsampledstats, Ass = assbasestats)
            opt$readstats <- data.frame(Read_type = c("Raw", "Trimmed", "NonHost", "Subsampled", "Assembled"), Base_counts = c(rawbases, trimbases, NAHSbases, subsampledbases, assbases), stringsAsFactors = FALSE)
        }
    } else {
        #No subsampling happened
        if (opt$host == "none" ){
            trimcount <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Count)
            trimbases <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Bases)

            #Compute number of bases which were assembled. Assembly was done from trimreads
            trimbasestats <- round(100 * (c((trimbases / rawbases), 1 - (trimbases / rawbases))), 1)
            assbases <- trimbases - NAssbases
            assbasestats <- round(100 * (c((assbases / trimbases), 1 - (assbases / trimbases))), 1)
            readsdf <- data.frame(Condition = c("Surviving", "Not-surviving"), QC = trimbasestats, Ass = assbasestats)
            opt$readstats <- data.frame(Read_type = c("Raw", "Trimmed", "Assembled"), Base_counts = c(rawbases, trimbases, assbases), stringsAsFactors = FALSE)

        } else {
            trimcount <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Count)
            trimbases <- sum(subset(opt$fastqstats, Reads %in% opt$trimreads)[]$Bases)
            NAHSbases <- sum(subset(opt$fastqstats, Reads %in% opt$nahsreads)[]$Bases)
            NAHSbasestats <- round(100 * (c((NAHSbases / trimbases), 1 - (NAHSbases / trimbases))), 1)

            #Compute number of bases which were assembled. Assembly was done from NAHSreads
            trimbasestats <- round(100 * (c((trimbases / rawbases), 1 - (trimbases / rawbases))), 1)
            assbases <- NAHSbases - NAssbases
            assbasestats <- round(100 * (c((assbases / NAHSbases), 1 - (assbases / NAHSbases))), 1)
            readsdf <- data.frame(Condition = c("Surviving", "Not-surviving"), QC = trimbasestats, NAHS = NAHSbasestats, Ass = assbasestats)
            opt$readstats <- data.frame(Read_type = c("Raw", "Trimmed", "NonHost", "Assembled"), Base_counts = c(rawbases, trimbases, NAHSbases, assbases), stringsAsFactors = FALSE)
        }
    }

    rownames(opt$readstats) <- opt$readstats$Read_type
    readsdf$survivingQC <- paste(paste(readsdf$Condition, readsdf$QC, sep = "="), "%", sep = "")
    p <- ggplot(readsdf, aes(x = "", y = QC, fill = survivingQC)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("red", "green"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Passing QC", sep=" - "))
    p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))
    pieQC <- p

   if (opt$host != "none"){
        readsdf$Alignment_to_Host <- paste(paste(c("Not Aligned", "Aligned"), readsdf$NAHS, sep = "="), "%", sep = "")
        p <- ggplot(readsdf, aes(x = "", y = NAHS, fill = Alignment_to_Host)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("red", "green"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Aligned to Host Species", sep=" - "))
        p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 14, face = "bold"))
        pieNAHS <- p
        flog.info(paste("Host contamination rate is", subset(readsdf, Condition == "Not-surviving")[]$NAHS, "%"))
    } else {
        pieNAHS <- NULL
    }

    if ("Subsampled" %in% colnames(readsdf)){
         readsdf$Subsampling <- paste(paste(c("Kept", "Eliminated"), readsdf$Subsampled, sep = "="), "%", sep = "")
         p <- ggplot(readsdf, aes(x = "", y = Subsampled, fill = Subsampling)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("green", "red"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Kept after Subsampling due to Assembly Input Cap", sep = " - "))
         p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 14, face = "bold"))
         pieSubsampling <- p
         #flog.info(paste0("After subsampling ", subset(readsdf, Condition == "Not-surviving")[]$Subsampled, "%"))
     } else {
         pieSubsampling <- NULL
     }

    readsdf$Assembled_bases <- paste(paste(c("Assembled", "Not Assembled"), readsdf$Ass, sep = "="), "%", sep = "")
    p <- ggplot(readsdf, aes(x = "", y = Ass, fill = Assembled_bases)) + geom_bar(width = 1, stat = "identity") + scale_fill_manual(values = c("green", "red"))  + coord_polar("y") + labs(title = paste(opt$prefix, "% Bases Assembled into Contigs", sep=" - "))
    p <- p + theme( axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 14, face = "bold"))
    pieAss <- p
    flog.info(paste("Contig assembly rate is", subset(readsdf, Condition == "Surviving")[]$Ass, "%"))

    #Transform counts to Giga Basepairs for bar plotting.
    if (opt$host != "none" ){
        if ("Subsampled" %in% colnames(readsdf)){
            GBases <- round((c(rawbases, trimbases, NAHSbases, subsampledbases, assbases)) / 1000000000, 2)
            Read_type <- c("Raw Reads", "Trimmed Reads", "Non-Aligned to Host", "Subsampled Reads", "Assembled")
        } else {
            GBases <- round((c(rawbases, trimbases, NAHSbases, assbases)) / 1000000000, 2)
            Read_type <- c("Raw Reads", "Trimmed Reads", "Non-Aligned to Host", "Assembled")
        }
    } else {
        if ("Subsampled" %in% colnames(readsdf)){
            GBases <- round((c(rawbases, trimbases, subsampledbases, assbases)) / 1000000000, 2)
            Read_type <- c("Raw Reads", "Trimmed Reads", "Subsampled Reads", "Assembled")
        } else {
            GBases <- round((c(rawbases, trimbases, assbases)) / 1000000000, 2)
            Read_type <- c("Raw Reads", "Trimmed Reads", "Assembled")
        }
    }
    Read_type <- factor(Read_type, levels = Read_type)
    readsdf2 <- data.frame(Read_type = Read_type, Gigabases = GBases)

    Bar.gigabases <- ggplot(readsdf2, aes(Read_type, GBases)) + geom_bar(stat = "identity", fill = "blue") + labs(title = paste(opt$prefix, "Total Basepairs per Read Type (Gbp)", sep = " - "))

    opt$readplots <- list(pieQC, pieNAHS, pieSubsampling, pieAss, Bar.gigabases)
    #Flush out empty elements from list
    opt$readplots <- opt$readplots[sapply(opt$readplots, function(x){ !(is.null(x)) })]

    return(opt)
}
