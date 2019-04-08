#' kraken_classify_taxonomy
#'
#' JAMSalpha function
#' @export

kraken_classify_taxonomy<-function(opt=NULL, fastafile=NULL, usekraken1=FALSE){

        if(missing(usekraken1) || (usekraken1 != TRUE)){
            #Default is kraken2
            #Count kmers
            flog.info("Counting k-mers of query sequences with kraken2")
            krakenargs<-c("--db", opt$workingkrakendb, "--threads", opt$threads, fastafile)
            kraken2taxid<-system2("kraken2", args=krakenargs, stdout = TRUE, stderr = FALSE)
            kraken2taxid<-strsplit(kraken2taxid, split="[\t]", fixed = FALSE)
            k2out<-plyr::ldply(kraken2taxid, rbind)
            colnames(k2out)<-c("ClassFlag", "Sequence", "Taxid", "Length", "kmers")
            k2out[] <- lapply(k2out, as.character)
            JAMStaxtablefile<-file.path(opt$workingkrakendb, "JAMS_taxtable.tsv")
            if(file.exists(JAMStaxtablefile)){
                JAMStaxtable<-read.table(file=JAMStaxtablefile, sep="\t", fill=TRUE, stringsAsFactors = FALSE, quote=NULL, header=TRUE)
            } else {
                #Fall back on generic taxonomy table and warn user
                flog.info("JAMS taxonomy table not found. Falling back on generic JAMS taxtable.")
                data(JAMStaxtable)
            }
            JAMStaxtable[] <- lapply(JAMStaxtable, as.character)
            krakendf<-k2out[,c("Sequence", "Taxid")]
            krakendf<-left_join(krakendf, JAMStaxtable)
            #Temporarily, fill in taxids which are NOT in the database with unclassifieds
            unclassdf<-JAMStaxtable[which(JAMStaxtable$Taxid == 0), 2:ncol(JAMStaxtable)]
            krakendf[which(is.na(krakendf$Domain==TRUE)), colnames(unclassdf)]<-unname(unclassdf)

        } else if (usekraken1 == TRUE) {
            #Use kraken1 specifically (deprecated)
            #Count kmers
            flog.info("Counting k-mers of query sequences with kraken")
            krakenargs<-c("--preload", "--db", opt$workingkrakendb, "--threads", opt$threads, fastafile, "> tmp.kraken")
            system2("kraken", args=krakenargs, stdout = TRUE, stderr = TRUE)

            #Compute labels
            flog.info("Naming taxonomy of query sequences with kraken-translate")
            krakenargs<-c("--db", opt$workingkrakendb, "--mpa-format", "tmp.kraken")
            krakenlabels<-system2("kraken-translate", args=krakenargs, stdout = TRUE, stderr = FALSE)
            krakenlabels<-strsplit(krakenlabels, split="[\\||\t]", fixed = FALSE)
            krakendf<-plyr::ldply(krakenlabels, rbind)

            #Obtain contig names to account for unclassified contigs
            flog.info("Accounting for unclassified reads")
            contignames<-system2('cat', args=c(fastafile, "|", "grep", "^\\>", "|", "tr", "-d", "\\>"), stdout = TRUE, stderr = FALSE)
            unclassctgs<-contignames[which(!(contignames %in% as.vector(krakendf[,1])))]
            unclassctgdata<-as.data.frame(matrix(ncol=ncol(krakendf), nrow=length(unclassctgs), data=NA))
            unclassctgdata[,1]<-unclassctgs
            colnames(krakendf)<-paste0("V", 1:ncol(krakendf))
            krakendf<-rbind(krakendf, unclassctgdata)
            flog.info("Fixing kraken output to JAMS format")
            krakendf<-fix_kraken(taxtable=krakendf, keepstrain=FALSE)

            #cleanup
            file.remove("tmp.kraken")
        }

        return(krakendf)
}

