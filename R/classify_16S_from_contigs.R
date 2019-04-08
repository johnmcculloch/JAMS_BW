#' classify_16S_from_contigs(opt=opt)
#' 
#' Classifies taxonomically 16S features found in contig annotations using Dada2
#' @export

classify_16S_from_contigs<-function(opt=NULL){
    setwd(opt$workdir)
    #Find 16S features to annotate
    #Get only rRNAs
    featuresofinterest<-opt$featuredata[(opt$featuredata$FeatType =="rRNA"), ]
    #Get only 16S features
    features16S<-featuresofinterest[grep("16S_ribosomal_RNA", featuresofinterest$Product), "Feature"]

    #Set 16S database to use
    db16S<-file.path(opt$dbdir, "16Sdbs")
    #make sure to only pick one if there are several different versions in the same directory.
    db16S_trainset<-sort(file.path(db16S, list.files(path=db16S, pattern="train")))[1]
    db16S_species<-sort(file.path(db16S, list.files(path=db16S, pattern="species")))[1]

    if(length(features16S) > 0){

        #Write 16S data to system for Dada2
        features16Ssequences<-filter_sequence_by_name(input_sequences=opt$genes, sequencenames=features16S, keep=TRUE)

        #make sure there are no ambiguous bases, else dada2 will stall.
        validbases<-c("A", "C", "T", "G")
        ambiguities<-sapply(1:length(features16Ssequences), function (x) { length(grep(paste(validbases, collapse="|"), as.character(features16Ssequences[[x]]), ignore.case = TRUE, invert = TRUE)) })
        offenders<-which(ambiguities != 0)
        if(length(offenders) > 0){
            ambiguousseqs<-names(features16Ssequences)[offenders]
            flog.info(paste("Skipping 16S rRNA sequences", paste(ambiguousseqs, collapse=", "), " because these contain ambiguous bases and cannot be identified by exact matches using dada2. "))
            features16Ssequences<-features16Ssequences[!(names(features16Ssequences) %in% ambiguousseqs)]
        }
        write.fasta(sequences = features16Ssequences, names = names(features16Ssequences), nbchar = 80, file.out = paste(opt$prefix, "16SrRNA.fna", sep= "_"))

        set.seed(100)
        #Classify from Phylum to Genus
        flog.info(paste("Found", length(features16Ssequences), "16S rRNA sequences in contigs. Classifying them via dada2 using database", db16S_species))
        flog.info("Step 1 - Classifying down to Genus level using a Naive Bayesian Classifier")
        dadataxa <- assignTaxonomy(paste(opt$prefix, "16SrRNA.fna", sep= "_"), db16S_trainset, multithread=TRUE, minBoot=80)
        taxa_16S<-as.data.frame(unname(dadataxa))
        taxa_16S[] <- lapply(taxa_16S, as.character)
        colnames(taxa_16S)<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
        taxa_16S[is.na(taxa_16S)] <- "Unclassified"
        taxa_16S$Feature<-names(features16Ssequences)

        #Classify Genus and species
        flog.info("Step 2 - Classifying down to Species level using exact matches.")
        taxa_16S_Gs<-assignSpecies(paste(opt$prefix, "16SrRNA.fna", sep= "_"), db16S_species)
        taxa_16S_Gs<-as.data.frame(unname(taxa_16S_Gs))
        taxa_16S_Gs[] <- lapply(taxa_16S_Gs, as.character)
        colnames(taxa_16S_Gs)<-c("Genus", "Species")
        taxa_16S_Gs[is.na(taxa_16S_Gs)] <- "Unclassified"
        taxa_16S_Gs$Feature<-names(features16Ssequences)

        #Join the tables
        taxa_16S_cons<-taxa_16S
        taxa_16S_cons$GenusGs<-taxa_16S_Gs$Genus[match(taxa_16S$Feature, taxa_16S_Gs$Feature)]
        taxa_16S_cons$Species<-taxa_16S_Gs$Species[match(taxa_16S$Feature, taxa_16S_Gs$Feature)]
        taxa_16S_cons[] <- lapply(taxa_16S_cons, as.character)

        #Keep Genus level from trainset if genus from Genus+Species db is Unclassified. 
        for (g in 1:length(taxa_16S_cons$GenusGs)){
            taxa_16S_cons$GenusGs[g]<-ifelse((taxa_16S_cons$GenusGs[g] == "Unclassified") && (taxa_16S_cons$Genus[g] != "Unclassified"),  taxa_16S_cons$Genus[g], taxa_16S_cons$GenusGs[g])
        }
        taxa_16S_cons$Genus<-NULL
        colnames(taxa_16S_cons)[which(colnames(taxa_16S_cons)=="GenusGs")]<-"Genus"

        #Reorder consolidated data frame
        taxlvls<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
        taxa_16S_cons<-taxa_16S_cons[c("Feature", taxlvls, "Species")]

        taxa_16S_cons$Genus<-unlist(lapply(1:length(taxa_16S_cons$Genus), function(x) { strsplit(taxa_16S_cons$Genus, split="_")[[x]][1] }))

        #Add level tags
        taxtags<-c("d__", "p__", "c__", "o__", "f__", "g__")

        for (t in 1:length(taxlvls)){
            taxa_16S_cons[ which(taxa_16S_cons[, which(colnames(taxa_16S_cons)==taxlvls[t])] != "Unclassified"), which(colnames(taxa_16S_cons)==taxlvls[t])]<-paste(taxtags[t], taxa_16S_cons[ which(taxa_16S_cons[, which(colnames(taxa_16S_cons)==taxlvls[t])] != "Unclassified"), which(colnames(taxa_16S_cons)==taxlvls[t])], sep="")
        }

        taxa_16S_cons<-infer_LKT(taxa_16S_cons)

        opt$featuredata$SixteenSid<-rep("none", nrow(opt$featuredata))
        opt$featuredata$SixteenSid[match(taxa_16S_cons$Feature, opt$featuredata$Feature)]<-taxa_16S_cons$LKT
 
        #Bank 16S files to project directory
        if(opt$workdir != opt$sampledir){
            #Write 16S features to system
            write.fasta(sequences = features16Ssequences, names = names(features16Ssequences), nbchar = 80, file.out = file.path(opt$sampledir, paste(opt$prefix, "16SrRNA.fna", sep= "_")))
        }

    } else {
        flog.info("No 16S rRNA features found in contigs.")
    }

    return(opt)
}
