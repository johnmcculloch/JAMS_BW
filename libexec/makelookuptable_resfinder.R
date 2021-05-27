library(JAMS)

#Read classes into list
atbclasses <- scan("antibiotic_classes.txt", what="", sep="\n")
atbclasses <- strsplit(atbclasses, "\t")
names(atbclasses) <- sapply(atbclasses, `[[`, 1)
atbclasses <- lapply(atbclasses, `[`, -1)

#Read phenotypes into table
atbphenotypes <- read.table(file = "phenotypes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = '')
colnames(atbphenotypes) <- c("Accession", "Class", "ResistancePhenotype", "PMID", "Mechanism", "Notes", "Required_gene")
atbphenotypes$Accession <- gsub("\'", "", atbphenotypes$Accession)

atbphenotypes$Notes <- NULL
atbphenotypes$Required_gene <- NULL

#Remove dupes
dupes <- atbphenotypes$Accession[duplicated(atbphenotypes$Accession)]
atbphenotypes <- subset(atbphenotypes, (!Accession %in% dupes))

#Fill in blanks with none
for (colm in 1:ncol(atbphenotypes)){
    atbphenotypes[,colm][which(atbphenotypes[,colm] %in% c("", "None"))] <- "none"
}

rownames(atbphenotypes) <- atbphenotypes$Accession

#Get a vector of possible antibiotics
antibiotics <- unlist(strsplit(unique(atbphenotypes$ResistancePhenotype), split = ","))
antibiotics <- trimws(antibiotics)
antibiotics <- unique(antibiotics)
antibiotics <- antibiotics[which(!(antibiotics %in% c("see Notes", "see notes", "unknown", "Unknown", "none")))]
antibiotics <- sort(antibiotics)
Accession <- atbphenotypes$Accession

antibiogram <- matrix(ncol = length(antibiotics), nrow = length(Accession), data = "S")
antibiogram <- as.data.frame(antibiogram, stringsAsFactors = FALSE)
rownames(antibiogram) <- Accession
colnames(antibiogram) <- antibiotics

#Loop round accessions and fix antibiogram
for (accession in rownames(antibiogram)){
    phenointerest <- atbphenotypes[accession, ]
    atbsR <- phenointerest$ResistancePhenotype
    atbsR <- unlist(strsplit(atbsR, split = ","))
    atbsR <- trimws(atbsR)
    atbsR <- unique(atbsR)
    atbsR <- atbsR[which(!(atbsR %in% c("see Notes", "see notes", "unknown", "Unknown", "none")))]
    atbsR <- sort(atbsR)
    antibiogram[accession, atbsR] <- "R"
}

antibiogram$Accession <- rownames(antibiogram)

lookuptable <- left_join(atbphenotypes, antibiogram, by = "Accession")
write.table(lookuptable, file = "lookup.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
