#' assemble_contigs
#'
#' JAMSalpha function
#' @export

assemble_contigs<-function(opt=NULL){
    #Obtain reads to use for assembly
    opt<-get_reads(opt=opt)
    save.image(file = opt$projimage)
    #Filter reads if necessary
    opt<-trim_reads(opt=opt)
    save.image(file = opt$projimage)
    #Eliminate host reads if applicable
    opt<-filter_host(opt=opt)
    #Check if reads are in tmpdir
    if(opt$workdir != opt$sampledir){
        readsworkdir<-file.path(opt$workdir, "reads")
    } else {
        readsworkdir<-opt$readsdir
    }
    save.image(file = opt$projimage)

    #Set working directory to workdir
    setwd(opt$workdir)

    #Choose input. Assemble with either NAHS or if not present, with trim. 
    if((opt$host=="none") || (opt$analysis %in% c("isolate", "isolaternaseq"))){
        inputreads<-opt$trimreads
    } else {
        inputreads<-opt$nahsreads
    }

    if(opt$analysis %in% c("metatranscriptome", "isolaternaseq")){
        #Disconsider unpaired reads from trimming as an input if RNA.
        if(length(inputreads) > 2){
            inputreads<-inputreads[1:2]
            flog.info("For assembling RNA with SPAdes, unpaired trimmed reads will be not be used.")        
        }
    }

    #Adjust input reads to absolute path
    inputreads<-file.path(readsworkdir, inputreads)

    #Check the files exist or die.
    if(!(all(file.exists(inputreads)))){
        flog.info("Could not find input reads to assemble. Aborting now.")
        q()
    }

    #Estimate number of reads to tag as deep sequencing
    numreadsarg<-c("-c", "^@", inputreads[1])
    estnumreads<-as.numeric(system2('grep', args=numreadsarg, stdout=TRUE))
    numreads<-estnumreads*length(inputreads)
    if(numreads > 10000000){
        opt$deepseq<-TRUE
    } else {
        opt$deepseq<-FALSE
    }

    #If project type is either isolate or rna, make sure the assembler is spades.
    opt$assembler<-switch(opt$analysis,  "metagenome" = opt$assembler, "metatranscriptome" = "spades", "isolate" ="spades", "isolaternaseq" = "spades")
    
    #If read type is Illumina mate-pair library, then default to SPAdes
    if(opt$seqtype %in% c("illuminamp")){
        opt$assembler<-"spades"
    }

    #If analysis is metagenome and reads are single end, switch to megahit.
    if((opt$analysis %in% c("metagenome", "metatranscriptome")) && (length(inputreads) < 2)){
        flog.info("Single end metagenomic reads can only be assembled with MEGAHIT.")
        opt$assembler<-"megahit"
    }
 
    membytes90<-(opt$totmembytes * 0.9)
    assemblercmd<-switch(opt$assembler, "megahit" = "megahit", "spades" = "spades.py")
    opt$assemblerversion<-system2(assemblercmd, args="--version", stdout = TRUE, stderr = TRUE)

    #Build command from options present
    if(opt$assembler == "megahit"){
        #Build arguments for MEGAHIT
        if(opt$deepseq == TRUE){
            presets<-"meta-large"
        } else {
            presets<-"meta-sensitive"
        }

        commonargs<-c("--presets", presets, "-m", membytes90, "--num-cpu-threads", opt$threads, "-o", paste(opt$prefix, "assembly", sep="_"), "--min-contig-len", "500", "--verbose")
        readargs<-list("-r", c("-1", "-2"), c("-1", "-2", "-r"))
        readinput<-as.vector(rbind(readargs[[length(inputreads)]], inputreads))
        assemblerargs<-c(commonargs, readinput)
        asscontigs<-file.path(paste(opt$prefix, "assembly", sep="_"), "final.contigs.fa")
        asslog<-file.path(paste(opt$prefix, "assembly", sep="_"), "log")

    } else if(opt$assembler=="spades"){
        #Build arguments for SPAdes
        commonargs<-c("-t", opt$threads, "-m", (as.integer(opt$totmembytes/1000000000)), "-o", paste(opt$prefix, "assembly", sep="_"))
        if(opt$seqtype == "illuminamp"){
            readargs<-list("-s", c("--hqmp1-1", "--hqmp1-2"), c("--hqmp1-1", "--hqmp1-2", "--hqmp1-s"))
        } else {
            readargs<-list("-s", c("--pe1-1", "--pe1-2"), c("--pe1-1", "--pe1-2", "--pe1-s"))
        }

        readinput<-as.vector(rbind(readargs[[length(inputreads)]], inputreads))
        if(length(inputreads)>1){
            if(opt$seqtype == "illuminamp"){
                flog.info("Sequence read input is Illumina High Quality Mate-Pair type.")
                #readinput<-c("--hqmp1-rf", readinput)            
            } else {
                readinput<-c("--pe1-fr", readinput)
            }
        }

        if(opt$analysis == "isolate"){
            otherargs<-c("--cov-cutoff", "auto")
            #if not too deep, do it carefully
            if(opt$deepseq == FALSE){
                flog.info("Using SPAdes --careful mode for assembly.")
                otherargs<-c("--careful", otherargs)
            }
        } else {
            otherargs<-NULL
        }

        if(opt$analysis == "metagenome"){
            moltypeargs<-c("--meta")
            SPAdescontigsname<-"contigs.fasta"
        } else if(opt$analysis %in% c("metatranscriptome", "isolaternaseq")){
            moltypeargs<-c("--rna")
            SPAdescontigsname<-"transcripts.fasta"
         } else {
            moltypeargs<-NULL
            SPAdescontigsname<-"contigs.fasta"
        }

        if(opt$seqtype=="iontorrent"){
            platformargs<-"--iontorrent"
        } else {
            platformargs<-NULL
        }

        assemblerargs<-c(commonargs, otherargs, moltypeargs, platformargs, readinput)

        asscontigs<-file.path(paste(opt$prefix, "assembly", sep="_"), SPAdescontigsname)
        asslog<-file.path(paste(opt$prefix, "assembly", sep="_"), "spades.log")
    }

    #Launch assembly
    flog.info(paste0("Assembling reads using ", opt$assembler, ". Please be patient..."))
    flog.info(paste0("Assembler command used: ", paste(assemblercmd, paste0(assemblerargs, collapse=" "))))

    system2(assemblercmd, args = assemblerargs, stdout = TRUE, stderr = TRUE)
    flog.info(paste("Assembly of reads using", opt$assembler, "complete."))

    #Copy assembly log to the system
    file.copy(asslog, file.path(opt$workdir, paste(paste(opt$prefix, opt$assembler, sep="_"), "log", sep=".")))

    opt<-prepare_contigs_for_JAMS(opt=opt, fastafile=asscontigs)

    #Bank assembly files to project directory
    NHcontigsfastaout<-paste(paste(opt$prefix, "contigs", sep="_"), "fasta", sep=".")
    write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$workdir, NHcontigsfastaout))

    if(opt$workdir != opt$sampledir){
        #Write contigs from opt directly to system
        write.fasta(sequences = opt$NHcontigs_sequence, names = names(opt$NHcontigs_sequence), nbchar = 80, file.out = file.path(opt$sampledir, NHcontigsfastaout))
        system(paste("cp", asslog, file.path(opt$sampledir, paste(paste(opt$prefix, opt$assembler, sep="_"), "log", sep="."))))
    }

    return(opt)
}
