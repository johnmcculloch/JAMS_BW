#' filter_host
#'
#' JAMSalpha function
#' @export

filter_host <- function(opt = NULL){
    if(opt$host != "none"){

        #Check if reads are in tmpdir
        if (opt$workdir != opt$sampledir){
            readsworkdir <- file.path(opt$workdir, "reads")
        } else {
            readsworkdir <- opt$readsdir
        }
        setwd(readsworkdir)

        if (!(all(file.exists(file.path(readsworkdir, opt$trimreads))))){
            flog.info("Could not find trimmed reads to align to host. Aborting now")
            q()
        }

        #Build an index if explicit species specified
        if (!(is.null(opt$host_accession_number))){
            dir.create("hostindex")
            opt$indexpath <- file.path(readsworkdir, "hostindex")
            setwd(opt$indexpath)
            #Download host genome
            flog.info(paste("Downloading", opt$host, "genome."))
            #Trying
            opt$host_assembly <- get_genomes_NCBI(organisms = NULL, assembly_accession = opt$host_accession_number, outputdir = opt$indexpath, ntop = 1, fileformat = "fasta", simulate = FALSE)

            #Gunzip it
            flog.info("Decompressing downloaded host genome.")
            system(paste("pigz -d", opt$host_assembly$destfn))
            hostfasta <- list.files(pattern="fna")
            flog.info("Building host genome index. Please be patient.")
            buildlog <- system2("bowtie2-build", args=c("-f", hostfasta, "--threads", opt$threads, opt$hostspecies), stderr=FALSE, stdout=TRUE)
            index <- file.path(opt$indexpath, opt$hostspecies)
            setwd(readsworkdir)
        } else {
            #Index is supplied
            #If on a tmppath, copy host genome to tempreadsir for speed
            if(opt$workdir != opt$sampledir){
                flog.info("Copying host genome indices to temporary directory.")
                file.copy(opt$relindexfiles, readsworkdir)
                index<-file.path(readsworkdir, opt$hostspecies)
            } else {
                index<-file.path(opt$indexpath, opt$hostspecies)
            }
        }

        #Set general options
        if(nchar(as.character(Sys.getenv("SLURM_JOB_ID"))) > 4){
            bowtiethreads<-(opt$threads-2) #Two threads must be free on Biowulf
        } else {
            bowtiethreads<-opt$threads
        }

        myfastqs<-opt$trimreads

        flog.info(paste("Aligning trimmed reads to the", opt$host, "genome."))
        #Build command from options present
        bowtiecommonargs<-c("-q", "--very-fast-local", "--ignore-quals", "--mm", "--threads", bowtiethreads, "-x", index)
        bowtiereadargs<-list("-U", c("-1", "-2"), c("-1", "-2", "-U"))
        bowtiereadinput<-as.vector(rbind(bowtiereadargs[[length(myfastqs)]], myfastqs))
        bowtiereadoutput<-c("--un", (paste(opt$prefix, "NAHS_SE.fastq", sep="_")), "--un-conc", (paste(opt$prefix, "NAHS_PE.fastq", sep="_")))[1:(2*length(opt$rawreads))]
        bowtiesam<-c("-S", "RvsHS.sam")
        bowtieargs<-c(bowtiecommonargs, bowtiereadinput, bowtiereadoutput, bowtiesam)
        opt$bowtie_host_output<-system2('bowtie2', args = bowtieargs, stdout = TRUE, stderr = TRUE)

        flog.info(paste("Alignment to host completed with", opt$bowtie_host_output[length(opt$bowtie_host_output)]))
        #Clean up
        file.remove(list.files(pattern="*.bt2$"))
        file.remove("RvsHS.sam")

        #Disconsider SE leftover if it is too small
        senahs<-list.files(pattern="*NAHS_SE.fastq")
        if(length(senahs) > 0){
            senahsfs<-file.size(senahs)
            if(senahsfs < 1000000){
                file.remove(senahs)
            }
        }
        opt$nahsreads<-sort(list.files(pattern="*_NAHS_[SP]E"))

        if(opt$workdir != opt$sampledir){
            #Bank NAHS reads to project directory
            file.copy(opt$nahsreads, opt$readsdir)
        }

    }
    return(opt)
}
