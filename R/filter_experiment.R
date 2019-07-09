#' filter_experiment(mgseqobj=NULL, featmaxatleastPPM=0, featcutoff=c(0, 0), samplesToKeep=NULL, featuresToKeep=NULL, asPA=FALSE, asPPM=TRUE, mgSeqnorm=FALSE)
#'
#' Filters a metagenomeSeq object by several criteria.
#' @export

filter_experiment <- function(mgseqobj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, asPA = FALSE, asPPM = TRUE, mgSeqnorm = FALSE){

        #Get only samples you asked for
        if(!(is.null(samplesToKeep))){
            mgseqobj=mgseqobj[ , samplesToKeep]
        }

        #Get only samples you asked for
        if(!(is.null(featuresToKeep))){
            mgseqobj=mgseqobj[featuresToKeep, ]
        }

        #If setting featmaxatleastPPM or featcutoff to anything other than the defaults, then return as PPM.
        if((featmaxatleastPPM != 0) || featcutoff != c(0, 0)){
            asPPM=TRUE
        }

        #Flush out empty rows
        mgseqobj=mgseqobj[(rowSums(MRcounts(mgseqobj)) > 0), ]
        #Flush out empty Samples
        emptysamples<-names(which(colSums(MRcounts(mgseqobj)) == 0) == TRUE)
        if(length(emptysamples) > 0){
            print(paste("Samples", paste0(emptysamples, collapse=", "), "are empty and will be discarded."))
            validsamples<-names(which(colSums(MRcounts(mgseqobj)) > 0) == TRUE)
            mgseqobj=mgseqobj[ , validsamples]
        }

        #Transform to PPM if applicable
        if(asPPM == TRUE){
            mat = MRcounts(mgseqobj, norm = mgSeqnorm, log = FALSE)
            mat=mat[!is.na(row.names(mat)),]
            #flush out empty rows
            #mat=mat[rowSums(mat)>0, ]
            #transform into PPM
            mat2<-sapply(1:ncol(mat), function(x){ mat[,x]<-((mat[,x]/(sum(mat[,x])))*1000000)} )
            colnames(mat2)<-colnames(mat)
            #Round to integer
            for(c in 1:ncol(mat2)){
                mat2[,c]<-round(mat2[,c], 0)
            }
            countmatrix<-mat2
            pheno2<-pData(mgseqobj)
            ftt<-fData(mgseqobj)
            ##Create a non-normalised class object in metagenomeseq (mgseq)
            phenotypeData = AnnotatedDataFrame(pheno2)
            ttdata = AnnotatedDataFrame(ftt)
            mgseqobj = newMRexperiment(countmatrix, phenoData=phenotypeData, featureData=ttdata)
        }

        #Discard features which do not match certain criteria
        if(!(missing(featcutoff))){
            if(length(featcutoff) != 2){
                stop("Please specify the minimum PPM in what percentage of the samples you want with a numerical vector of size two. For example, featcutoff=c(2000,10) would discard features which are not at least 2000 PPM in at least 10% of samples.")
            }
            thresholdPPM<-featcutoff[1]
            sampcutoffpct<-min(featcutoff[2], 100)

            featuresToKeep=which((lapply(1:nrow(MRcounts(mgseqobj, norm = FALSE, log = FALSE)), function(x) length(which(MRcounts(mgseqobj, norm = FALSE, log = FALSE)[x,] >= thresholdPPM))/ncol(MRcounts(mgseqobj, norm = FALSE, log = FALSE)))) > (sampcutoffpct/100))
            mgseqobj=mgseqobj[featuresToKeep, ]
            cutoffmsg<-paste("Feature must be >", thresholdPPM, "PPM in at least", sampcutoffpct, "% of samples", sep=" ")
            print(cutoffmsg)
        }

        if(!(missing(featmaxatleastPPM))){
            featuresToKeep = which(rowMax(MRcounts(mgseqobj)) >= featmaxatleastPPM)
            mgseqobj=mgseqobj[featuresToKeep, ]
            minPPMmsg<-paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep=" ")
            print(minPPMmsg)
        }

        if(asPA==TRUE){
            countmatrix<-MRcounts(mgseqobj)
            countmatrix[which(countmatrix != 0)]<-1
            pheno2<-pData(mgseqobj)
            ftt<-fData(mgseqobj)
            ##Create a non-normalised class object in metagenomeseq (mgseq)
            phenotypeData = AnnotatedDataFrame(pheno2)
            ttdata = AnnotatedDataFrame(ftt)
            mgseqobj = newMRexperiment(countmatrix, phenoData=phenotypeData, featureData=ttdata)
        }

        return(mgseqobj)
}
