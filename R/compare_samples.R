#' compare_samples(mgseqobj=NULL, mgSeqnorm=FALSE, glomby=NULL, stattype=c("variance", "binary"), samplesToKeep=NULL, featuresToKeep=NULL, binaryclass=NULL, invertbinaryorder=FALSE, subsetby=NULL, featmaxatleastPPM=0, featcutoff=c(0, 0), maxl2fc=NULL, minl2fc=NULL, genomecompleteness=NULL, list.data=NULL)
#'
#' Returns a list of dataframes with relative abundance statistics
#' @export

compare_samples<-function(mgseqobj=NULL, mgSeqnorm=FALSE, glomby=NULL, stattype=c("variance", "binary"), samplesToKeep=NULL, featuresToKeep=NULL, binaryclass=NULL, invertbinaryorder=FALSE, subsetby=NULL, featmaxatleastPPM=0, featcutoff=c(0, 0), maxl2fc=NULL, minl2fc=NULL, genomecompleteness=NULL, list.data=NULL){

    #Get appropriate object to work with
    obj<-mgseqobj

    #Exclude samples and features if specified
    if(!(is.null(samplesToKeep))){
        obj<-obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj<-obj[featuresToKeep, ]
    }

    analysis<-attr(obj, "analysis")

    if(!(is.null(glomby))){
        obj <- aggTax(obj, lvl = glomby, out = 'MRexperiment', norm = FALSE)
    }
    
    #Test for silly stuff
    if((stattype != "variance") && (is.null(binaryclass))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *binaryclass*")
    }
    
    if((analysis != "LKT") && (!(is.null(genomecompleteness)))){
        stop("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
    }

    if((!(is.null(genomecompleteness))) && (is.null(list.data))){
        stop("Genome completeness can only be obtained by supplying a list.data object.")
    }


    if (!(is.null(subsetby))){
        subset_points<-sort(unique((pData(obj)[, which(colnames(pData(obj))==subsetby)])))
    }else{
        subset_points<-"none"
    }

    svec<-NULL
    svec<-vector("list", length=1000)
    n<-1

    #subset by metadata column
    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep = which((pData(obj)[,which(colnames(pData(obj))==subsetby)])==subset_points[sp])
            print(paste("Calculating", subset_points[sp]))
            subsetname<-subset_points[sp]
        } else {
            samplesToKeep = rownames(pData(obj))
            subsetname<-"no_subsetting"
        }

        #Discard features which do not match certain criteria
        if(!(is.null(featcutoff))){
            thresholdPPM<-featcutoff[1]
            sampcutoffpct<-min(featcutoff[2], 100)
            cutoffmsg<-paste("Feature must be >", thresholdPPM, "PPM in at least ", sampcutoffpct, "% of samples", sep="")
        } else {
            cutoffmsg<-"Feature must be > 0 PPM in at least 0% of samples"
            featcutoff<-c(0,0)
        }

        if(!(is.null(featmaxatleastPPM))){
            minPPMmsg<-paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep=" ")
        } else {
            minPPMmsg<-"Highest feature must be > 0 PPM"
            featmaxatleastPPM<-0
        }

        currobj<-filter_experiment(mgseqobj=obj, featmaxatleastPPM=featmaxatleastPPM, featcutoff=featcutoff, samplesToKeep=samplesToKeep, asPA=asPA, asPPM=TRUE, mgSeqnorm=mgSeqnorm)

        #Get PPM matrix
        mat = MRcounts(currobj, norm = FALSE, log = FALSE)

        #Protect against rows with empty data
        rowsToKeep = which(rowSums(mat) > 0 & rownames(mat) !="")
        mat = mat[rowsToKeep, ]

        #Rename rows to include description if not taxonomic data
        if(analysis != "LKT"){
            feattable<-fData(obj)
            feattable$Feature<-paste(feattable$Accession, feattable$Description, sep="-")
            rownames(mat)<-feattable$Feature[match(rownames(mat), feattable$Accession)]
        }
        matrixSamples<-colnames(mat)
        matrixRows<-rownames(mat)

        #Discard taxa below required level of completeness
        if(!(is.null(genomecompleteness))){
            genomecompletenessdf<-get_genome_completeness(pheno=pData(currobj), list.data=list.data)
            featuresToKeep2 = rownames(genomecompletenessdf)[which(rowMax(as.matrix(genomecompletenessdf)) >= genomecompleteness)]
            mat<-mat[(rownames(mat)[(rownames(mat) %in% featuresToKeep2)]), ]
            completenessmsg<-paste("Genome completeness >", genomecompleteness)
        }
        topcats <- nrow(mat)

        #Calculate matrix stats and get new matrix.
        if(stattype == "variance"){
            matstats<-calculate_matrix_stats(countmatrix=mat, stattype="variance")
            matstats[]<-round(matstats[], 0)
            colnames(matstats)<-c("Standard_Deviation(PPM)","Median_Absolute_Deviation(PPM)")
            statmsg<-stattype
        } else if (stattype == "binary"){
            cl = pData(currobj)[ , which(colnames(pData(currobj))==binaryclass)]
            discretenames<-sort(unique(cl))
            #Checkpoint: classes must be binary
            if(length(discretenames)!=2){
                stop("The binaryclass argument can only admit categories with exactly two classes.")
            }
            orderrowstatby<-"pval"
            matstats<-calculate_matrix_stats(countmatrix=mat, uselog=FALSE, stattype="binary", classesvector = cl, invertbinaryorder=invertbinaryorder)
            matstats$absl2fc<-abs(matstats$l2fc)

            if(!(missing(minl2fc))){
                matstats<-subset(matstats, absl2fc > minl2fc)
                minl2fcmsg<-paste("log2foldchange >", minl2fc)
            } else {
                minl2fcmsg<-"log2foldchange > 0"
            }

            if(!(missing(maxl2fc))){
                matstats<-subset(matstats, absl2fc < maxl2fc)
                maxl2fcmsg<-paste("log2foldchange <", maxl2fc)
            } else {
                maxl2fcmsg<-"log2foldchange < Inf"
            }
            statmsg<-paste(stattype, binaryclass, sep="_")
        }
        svec[[n]]<-as.data.frame(matstats)
        if(!(is.null(glomby))){
            analysisname<-glomby
        } else {
            analysisname<-analysis
        }
        stattitle<-paste(analysisname, statmsg, subsetname, sep="_")
        names(svec)[n]<-stattitle
        n <- n + 1
    }
    #Redefine graphics list as ones only containing plots
    svec<-svec[sapply(svec, function(x){!(is.null(x))})]

    return(svec)
}
