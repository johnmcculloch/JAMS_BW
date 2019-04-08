#' plot_feature_relabund<-function(mgseqobj=NULL, mgSeqnorm=FALSE, feature=NULL, glomby=NULL, accession=NULL, groupby=NULL, colourby=NULL, shapeby=NULL, subsetby=NULL, uselog=TRUE, statmeth="wilcox.test", samplesToKeep=NULL, featuresToKeep=NULL, signiflabel=c("p.signif", "p.format"), plottitle=NULL, asPA=FALSE, subsetbytaxlevel=NULL, taxtable=NULL, list.data=NULL, cdict=NULL, ...)
#'
#' Plots the relative abundance of an entity present in a sample or group of samples and subsets by taxonomy.
#' @export

plot_feature_relabund<-function(mgseqobj=NULL, mgSeqnorm=FALSE, feature=NULL, glomby=NULL, accession=NULL, groupby=NULL, colourby=NULL, shapeby=NULL, subsetby=NULL, uselog=TRUE, statmeth="wilcox.test", samplesToKeep=NULL, featuresToKeep=NULL, signiflabel="p.format", plottitle=NULL, asPA=FALSE, subsetbytaxlevel=NULL, taxtable=NULL, list.data=NULL, cdict=NULL,max_categories = 3, ...){

    require(reshape2)
    obj<-mgseqobj

    #Exclude samples and features if specified
    if(!(is.null(samplesToKeep))){
        obj<-obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj<-obj[featuresToKeep, ]
    }

    analysis<-attr(obj, "analysis")

    if(!(missing(glomby))){
        obj <- aggTax(obj, lvl = glomby, out = 'MRexperiment', norm = FALSE)
    }

    if((analysis == "LKT") && (!(is.null(subsetbytaxlevel)))){
        stop("Plot taxa is only for functional (not taxonomic) analyses. It will additionally plot which taxa bear which of the funcitons of interest. If you are trying to plot only a certain taxon itself, then use the feature argument with a taxonomical metagenomeSeq experiment.")
    }

    #Set asPA as FALSE if it is missing
    if(missing(asPA)){
        asPA=FALSE
    }

    if (!(missing(subsetby))){    
        subset_points<-sort(unique((pData(obj)[, which(colnames(pData(obj))==subsetby)])))
    }else{
        subset_points<-"none"
    }

    if (is.null(feature) && is.null(accession)){
        stop("You must specify either the feature name (like a taxon) or an accession number for functional analyses.")
    }

    if (!is.null(accession)){
        feature<-accession
    }

    if(is.null(plottitle)){
        featdesc<-paste(feature, fData(obj)[feature,"Description"], sep="-")
        if(length(feature) < 3){
            plottitle<-paste(featdesc, collapse = ' AND ')
        } else {
            plottitle<-paste("Several features of", analysis)
        }
    }

    gvec<-NULL
    gvec<-vector("list",length=1000)
    plotcount=1

    for (sp in 1:length(subset_points)){

        if (!(is.null(subsetby))){ 
            samplesToKeep = which((pData(obj)[,which(colnames(pData(obj))==subsetby)])==subset_points[sp])
            currobj=obj[ , samplesToKeep]
            maintit<-paste(plottitle, paste("Subset:", subset_points[sp]), sep="\n")
        } else {
            currobj = obj
            maintit <- plottitle
        }

        currobj<-filter_experiment(mgseqobj=currobj, asPA=asPA, asPPM=TRUE, mgSeqnorm=mgSeqnorm)

        if (is.factor(pData(currobj)[,groupby])) { # use order if factor
          discretenames <- levels(pData(currobj)[,groupby])
        } else {
          discretenames<-sort(unique(as.character(pData(currobj)[,which(colnames(pData(currobj))==groupby)])))
        }
        classIndex<-NULL
        classIndex<-list()
        for (n in 1:length(discretenames)){
            classIndex[[n]]=which(pData(currobj)[,which(colnames(pData(currobj))==groupby)]==discretenames[n])
            names(classIndex)[[n]] <- discretenames[n]
        }
        par(cex.axis=0.5, cex.main=1, cex.sub=0.5, las=2)
        
        mat = MRcounts(currobj, norm = FALSE, log = uselog)

        #Protect against feature being absent
        if(all(feature %in% rownames(mat))){

            #if there is more than a single feature, add their relabunds up in an appropriate manner. 
            if(length(feature) == 1){
                l = lapply(classIndex, function(j) { mat[feature, j] })
            } else {
                mat2<-mat[feature, ]
                if(uselog==TRUE){
                    mat3<-(2^(mat2)-1)
                    mat2<-t(colSums(mat3))
                    mat2<-log2(mat2+1)
                } else {
                    mat2<-t(colSums(mat2))
                }
                l = lapply(classIndex, function(j) { mat2[, j] })
            }
            
            y = unlist(l)
            x = rep(seq(along = l), sapply(l, length))
            xnam = rep(names(l), sapply(l, length))
            dat<-data.frame(xnam=xnam, x=x, y=y)
            
            if(length(discretenames) < nrow(pData(currobj))){
                jitfact=-(0.3/nrow(pData(currobj)))*(length(discretenames)) + 0.25
            } else {
                jitfact=0
            }
            
            if (!(is.null(shapeby))){
                shapeindex = lapply(classIndex, function(j) { pData(currobj)[j, which(colnames(pData(currobj))==shapeby)] })
                shp=unlist(shapeindex)
                dat$shape<-shp
            }
            
            if(!(is.null(colourby))){
                colourindex = lapply(classIndex, function(j) { pData(currobj)[j, which(colnames(pData(currobj))==colourby)] })
                col=unlist(colourindex)
                dat$colours<-col
            }

            #I know, I know, call me inelegant, but it works.
            if(!(is.null(colourby))){
                p <- ggplot(dat, aes(pData(currobj)[,groupby], y, colour = colours))
                if(is.numeric(dat$colours)){
                    p <- p + scale_color_gradient(low="blue", high="red")
                } else {
                    #if there is a colour dictionary, then use that
                    if(!(is.null(cdict))){
                        ct<-cdict[[colourby]]
                        groupcols<-setNames(as.character(ct$Colour), as.character(ct$Name))
                        p <- p + scale_color_manual(values = groupcols)
                    }
                }
                
            } else {
                p <- ggplot(dat, aes(factor(xnam), y))
            }
            
            if(asPA!=TRUE){
                p <- p + geom_boxplot(outlier.shape=NA)
            }
            
            if (!(is.null(shapeby))){
                p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape=shape))
            } else {
                p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
            }
            p <- p + theme_minimal()
            
            if((length(discretenames) > 1) && (length(discretenames) <= max_categories)){
                if(missing(signiflabel)){
                    signiflabel="p.format"
                }
                #Add pval
                my_comparisons<-combn(discretenames, m=2, simplify=FALSE)
                p <- p + stat_compare_means(method = statmeth, comparisons = my_comparisons, paired=FALSE, label=signiflabel)
            } else {
                print("There are too many combinations to plot significance.")
            }

            p <- p + ggtitle(maintit)

            if(!(is.null(colourby))){ 
                p <- p + labs(colour = colourby)
            }

            if(!(is.null(shapeby))){ 
                p <- p + labs(shape = shapeby)
            }
            
            if(uselog==TRUE){
                ytit<-"log2 of (Relative Abundance in PPM)"
                #p <- p + scale_y_continuous(expand = c(0, 0),  sec.axis = sec_axis(trans=~expm1(.*log(2)), name="Percentage relative abundance"))
            } else {
                if(asPA==TRUE){
                    ytit<-"Absence/Presence"
                    p<- p + scale_y_discrete(limits=c(0,1))
                } else {
                    ytit<-"Relative Abundance in PPM"
                }
            }
            p <- p + labs(x=groupby, y=ytit)
            p <- p + theme(axis.text.x=element_text(angle=90, size = rel(1), colour="black"))
            p <- p + theme(plot.title = element_text(size=10))

            gvec[[plotcount]]<-p
            plotcount = plotcount + 1

            if(!(is.null(subsetbytaxlevel))){
                #split production of that feature or features by taxon.
                #Amass the counts for that taxonomy
                print(paste("Subsetting counts by", subsetbytaxlevel))
                #Make a colour dictionary by phylum
                tt<-taxtable
                tt[] <- lapply(tt, as.character)
                tt$Colour<-rainbow(length(tt$Phylum), start=.5, end=.1, alpha=1)[rank(tt$Phylum)]
                ctit<-"black"
                tt$Colour[which(tt$Phylum=="Unclassified")]<-"#000000"

                #Get data for features 
                wantedSamples<-rownames(pData(currobj))
                featureobjects<-paste(wantedSamples, "featuredose", sep="_")
                featuredoses<-list.data[featureobjects]
                names(featuredoses)<-wantedSamples

                #subset doses to contain only the analysis wanted
                analysisdoses<-NULL
                analysisdoses<-list()
                for(a in 1:length(featuredoses)){
                    totalCDSbases<-as.numeric(featuredoses[[a]][(which(featuredoses[[a]]$Accession == "CDS")), (which(colnames(featuredoses[[a]]) == "NumBases"))])
                    LKTindf<-grep("LKT__", colnames(featuredoses[[a]]))
                    doses<-subset(featuredoses[[a]], Analysis == analysis)
                    doses<-subset(doses, Accession %in% feature)
                    doses<-as.data.frame(doses[ , LKTindf])
                    counts<-colSums2(as.matrix(doses))
                    #Transform to PPM
                    counts<-round((counts/totalCDSbases)*1000000, 0)
                    names(counts)<-names(doses)
                    counts<-counts[which(counts != 0)]
                    counts<-as.data.frame(counts)
                    colnames(counts)[1]<-"PPM"
                    counts$Taxon<-rownames(counts)
                    analysisdoses[[a]]<-counts
                    names(analysisdoses)[a]<-names(featuredoses[a])
                }

                featureall<-bind_rows(analysisdoses, .id = "id")
                featureall[is.na(featureall)]<-0
                colnames(featureall)[1]<-"Sample"

                submat <- featureall %>% tidyr::spread(Sample, PPM)
                submat[is.na(submat)]<-0
                rownames(submat)<-submat$Taxon
                submat$Taxon<-NULL
                submat<-as.matrix(submat)
                
                #Add missing data for samples with 0
                samplesnotzero<-wantedSamples[(wantedSamples %in% colnames(submat))]
                sampleszero<-wantedSamples[!(wantedSamples %in% samplesnotzero)]

                #If sample is zero, add a column for that sample with zero in order to get the medians right.
                if(length(sampleszero) > 0){
                    extradata<-matrix(nrow=nrow(submat), ncol=length(sampleszero), data=0)
                    rownames(extradata)<-rownames(submat)
                    colnames(extradata)<-sampleszero
                    submat<-cbind(submat, extradata)
                }

                #Split by discrete group
                classIndex2<-NULL
                classIndex2<-list()
                for (n in 1:length(discretenames)){
                    classIndex2[[n]]=rownames(pData(currobj))[which(pData(currobj)[,which(colnames(pData(currobj))==groupby)]==discretenames[n])]
                    names(classIndex2)[[n]] <- discretenames[n]
                }

                for(ds in 1:length(discretenames)){
                    #Get a matrix with only the samples applicable
                    samplesingroup<-classIndex2[[discretenames[ds]]]
 
                    #Fix dimensions if a single feature
                    if(nrow(submat) != 1){
                        groupstratmat<-submat[ , samplesingroup]
                    } else {
                        groupstratmat<-t(submat[ , samplesingroup])
                        rownames(groupstratmat)<-rownames(submat)
                        colnames(groupstratmat)<-samplesingroup
                    }

                    #Eliminate empty rows, if there is anything to begin with.
                    skipplot = FALSE
                    if(nrow(groupstratmat) > 1){
                        notzero<-names(which(rowSums(groupstratmat) > 0))
                        if(length(notzero) > 1 ){                        
                            groupstratmat<-groupstratmat[notzero, ]
                        } else if (length(notzero) == 1)  {
                            groupstratmat<-t(groupstratmat[notzero, ])
                            rownames(groupstratmat)<-notzero
                        } else {
                            #nothing to see here
                            skipplot = TRUE
                        }
                    } else if (nrow(groupstratmat) == 0){
                        stop("No Rows")
                    }

                    #Plot if there is enough data to do so
                    #Order by median
                    if(skipplot == FALSE){
                        Meds<-rowMedians(groupstratmat)
                        groupstratmat<-as.data.frame(groupstratmat)
                        groupstratmat$Medians<-Meds
                        groupstratmat<-groupstratmat[order(groupstratmat$Medians, decreasing=TRUE), ]
                        #Show only top 35, or else will not fit in plot.
                        groupstratmat<-groupstratmat[1:min(nrow(groupstratmat), 35), ]
                        groupstratmat$Medians <- NULL
                        groupstratmat$Taxon <- rownames(groupstratmat)
                        groupstratmat$Taxon <- factor(groupstratmat$Taxon, levels=groupstratmat$Taxon)

                        mdata2 <- melt(groupstratmat, id="Taxon")

                        if(uselog==TRUE){
                            mdata2$value<-log2((mdata2$value)+1)
                        }

                        if(subsetbytaxlevel %in% c("Class", "Order", "Family", "Genus", "LKT")){
                            mdata2$Colour <- tt$Colour[match(mdata2$Taxon, tt[, subsetbytaxlevel])]
                            mdata2$PhylumInfo<-tt[ , "Phylum"][match(mdata2$Taxon, tt[ , subsetbytaxlevel])]
                            PhyColours<-unique(mdata2$Colour)
                            names(PhyColours)<-mdata2$PhylumInfo[match(PhyColours, mdata2$Colour)]
                            p <- ggplot(mdata2, aes(factor(Taxon), value, fill=PhylumInfo))
                            p <- p + geom_boxplot(outlier.shape=NA)
                            p <- p + scale_fill_manual(values = PhyColours)
                        } else {
                            p <- ggplot(mdata2, aes(factor(Taxon), value))
                            p <- p + geom_boxplot(outlier.shape=NA)
                        }

                        p <- p + theme_minimal()

                        subplottitle<-paste(analysis, plottitle, paste("Only within group", discretenames[ds], sep="= "), sep="\n")
                        p <- p + ggtitle(subplottitle)

                        if(uselog==TRUE){
                            yexp<-"log2 PPM in Sample"
                        } else {
                            yexp<-"PPM in Sample"
                        }

                        p <- p + labs(x=analysis, y=yexp)
                        p <- p + theme(axis.text.x=element_text(angle=90, size = rel(0.9), colour="black"))
                    } else {
                        #Give message that there is nothing in this group
                        plot.new()
                        grid.table(c(feature, "not found in group", discretenames[ds]), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                        p <- recordPlot()
                    }
                    gvec[[plotcount]]<-p
                    plotcount = plotcount + 1
                }
            }
        }
    }

    #Redefine graphics list as ones only containing plots
    gvec<-gvec[sapply(gvec, function(x){!(is.null(x))})]
    
    return(gvec)
}
