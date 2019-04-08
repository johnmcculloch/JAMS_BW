#' plot_alpha_diversity(mgseqobj=NULL, glomby=NULL,  groupby=NULL, colourby=NULL, shapeby=NULL, subsetby=NULL, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"), ignoresamples=NULL, signiflabel=c("p.signif", "p.format"), plottitle=NULL, cdict=NULL, ...)
#'
#' Plots the alpha diversity in a sample or group of samples.
#' @export

plot_alpha_diversity<-function(mgseqobj=NULL, glomby=NULL,  groupby=NULL, colourby=NULL, shapeby=NULL, subsetby=NULL, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"), statmeth="wilcox.test", samplesToKeep=NULL, featuresToKeep=NULL, signiflabel="p.format", plottitle=NULL, cdict=NULL, max_categories = 3,...){

    require(vegan)

    #Get appropriate object to work with
    obj<-mgseqobj
    
    #Exclude samples and features if specified
    if(!(is.null(samplesToKeep))){
        obj<-obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj<-obj[featuresToKeep, ]
    }

    #Define analysis type
    analysis<-attr(obj, "analysis")
    analysisname<-analysis

    if(!(is.null(glomby))){
        obj <- aggTax(obj, lvl = glomby, out = 'MRexperiment', norm = FALSE)
        analysisname<-glomby
    }

    if (!(is.null(subsetby))){
        subset_points<-sort(unique((pData(obj)[, which(colnames(pData(obj))==subsetby)])))
    } else {
        subset_points<-"none"
    }

    #Create list vector to hold plots
    gvec<-NULL
    gvec<-vector("list", length=1000)
    gvn<-1

    maintit<-paste("Alpha diversity of", analysisname, "by", groupby)

    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){ 
            samplesToKeep = which((pData(obj)[,which(colnames(pData(obj))==subsetby)])==subset_points[sp])
            currobj=obj[ , samplesToKeep]
            subsetmaintit<-paste(maintit, paste("Subset:", subset_points[sp]), sep="\n")
        }else{
            currobj=obj
            subsetmaintit<-maintit
        }

        if (is.factor(pData(currobj)[,groupby])) { # use order if factor
          groupbynames <- levels(pData(currobj)[,groupby])
        } else {
          groupbynames<-sort(unique(as.character(pData(currobj)[,which(colnames(pData(currobj))==groupby)])))
        }
        classIndex<-NULL
        classIndex<-list()
        for (n in 1:length(groupbynames)){
            classIndex[[n]]=which(pData(currobj)[,which(colnames(pData(currobj))==groupby)]==groupbynames[n])
            names(classIndex)[[n]] <- groupbynames[n]
        }
        par(cex.axis=0.5, cex.main=1, cex.sub=0.5, las=2)

        mat = MRcounts(currobj, norm = FALSE, log = FALSE)
        #calculate alpha diversity measures 
        tmat<-t(mat)
        alphadiv<-estimateR(tmat)
        for(meas in c("invsimpson", "simpson", "shannon")){
            alphadiv2<-diversity(tmat, index=meas)
            alphadiv<-rbind(alphadiv, alphadiv2[colnames(alphadiv)])
            rownames(alphadiv)[nrow(alphadiv)]<-meas
        }
        alphadiv[c(3,5,6,7,8), ]<-round(alphadiv[c(3,5,6,7,8), ], 2)
        alphadiv[c(1,2,4), ]<-round(alphadiv[c(1,2,4), ], 0)

        #Loop by measure
        for(measurename in measures){
            wantedmeasure<-switch(measurename, "Observed"="S.obs", "Chao1"="S.chao1", "ACE"="S.ACE", "Shannon"="shannon", "Simpson"="simpson", "InvSimpson"="invsimpson")
        
            l = lapply(classIndex, function(j) { alphadiv[wantedmeasure, j] })
            y = unlist(l)
            x = rep(seq(along = l), sapply(l, length))
            xnam = rep(names(l), sapply(l, length))
            dat<-data.frame(xnam=xnam, x=x, y=y)

            if(length(groupbynames) < nrow(pData(currobj))){
                jitfact=-(0.3/nrow(pData(currobj)))*(length(groupbynames)) + 0.25
            } else {
                jitfact=0
            }

            if(!(is.null(shapeby))){
                shapeindex = lapply(classIndex, function(j) { pData(currobj)[j, which(colnames(pData(currobj))==shapeby)] })
                shp=unlist(shapeindex)
                dat$shape<-shp
            }           

            if(!(is.null(colourby))){
                colourindex = lapply(classIndex, function(j) { pData(currobj)[j, which(colnames(pData(currobj))==colourby)] })
                col=unlist(colourindex)
                dat$colours<-col
            }

            #Make basic violin plot
            #p <- ggplot(dat, aes(pData(currobj)[,groupby], y))
            p <- ggplot(dat, aes(xnam, y))
            p <- p + geom_violin() + geom_boxplot(width=.1, outlier.shape=NA)

            #I know, I know, call me inelegant, but it works.
            if(!(is.null(colourby))){
                if(!(is.null(shapeby))){                
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0),  aes(shape=shape, colour=colours))
                } else {
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0),  aes(colour=colours))
                }
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
                if(!(is.null(shapeby))){
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0),  aes(shape=shape))
                } else {
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                }
            }

            p <- p + theme_minimal()

            if((length(groupbynames) > 1) && (length(groupbynames) <= max_categories)){
                if(is.null(signiflabel)){
                    signiflabel="p.format"
                }
                #Add pval
                my_comparisons<-combn(groupbynames, m=2, simplify=FALSE)
                p <- p + stat_compare_means(method = statmeth, comparisons = my_comparisons, paired=FALSE, label=signiflabel)
            } else {
                print("There are too many combinations to plot significance.")
            }

            measureexpl<-switch(measurename, "Observed"="Number of observed features", "Chao1"="Chao1 index", "ACE"="Abundance Based Coverage Estimator", "Shannon"="Shannon Index", "Simpson"="Simpson Index", "InvSimpson"="Inverse Simpson Index")

            plottit<-paste(subsetmaintit, paste("Measure:", measureexpl), sep="\n")
            p <- p + ggtitle(plottit)

            if(!(is.null(colourby))){ 
                p <- p + labs(colour = colourby)
            }

            if(!(is.null(shapeby))){ 
                p <- p + labs(shape = shapeby)
            }

            p <- p + labs(x=groupby, y=measureexpl)
            p <- p + theme(axis.text.x=element_text(angle=90, size = rel(1), colour="black"))
            p <- p + theme(plot.title = element_text(size=10))

            gvec[[gvn]]<-p
            gvn <- gvn + 1
        }
    }
    #Redefine graphics list as ones only containing plots
    gvec<-gvec[sapply(gvec, function(x){!(is.null(x))})]

    return(gvec)
}
