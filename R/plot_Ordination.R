#' plot_Ordination(mgseqobj=NULL, glomby=NULL, subsetby=NULL, samplesToKeep=NULL, featuresToKeep=NULL, mgSeqnorm=FALSE, featmaxatleastPPM=0, featcutoff=c(0, 0), algorithm="PCA", colourby=NULL, shapeby=NULL, log2tran=TRUE, transp=TRUE, perplx=NULL, permanova=FALSE, ellipse=c(TRUE, FALSE, "auto"), plottit=NULL, plot3D=FALSE, theta=130, phi=60, cdict=NULL, grid = TRUE, forceaspectratio=NULL, ...)
#'
#' Creates ordination plots based on tSNE or PCA
#' @export

plot_Ordination<-function(mgseqobj=NULL, glomby=NULL, subsetby=NULL, samplesToKeep=NULL, featuresToKeep=NULL, mgSeqnorm=FALSE, featmaxatleastPPM=0, featcutoff=c(0, 0), algorithm="PCA", colourby=NULL, shapeby=NULL, sizeby=NULL, log2tran=TRUE, transp=TRUE, perplx=NULL, permanova=FALSE, ellipse=FALSE, plottit=NULL, plot3D=FALSE, theta=130, phi=60, cdict=NULL, grid = TRUE, forceaspectratio=NULL, ...){
    #Get appropriate object to work with
    obj<-mgseqobj
    
    #Exclude samples and features if specified
    if(!(is.null(samplesToKeep))){
        obj<-obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj<-obj[featuresToKeep, ]
    }

    #Aggregate if required
    if(!(is.null(glomby))){
        obj <- aggTax(obj, lvl = glomby, out = 'MRexperiment', norm = FALSE)
    }

    #Define analysis type
    analysis<-attr(obj, "analysis")

    if (!(is.null(subsetby))){
        subset_points<-sort(unique((pData(obj)[, which(colnames(pData(obj))==subsetby)])))
    } else {
        subset_points<-"none"
    }

    #Create list vector to hold plots
    gvec<-NULL
    gvec<-vector("list", length=length(subset_points))
    pcatitbase<-paste(algorithm, "of", analysis)

    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep = rownames(pData(obj))[which((pData(obj)[,which(colnames(pData(obj))==subsetby)])==subset_points[sp])]
            pcatit<-paste(pcatitbase, "within", subset_points[sp])
            currobj=obj[,samplesToKeep]
        } else {
            currobj=obj
            samplesToKeep = rownames(pData(currobj))
            pcatit<-pcatitbase
        }

        currobj<-filter_experiment(mgseqobj=currobj, featmaxatleastPPM=featmaxatleastPPM, featcutoff=featcutoff, samplesToKeep=samplesToKeep, asPA=FALSE, asPPM=TRUE, mgSeqnorm=mgSeqnorm)

        mat = MRcounts(currobj)

        #log2 transform if applicable
        if(log2tran == TRUE){
            mat2<-sapply(1:ncol(mat), function(x){ mat[,x]<-(log2(mat[,x] + 1))} )
            colnames(mat2)<-colnames(mat)
            mat<-mat2
        }

        n = min(nrow(mat), 1000)
        comp = 1:3
        rowsToKeep <- which(rowSums(mat) > 0)
        rowVars <- rowSds(mat[rowsToKeep, ])
        rowIndices <- rowsToKeep[order(rowVars, decreasing = TRUE)[seq_len(n)]]
        mat <- mat[rowIndices, ]
        if (transp == TRUE) {
            mat = t(mat)
        }

        if (algorithm == "tSNE"){
            #tSNE algorithm
            permanova=FALSE
            if(is.null(perplx)){
                perplx<-round(nrow(pData(currobj))*0.3, 0)
            }
            set.seed(perplx)
            tsne_out <- Rtsne(mat, dims = 3, initial_dims = 500, perplexity = perplx, theta = 0.5, check_duplicates = TRUE, pca = TRUE, max_iter = 1000)
            dford<-as.data.frame(tsne_out$Y)
            rownames(dford)<-rownames(pData(currobj))
            colnames(dford)[1:3]<-c("PC1", "PC2", "PC3")
            xl<-"tSNE 1"
            yl<-"tSNE 2"
            zl<-"tSNE 3"
        } else {
            #Not tSNE, so use PCA
            distfun = stats::dist
            d <- distfun(mat, method = "euclidian")
            pcaRes <- prcomp(d)
            ord <- pcaRes$x
            vars <- pcaRes$sdev^2
            vars <- round(vars/sum(vars), 5) * 100

            pnmetadata<-pData(currobj)
            cats<-pnmetadata[,colourby]

            if(!(is.numeric(cats))){
                permanovap <- vegan::adonis(as.formula(paste("d ~ ", colourby)), data = pData(currobj))$aov.tab$`Pr(>F)`[1]
            } else {
                #print("Impossible to get permanova because colourby is continuous")
                permanova <- FALSE
                permanovap <- 1
            }

            xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
            yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])
            zl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[3]], vars[comp[3]])
            dford<-as.data.frame(ord[,comp])
        }
        #Add colour, size, shape
        dford$Colours<-pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj))==colourby)]
        dford$Size<-pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj))==sizeby)]
        dford$Shape<-pData(currobj)[match(rownames(dford), rownames(pData(currobj))), which(colnames(pData(currobj))==shapeby)]
        
       
        if (plot3D==TRUE) {
          aesthetic <- aes(x=PC1, y=PC2, z=PC3)
        }
        else {
          aesthetic <- aes(x=PC1, y=PC2)
        }
        p<-ggplot(dford, aesthetic)

        if(!is.null(colourby)) {
          p <- p + aes(col=Colours)
        }
  
        if (!(is.null(shapeby))){
          p <- p + aes(shape=Shape)
          numshapes<-length(unique(dford$Shape))
           p <- p + scale_shape_manual(values = 0:numshapes)
        }
        
        if (!(is.null(sizeby))){
          p <- p + aes(size=Size)
          numsizes<-length(unique(dford$Size))
          p <- p + scale_shape_manual(values = 0:numsizes)
        }

        if(is.numeric(dford$Colours)){
            #Check if there is enough variance in the continuous data to plot a gradient
            if((max(dford$Colours)- min(dford$Colours)) > 0){
                p <- p + scale_color_gradient(low="blue", high="red")
            }
        } else {
            #if there is a colour dictionary, then use that
            if(!(is.null(cdict))){
                ct<-cdict[[colourby]]
                groupcols<-setNames(as.character(ct$Colour), as.character(ct$Name))
                p <- p + scale_color_manual(values = groupcols)
            }
        }
        if(!(missing(plottit))){
            pcatit<-plottit
        }
        if ( (permanova != FALSE) && (algorithm != "tSNE") ) {
            pcatit <- paste(pcatit, paste("p <", permanovap))
        }
        if(mgSeqnorm == TRUE){
            pcatit <- paste(pcatit, (paste0("MetagenomeSeq normalization = ", as.character(mgSeqnorm))), sep="\n")
        }
        p <- p + ggtitle(pcatit) 
        p <- p + labs(colour = colourby)
        if (!(is.null(shapeby))){
            p <- p + labs(shape = shapeby)
        }

        if ( (!isFALSE(ellipse)) ) {
            if(ellipse == "auto" && algorithm != "tSNE"){
                if(permanovap < 0.05){
                    p <- p + stat_ellipse(type="norm")
                }
            } else if (ellipse == TRUE) {
                p <- p + stat_ellipse(type="norm")
            }
        }

        p <- p + theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=10))
        if (plot3D==TRUE) {
            p <- p + axes_3D(theta=theta, phi=phi) + stat_3D()
            p <- p + labs_3D(labs=c(xl,yl,zl), hjust=c(-0.25, 1,0.25)) + labs(x="", y="")
            p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) 
        } else {
            p <- p + geom_point() + labs(x=xl, y=yl) 
            if(!(is.null(forceaspectratio))){
                p <- p + theme(aspect.ratio=(1/forceaspectratio))
            }
        }

        if (grid == FALSE ){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))
        }


        gvec[[sp]]<-p
    }
    
    return(gvec)
}
