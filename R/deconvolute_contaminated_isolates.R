#' deconvolute_contaminated_isolates
#' This function is experimental.
#' Function to deconvolute contaminated isolates.
#'
#' JAMSalpha function
#' @export

deconvolute_contaminated_isolates<-function(opt=NULL, contigsdata=NULL, assemblystats=NULL, outputdir=NULL){

    #Get data from somewhere
    if(is.null(contigsdata)){
        contigsdata<-opt$contigsdata
    }
    #Get data from somewhere
    if(is.null(assemblystats)){
        assemblystats<-opt$assemblystats
    }

    taxlvls<-c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %in% unique(assemblystats$TaxLevel)]
    taxtags<-c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")[c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %in% unique(assemblystats$TaxLevel)]

    #Evaluate the current state of affairs
    LKTassembly<-subset(assemblystats, TaxLevel == "LKT")
    totalLKTgenomes<-sum(LKTassembly$ProbNumGenomes)
    if(totalLKTgenomes < 1.25){
        print("This sample does not look contaminated with more than a single known species.")
        wantedassemblies<-NULL
    } else {
        print(paste("Total genome size at the LKT level is currently at", (totalLKTgenomes * 100), " %. This suggests contamination with more than a single species. Will try and deconvolute species as best as possible. Keep in mind this algorithm is experimental and subject to change."))
        
        #Eliminate all assembly stats below 1000 bp as this is confounding non informative.
        assemblystatsdenoised<-subset(assemblystats, ContigSum > 1000)

        #How close is genome completeness to an integer (i.e. multiple of a full genome)?
        closesnesstoint<-function(x){
            closeness<-abs((round(x, 0))-(x))

            return(closeness)
        }
        assemblystatsdenoised$ClosenessToIntegrity<-lapply(1:nrow(assemblystatsdenoised), function (x) { closesnesstoint(assemblystatsdenoised$ProbNumGenomes[x] )} )

        #Find out which is the Last Common ancestor. This is the level at which there is a single taxon which is not unclassified.
        num_class_taxa<-function(taxlevel){
            asdf<-subset(assemblystatsdenoised, TaxLevel == taxlevel)
            dunno<-paste0(taxtags[which(taxlvls==taxlevel)], "Unclassified")
            asdf<-subset(asdf, Taxon != dunno)
            numtax<-length(unique(asdf$Taxon))
            return(numtax)
        }
        richness_levels<-sapply(taxlvls, function (x) { num_class_taxa(x) } )        
        highest_level_with_info<-names(richness_levels[which(richness_levels > 1)])[1]
        LCAtl <- names(richness_levels)[(which(names(richness_levels) == highest_level_with_info) - 1)]
        asdf<-subset(assemblystatsdenoised, TaxLevel == LCAtl)
        dunno<-paste0(taxtags[which(taxlvls==LCAtl)], "Unclassified")
        asdf<-subset(asdf, Taxon != dunno)
        LCA<-as.character(asdf$Taxon)

        print(paste("The Last Common Ancestor of the taxa in the sample is", LCAtl, LCA))
        
        assemblystatsofinterest<-subset(assemblystatsdenoised, TaxLevel == highest_level_with_info)
        #remove unclassified, as these contigs cannot be attributed to either component.
        assemblystatsofinterest<-subset(assemblystatsofinterest, Taxon != paste0(taxtags[which(taxlvls==highest_level_with_info)], "Unclassified"))
        #Remove small N90s
        assemblystatsofinterest<-subset(assemblystatsofinterest, N90 > 2000)
        
        wantedassemblies<-subset(assemblystatsofinterest, ProbNumGenomes > 0.9)
        wantedassemblies<-subset(wantedassemblies, ProbNumGenomes < 1.25)
    }
    opt$wantedassemblies<-wantedassemblies

    return(opt)
}
