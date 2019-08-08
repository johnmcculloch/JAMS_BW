#' make_colour_dictionary(variable_list=NULL, pheno=NULL)
#'
#' Returns a dictionary of colours for each class within each variable_list
#' @export

make_colour_dictionary <- function(variable_list = NULL, pheno = NULL){

    JAMSpalettenames <- c("royalblue4", "red1", "purple2", "orange2", "chocolate4", "navyblue", "green", "violet", "turquoise2", "tomato", "tan3", "steelblue3", "springgreen2", "violetred1", "slategray3", "wheat3", "violetred4",  "dodgerblue2", "black")
    JAMSpalette <- rep(JAMSpalettenames, 10)
    #Get list of variables which are discrete
    discretes <- unique(unname(unlist(variable_list[which(!(names(variable_list) %in% c("sample", "continuous")))])))
    cdict<-NULL
    cdict<-list()
    s=0
    for(d in 1:length(discretes)){
        classnames<-sort(unique(pheno[,discretes[d]]))
        classcols<-JAMSpalette[(s+1):(s+length(classnames))]
        classhex<-col2hex(classcols)
        cdict[[d]]<-data.frame(Name=classnames, Colour=classcols, Hex=classhex, stringsAsFactors = FALSE)
        names(cdict)[d]<-discretes[d]
        s=(s+length(classnames))
    }

    return(cdict)
}
