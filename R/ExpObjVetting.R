#' ExpObjVetting(ExpObj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, glomby = NULL, variables_to_fix = NULL, only_allow_CSBs = FALSE, class_to_ignore = NULL, featuresToKeep_warn_frac = 0.2, featuresToKeep_error_frac = 0.5, featuresToKeep_size_floor = 10)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions.
#'
#' featuresToKeep NAMESPACE CONTRACT:
#'   featuresToKeep is always interpreted in the POST-agglomeration namespace. That is,
#'   if you pass glomby = "Family", featuresToKeep must contain family-level names (e.g.
#'   "f__Enterobacteriaceae"), NOT the underlying LKT/CGB names that compose those families.
#'   If glomby is NULL, featuresToKeep are matched against the object's native feature names.
#'
#'   How absence is handled (see the three threshold arguments below): because we cannot, from
#'   the object alone, tell "a legitimate sparse subset" from "a wrong-namespace vector that
#'   happens to contain real taxon tokens", we DO NOT consult any external taxonomy (that would
#'   be both expensive and unsound, since a local JAMStaxtable may not match the Kraken database
#'   that built these .jams files). Instead we judge by how many supplied names are actually
#'   present in THIS object, scaled by how many were asked for:
#'     - 0 present                                   -> ERROR (nothing to plot; clearest signal of a bad vector)
#'     - small vector (<= size_floor), >=1 present   -> never error on absence; warn per missing name
#'     - large vector, absent fraction >= error_frac -> ERROR (almost certainly a namespace mistake / stale list)
#'     - large vector, absent fraction >= warn_frac  -> strong WARNING, but still proceeds
#'     - otherwise                                   -> proceed (quietly note any few absentees)
#'   The goal is explicitly NOT to prevent every misuse (the user must know roughly what they
#'   are asking for), but to guarantee we never SILENTLY dump a wrong/degenerate heatmap into a PDF.
#'
#' @param featuresToKeep_warn_frac Fraction of absent requested features above which a strong warning is issued (large vectors only). Default 0.2.
#' @param featuresToKeep_error_frac Fraction of absent requested features above which this stops with an error (large vectors only). Default 0.5.
#' @param featuresToKeep_size_floor Vector-length floor below which the error rule is suspended, so small, legitimate requests (e.g. 4 species, 2 absent) only ever warn. Default 10.
#' @export

ExpObjVetting <- function(ExpObj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, glomby = NULL, variables_to_fix = NULL, only_allow_CSBs = FALSE, class_to_ignore = NULL, featuresToKeep_warn_frac = 0.2, featuresToKeep_error_frac = 0.5, featuresToKeep_size_floor = 10){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input.")
        }

        if (only_allow_CSBs){
            obj <- filter_experiment(SEobj = ExpObj, only_allow_CSBs = only_allow_CSBs, give_info = TRUE)
        } else {
            obj <- ExpObj
        }

        if (!(is.null(glomby))){
            obj <- agglomerate_features(ExpObj = obj, glomby = glomby)
        }

        #Exclude samples and features if specified
        if (!(is.null(samplesToKeep))){
            samplesToKeep <- unique(samplesToKeep)
            samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(obj)]
            obj <- obj[, samplesToKeep]
        }

        if (!(is.null(featuresToKeep))){
            featuresToKeep <- unique(featuresToKeep)

            #### featuresToKeep presence check (object-intrinsic only) ####
            #featuresToKeep is interpreted in the POST-glom namespace (see contract above).
            #We judge purely by presence in THIS object, scaled by how many were requested.
            #We deliberately do NOT validate against any external taxonomy: it would be costly
            #(JAMStaxtable is large and this runs on every plot, often in loops) and unsound
            #(a local taxonomy may not match the Kraken database that built these .jams files,
            #so a genuinely-correct name could be wrongly rejected). The sole aim here is to
            #never silently emit a wrong/degenerate plot - not to catch every possible misuse.

            present_features <- rownames(obj)
            matched  <- featuresToKeep[featuresToKeep %in% present_features]
            absent   <- featuresToKeep[!(featuresToKeep %in% present_features)]

            n_req     <- length(featuresToKeep)
            n_match   <- length(matched)
            n_absent  <- length(absent)
            absent_frac <- n_absent / n_req

            #A helpful namespace hint for the louder messages.
            namespace_hint <- function(){
                if (!is.null(glomby)){
                    paste0("featuresToKeep is interpreted in the POST-agglomeration namespace. ",
                           "Because glomby = \"", glomby, "\", featuresToKeep must contain ", glomby,
                           "-level names (the values in the \"", glomby, "\" column of the feature table, ",
                           "i.e. what rownames() would show AFTER agglomerating to ", glomby, "). ",
                           "A very common cause of a high miss rate here is passing names from a different ",
                           "level (e.g. LKT or CGB names) with glomby = \"", glomby, "\".")
                } else {
                    paste0("featuresToKeep is matched against the object's native feature names (rownames). ",
                           "Check that you are passing feature names from the correct analysis space and namespace.")
                }
            }

            show_absent <- function(){
                show_n <- min(n_absent, 10)
                paste0("Absent examples: ", paste0(utils::head(absent, show_n), collapse = ", "),
                       if (n_absent > show_n) ", ..." else "")
            }

            #### Decision tree ####
            if (n_match == 0){
                #Nothing to plot, and the clearest possible signal of a bad vector. Always stop.
                flog.warn(paste0("featuresToKeep: 0 of ", n_req, " requested feature(s) are present in this object."))
                if (n_absent > 0) flog.warn(show_absent())
                stop(paste0("featuresToKeep matched NONE of the features present in this object. ", namespace_hint()))

            } else if (n_req <= featuresToKeep_size_floor){
                #Small request: absence is plausible and benign. Never error; just note misses.
                if (n_absent > 0){
                    flog.info(paste0("featuresToKeep: ", n_match, " of ", n_req,
                                     " requested feature(s) are present and will be kept; ",
                                     n_absent, " not present in this object and ignored. ", show_absent()))
                }

            } else if (absent_frac >= featuresToKeep_error_frac){
                #Large request, majority missing: almost certainly a mistake. Stop.
                flog.warn(paste0("featuresToKeep: only ", n_match, " of ", n_req, " requested feature(s) (",
                                 round((1 - absent_frac) * 100, 1), "%) are present in this object; ",
                                 n_absent, " (", round(absent_frac * 100, 1), "%) are absent."))
                flog.warn(show_absent())
                stop(paste0("featuresToKeep: ", round(absent_frac * 100, 1),
                            "% of the (", n_req, ") requested features are absent from this object, ",
                            "which is above the ", round(featuresToKeep_error_frac * 100, 1),
                            "% error threshold for a request of this size. This usually means a wrong ",
                            "namespace or a stale feature list. ", namespace_hint(),
                            " If you intend to request a sparse set and want this to be a warning rather ",
                            "than an error, raise featuresToKeep_error_frac."))

            } else if (absent_frac >= featuresToKeep_warn_frac){
                #Large request, substantial-but-not-majority missing: loud warning, still proceed.
                flog.warn(paste0("featuresToKeep: ", n_match, " of ", n_req, " requested feature(s) present (",
                                 round((1 - absent_frac) * 100, 1), "%); ", n_absent, " (",
                                 round(absent_frac * 100, 1), "%) absent. Proceeding, but this miss rate is high - ",
                                 "verify you are passing names in the correct namespace. ", namespace_hint()))
                flog.warn(show_absent())

            } else if (n_absent > 0){
                #Few absentees: quietly informative.
                flog.info(paste0("featuresToKeep: ", n_match, " of ", n_req,
                                 " requested feature(s) present and kept; ", n_absent,
                                 " absent and ignored. ", show_absent()))
            }

            #Need at least 2 features to build anything meaningful downstream.
            if (n_match < 2){
                stop(paste0("After matching featuresToKeep, fewer than 2 features remain (", n_match,
                            "). Impossible to build a SummarizedExperiment with fewer than 2 features. ",
                            "Check featuresToKeep against the features actually present in this object/subset."))
            }

            obj <- obj[matched, ]
        }

        obj <- suppressWarnings(filter_sample_by_class_to_ignore(SEobj = obj, variables = variables_to_fix, class_to_ignore = class_to_ignore))

    return(obj)
}