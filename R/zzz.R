# Function for detecting if the user has an outdated version of JAMS
.onAttach <- function(libname, pkgname) {
    
    # find dir and run git fetch
    jams_dir <- packageDescription("JAMS")$RemoteUrl
    
    if (is.null(jams_dir)) {
        jams_dir <- file.path(Sys.getenv("HOME"), "bin", "JAMS_BW")
        if (!file.exists(jams_dir)) {
            jams_dir <- file.path(Sys.getenv("HOME"), "JAMS_BW")
            if (!file.exists(jams_dir)) {
                return(invisible(NULL))
        }
    }
    }
    
    
    if (!is.null(jams_dir)) {
        fetch_output <- system2('git', c('-C', jams_dir, 'fetch'), stdout=TRUE, stderr=TRUE)
        current_hash <- system2('git', c('-C', jams_dir, 'rev-parse', 'HEAD'), stdout=TRUE, stderr=TRUE)
        master_hash <- system2('git', c('-C', jams_dir, 'rev-parse', 'origin/master'), stdout=TRUE, stderr=TRUE)

        # Check if the two hashes are the same

        if (current_hash != master_hash) {
            packageStartupMessage("An updated version of JAMS is available! Please follow the instructions here to update https://github.com/johnmcculloch/JAMS_BW/wiki/JAMS2_recipe_and_best_practices. If you don't, the angels shall weep for you.")
        }
    }

}

