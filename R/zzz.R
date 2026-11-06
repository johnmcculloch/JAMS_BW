# Function for detecting if the user has an outdated version of JAMS.
# This runs every time library(JAMS) is called, so it must never fail and never hang:
#  * laptops (Windows, or macOS without the Xcode command line tools) may not have git at all;
#  * HPC compute nodes often have no outbound internet, so "git fetch" could block for minutes;
#  * JAMS may have been installed from a tarball or GitHub URL rather than from a local clone.
# Any of these situations silently skips the check. Set options(JAMS.check_updates = FALSE) or the
# environment variable JAMS_NO_UPDATE_CHECK=1 to switch it off altogether.
.onAttach <- function(libname, pkgname) {

    if (isFALSE(getOption("JAMS.check_updates", TRUE)) || nzchar(Sys.getenv("JAMS_NO_UPDATE_CHECK"))) {
        return(invisible(NULL))
    }

    # never let git sit waiting for a password
    old_prompt <- Sys.getenv("GIT_TERMINAL_PROMPT", unset = NA)
    Sys.setenv(GIT_TERMINAL_PROMPT = "0")
    on.exit(if (is.na(old_prompt)) Sys.unsetenv("GIT_TERMINAL_PROMPT") else Sys.setenv(GIT_TERMINAL_PROMPT = old_prompt), add = TRUE)

    tryCatch({
        git <- Sys.which("git")
        if (!nzchar(git)) {
            return(invisible(NULL))
        }

        # find the local clone JAMS was installed from
        candidate_dirs <- c(packageDescription(pkgname)$RemoteUrl,
                            file.path(Sys.getenv("HOME"), "bin", "JAMS_BW"),
                            file.path(Sys.getenv("HOME"), "JAMS_BW"))
        candidate_dirs <- candidate_dirs[!is.na(candidate_dirs) & nzchar(candidate_dirs)]
        candidate_dirs <- candidate_dirs[dir.exists(file.path(candidate_dirs, ".git"))]
        if (length(candidate_dirs) == 0) {
            return(invisible(NULL))
        }
        jams_dir <- candidate_dirs[1]

        # git fetch with a short timeout so that offline nodes do not hang library(JAMS)
        fetch_status <- suppressWarnings(system2(git, c("-C", shQuote(jams_dir), "fetch", "--quiet"), stdout = FALSE, stderr = FALSE, timeout = 10))
        if (!identical(as.integer(fetch_status), 0L)) {
            return(invisible(NULL))
        }
        current_hash <- suppressWarnings(system2(git, c("-C", shQuote(jams_dir), "rev-parse", "--verify", "--quiet", "HEAD"), stdout = TRUE, stderr = FALSE, timeout = 10))
        master_hash <- suppressWarnings(system2(git, c("-C", shQuote(jams_dir), "rev-parse", "--verify", "--quiet", "origin/master"), stdout = TRUE, stderr = FALSE, timeout = 10))
        # a failed rev-parse carries a non-zero "status" attribute; skip the message then
        if (!is.null(attr(current_hash, "status")) || !is.null(attr(master_hash, "status"))) {
            return(invisible(NULL))
        }

        # Check if the two hashes are the same
        if ((length(current_hash) == 1) && (length(master_hash) == 1) && (current_hash != master_hash)) {
            packageStartupMessage("An updated version of JAMS is available! Please follow the instructions here to update https://github.com/johnmcculloch/JAMS_BW/wiki/JAMS2_recipe_and_best_practices. If you don't, the angels shall weep for you.")
        }
    }, error = function(e) invisible(NULL))

    invisible(NULL)
}
