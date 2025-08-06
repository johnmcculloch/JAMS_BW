# Function for detecting system RAM
get_system_ram <- function() {
    os <- as.character(Sys.info()["sysname"])

    if (os == "Linux") {
        return(get_linux_ram())
    } else if (os == "Darwin") {
        return(get_macos_ram())
    } else if (os == "Windows") {
        return(get_windows_ram())
    } else if (os %in% c("FreeBSD", "OpenBSD", "NetBSD", "SunOS")) {
        # These might work with Linux-style /proc/meminfo
        return(get_linux_ram())  # Try Linux method as fallback
    } else {
        stop("Unsupported OS, please put in a bug report at https://github.com/johnmcculloch/JAMS_BW/issues. OS:", os)
    }
}

# Detect RAM
get_linux_ram <- function() {
    ram_kb <- system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)
    return(as.numeric(ram_kb) * 1024L)
}

get_macos_ram <- function() {
    ram_bytes <- system("sysctl -n hw.memsize", intern = TRUE)
    return(as.numeric(ram_bytes))
}

get_windows_ram <- function() {
    stop("Windows functionality will be coming shortly! Apologies for the delay, and please contact the developers to use JAMS on windows. (https://github.com/johnmcculloch/JAMS_BW)")
}

