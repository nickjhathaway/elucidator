install_cran <- function(pkgs, ncpus = parallel::detectCores(),
                         repos = getOption("repos")) {
  install.packages(pkgs, Ncpus = ncpus, repos = repos)
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    cat("ERROR: Failed to install CRAN packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
  }
  cat("Successfully installed:", paste(pkgs, collapse = ", "), "\n")
}

install_bioconductor <- function(pkgs, version = NULL) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!is.null(version)) BiocManager::install(version = version, ask = FALSE)
  BiocManager::install(pkgs, ask = FALSE)
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    cat("ERROR: Failed to install Bioconductor packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
  }
  cat("Successfully installed:", paste(pkgs, collapse = ", "), "\n")
}

install_github <- function(repos) {
  if (!require("remotes", quietly = TRUE)) install.packages("remotes")
  # Extract expected package names from "user/pkg@ref" or "user/pkg"
  pkg_names <- sapply(strsplit(sapply(strsplit(repos, "/"), `[`, 2), "@"), `[`, 1)
  remotes::install_github(repos)
  missing <- pkg_names[!pkg_names %in% rownames(installed.packages())]
  if (length(missing) > 0) {
    cat("ERROR: Failed to install GitHub packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
  }
  cat("Successfully installed:", paste(pkg_names, collapse = ", "), "\n")
}
