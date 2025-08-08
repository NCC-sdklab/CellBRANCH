# ------------------------------------------------------------------------
# R Package Installation Script (for conda-based or isolated environments)
# ------------------------------------------------------------------------
#Tested on: 2025.07.04


# 1. Set CRAN mirror (currently set to Korea)
options(repos = c(CRAN = "https://cran.yu.ac.kr"))

# 2. Ensure BiocManager is installed
if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

# 3. Update already-installed Bioconductor packages only
BiocManager::install(ask = FALSE, update = TRUE)

# 4. Load list of dependencies
source("R_package_dependency.R")  # defines `bioc_packages` and `cran_packages`

# 5. Install only missing Bioconductor packages
missing.bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing.bioc) > 0) {
  BiocManager::install(missing.bioc, ask = FALSE, update = FALSE)
}

# 6. Install only missing CRAN packages
missing.cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing.cran) > 0) {
  install.packages(missing.cran)
}

# 7. Workaround for `Matrix`-related compatibility (see Signac issue #1596)
install.packages("irlba", type = "source")

# 8. Github install
if(!require(devtools)) install.packages("devtools")
devtools::install_version("dlm", version = "1.1.5", repos = "http://cran.us.r-project.org")
devtools::install_version("Rmagic", version = "2.0.3", repos = "http://cran.us.r-project.org")
if(!require(cutoff.scATOMIC)) devtools::install_github("inofechm/cutoff.scATOMIC", force = T)
if(!require(scATOMIC)) devtools::install_github("abelson-lab/scATOMIC")