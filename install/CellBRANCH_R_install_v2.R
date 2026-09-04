# Fresh conda env: pak-first installer
options(repos = c(CRAN = "https://cran.yu.ac.kr"))

# 1) pak
install.packages("pak")
library(pak)

# 2) Load lists
source("R_package_dependency.R")  # defines bioc_packages, cran_packages

# 3) Bioconductor
bioc_missing <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(bioc_missing)) {
  pak::pkg_install(paste0("bioc::", bioc_missing), ask = FALSE, upgrade = FALSE)
}

# 4) CRAN
cran_missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_missing)) {
  pak::pkg_install(cran_missing, ask = FALSE, upgrade = FALSE)
}

# 5) Pins / special cases
# 5.1 irlba (force source build)
if (!requireNamespace("withr", quietly = TRUE)) install.packages("withr")
withr::with_envvar(c(PAK_NO_BINARY = "true"), {
  pak::pkg_install("irlba@2.3.5.1", ask = FALSE, upgrade = FALSE)
})

# 5.2 dlm (newer works with scATOMIC on R 4.4)
pak::pkg_install("dlm@1.1-6.1", ask = FALSE, upgrade = FALSE)

# 5.3 Rmagic (archived)
pak::pkg_install(
  "Rmagic=url::https://cran.r-project.org/src/contrib/Archive/Rmagic/Rmagic_2.0.3.tar.gz",
  ask = FALSE, upgrade = FALSE
)

# 6) GitHub pins
pak::pkg_install(c(
  "inofechm/cutoff.scATOMIC@a43fd5e",
  "abelson-lab/scATOMIC@v2.0.3"
), ask = FALSE, upgrade = FALSE)

# 7) (Optional) version manifest for reproducibility
# ver <- as.data.frame(installed.packages()[, c("Package", "Version")])
# write.csv(ver, "cellbranch_base_versions.csv", row.names = FALSE)

# Functioning as of 2025-08-12 on Linux.