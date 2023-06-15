# Installs the required packages from CRAN and GitHub.

cran <- c("dplyr", "data.table", "remotes", "openxlsx", "rmarkdown", "plotly")
for (pkg in cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

github <- c("abbvie-external/OmicNavigator@*release")
remotes::install_github(github, dependencies = TRUE, upgrade = FALSE)
OmicNavigator::installApp()