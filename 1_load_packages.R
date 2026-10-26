# Package names

# 1. Define all packages
packages <- c(
  "plot3D", "fda", "mgcv", "HMDHFDplus", 
  "demography", "ftsa", "hdftsa", "sde", 
  "xtable", "dplyr", "quantreg"
)

# 2. Install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# 3. Load packages
# require(plot3D)
require(fda)
require(mgcv)
require(HMDHFDplus)
require(demography)
require(ftsa)
require(hdftsa)
require(sde)
require(xtable)
require(dplyr)
require(quantreg)
