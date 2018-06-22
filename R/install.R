sapply(
    c("callr","raster","data.table","maptools","readr",
      "RColorBrewer","ggplot2","mgcv","snow","dplyr","ape",
      "parallel","taxonlookup","xtable","treemap","ggtree"),
    function(x) if(!require(x)) install.packages(x)
)


source("R/data_manipulation_functions.R")
source("R/family_analysis.R")
source("R/gbif_process.R")
source("R/overlap.R")
source("R/treemap.R")
