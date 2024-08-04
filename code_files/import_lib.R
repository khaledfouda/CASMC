



packages <-
 c(
  "tidyverse",
  "kableExtra",
  "tidyverse",
  "doParallel",
  "foreach",
  "hexbin",
  "patchwork",
  "ggplot2",
  "softImpute",
  "irlba", # propack
  "MASS", # ginv (alternative to solve())
  "svd", # another propack. i forgot which one i ended up using :)
  "corpcor", # fast.svd
  "magrittr",
  "RSpectra",
  "Matrix",
  "knitr",
  "roxygen2" # to create help pages
 )

# unload the packages if they're loaded. 
for (pkg in packages) {
 name = paste0("package:",pkg)
 if (name %in% search()) {
  suppressMessages(suppressWarnings(detach(
   name,
   unload = TRUE,
   character.only = TRUE,
   force = F
  )))
 }
}
# load the packages
for (pkg in packages) {
 print(paste("Package", pkg, "loaded."))
 suppressMessages(suppressWarnings(library(
  pkg, character.only = TRUE, quietly = TRUE
 )))
}



path_to_proj = "/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/"
# setwed(path_to_proj)
path_to_code = paste0(path_to_proj, "functions2/")
path_to_data = paste0(path_to_proj, "saved_data/")

file_list <-
 list.files(
  path = path_to_code,
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = T
 )
file_list <- file_list[!grepl("main\\.R$", file_list)]
file_list <- file_list[!grepl("main[0-9]*\\.R$", file_list)]
file_list <- file_list[!grepl("/old/", file_list)]
file_list <- file_list[!grepl("/old.*/", file_list)]

. = lapply(file_list, source)




select <- dplyr::select
solve <- MASS::ginv









