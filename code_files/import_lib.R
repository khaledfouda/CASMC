



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
path_to_code = paste0(path_to_proj, "functions/")
path_to_data = paste0(path_to_proj, "saved_data/")

file_list <-
 list.files(
  path = path_to_code,
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = T
 )
file_list <- file_list[!grepl("main\\.R$", file_list)]
file_list <- file_list[!grepl("/old/", file_list)]
file_list <- file_list[!grepl("/old.*/", file_list)]

. = lapply(file_list, source)




select <- dplyr::select
solve <- MASS::ginv


# source(paste0(path_to_code,"Mao_fit.R"))
# source(paste0(path_to_code,"Mao_cv.R"))
# source(paste0(path_to_code,"Mao_sim.R"))
# source(paste0(path_to_code,"Ysf_sim.R"))
# source(paste0(path_to_code,"graph_estim.R"))
# source(paste0(path_to_code,"graph_estim_2methods.R"))
# source(paste0(path_to_code,"matrix_split_train_valid.R"))
# source(paste0(path_to_code,"SoftImpute_fit_covariates.R"))
# source(paste0(path_to_code,"SoftImputeALS_fit_covariates.R"))
# source(paste0(path_to_code,"SoftImpute_cv_covariates.R")) # fit to lambda2 only with lambda1 fixed
# source(paste0(path_to_code,"SoftImpute_cv_covariates_v4.R")) # fit to lambda 1 and 2 at the same time
# source(paste0(path_to_code,"SoftImpute_cv_covariates_v5.R")) # fit to lambda1 after knowing lambda2
# source(paste0(path_to_code,"SoftImpute_cv_covariates_v6.R")) # all previous with k-fold
