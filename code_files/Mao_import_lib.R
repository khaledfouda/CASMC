library(knitr)
library(kableExtra)
library(tidyverse)
library(magrittr)
require(foreach)
require(doParallel)
library(ggplot2)
library(hexbin)
library(patchwork)


#setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/")
path_to_code = "/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/code_files/"

source(paste0(path_to_code,"Mao_fit.R"))
source(paste0(path_to_code,"Mao_cv.R"))
source(paste0(path_to_code,"Mao_sim.R"))
source(paste0(path_to_code,"Ysf_sim.R"))
source(paste0(path_to_code,"graph_estim.R"))