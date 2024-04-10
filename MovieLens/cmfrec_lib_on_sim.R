
# setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
# suppressMessages(suppressWarnings(source("./code_files/import_lib.R", local = FALSE)))
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./MovieLens/load_data.R")
source("./code_files/import_lib.R")



# rm(list = ls(all.names = TRUE))
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./MovieLens/load_data.R")
source("./code_files/import_lib.R")
source("./MovieLens/apply_cmref_on_simulation.R")

dat <-generate_simulation_data_mao(800,800,10,10, "MAR", 2024)
apply_to_sim_dat(dat)

# rm(list = ls(all.names = TRUE))

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./MovieLens/load_data.R")
source("./code_files/import_lib.R")
source("./MovieLens/apply_cmref_on_simulation.R")

dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.9,coll=T)
apply_to_sim_dat(dat)
