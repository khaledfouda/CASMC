setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")
source("./BIXI/data/transform_data_into_mat_other.R")

results <- data.frame()

missp = 0.2
all_miss = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.9)
# all_miss = all_miss[1:2]
for(missp in all_miss)
{
dat <- load_model_bixi_dat2(time_cov = TRUE,2023,.2, missp, .2) 

out0 <- SImpute_Sim_Wrapper(dat) 

out1 <- CAMC_Ridge_Sim_Wrapper(dat, trace = F)

hpar = CAMC_Nuclear_hparams
#hpar$beta$n.lambda = 80
out2 <- CAMC_Nuclear_Sim_Wrapper(dat, trace=F, hpar=hpar)


hpar = CAMC_Lasso_hparams
hpar$beta$n.lambda = 40
hpar$beta$lambda.max = 3
out3 <- CAMC_Lasso_Sim_Wrapper(dat, trace=T, hpar = hpar)

out4 <- Naive_Sim_Wrapper(dat)

rbind(out0$results, out1,out2, out3, out4) %>% as.data.frame() %>% 
 mutate(miss = sum(is.na(dat$Y)|dat$Y==0)/length(dat$Y)) %>% 
 dplyr::select(model, miss, lambda.beta, lambda.M,
                          error.test, corr.test, error.train,
                          rank_M, rank_beta, sparse_all) %>% 
 mutate(error.test.diff = 
         (out0$results$error.test - as.numeric(error.test)) %>% 
         round(3)) %>% 
 mutate_at(vars(-model), function(x) round(as.numeric(x),3) )->
  all.out



results <- rbind(results, all.out)
print(results)
}

# for(missrows in  (0:length(all_miss) * 4)+1 ){
#  results[missrows:(missrows+3), ] %>% 
#   arrange(desc(corr.test)) %>% 
#  kable(format = "html", escape = FALSE) %>%
#   kable_styling(full_width = F, position = "left") %>%
#   print()
#  cat("<hr style='border:1px solid black;'>\n")
# }

unique_miss <- unique(results$miss)
color_palette <- scales::hue_pal()(length(unique_miss))
miss_color_map <- setNames(color_palette, unique_miss)



results %>% 
 arrange(miss, desc(corr.test)) %>%
 kable(format = "html", escape = FALSE) %>%
 kable_styling(full_width = F, position = "left") %>%
 row_spec(0, bold = TRUE) ->
 kable_table

for (val in unique_miss) {
 rows_to_color <- which(results$miss == val)
 kable_table <- kable_table %>% 
  row_spec(rows_to_color, background = miss_color_map[as.character(val)])
}

kable_table
