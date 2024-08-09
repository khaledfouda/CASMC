setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")
source("./BIXI/data/transform_data_into_mat_other_2.R")


dat <- load_model_bixi_dat3(time_cov = TRUE,2023,.2) 

out0 <- SImpute_Sim_Wrapper(dat) 

out1 <- CASMC_Ridge_Sim_Wrapper(dat, trace = T, return_fit = T, max_cores = 20)

hpar = CASMC_Nuclear_hparams
#hpar$beta$n.lambda = 80
out2 <- CASMC_Nuclear_Sim_Wrapper(dat, trace=T, hpar=hpar, return_fit = T)


hpar = CASMC_Lasso_hparams
hpar$beta$n.lambda = 40
hpar$beta$lambda.max = 3
out3 <- CASMC_Lasso_Sim_Wrapper(dat, trace=T, hpar = hpar, return_fit = T)

out4 <- Naive_Sim_Wrapper(dat)

rbind(out0$results, out1$results,out2$results, out3$results, out4) %>% 
 as.data.frame() %>% 
 mutate(miss = sum(is.na(dat$Y)|dat$Y==0)/length(dat$Y)) %>% 
 dplyr::select(model, miss, lambda.beta, lambda.M,
                          error.test, corr.test, error.train,
                          rank_M, rank_beta, sparse_all) %>% 
 mutate(error.test.diff = 
         out0$results$error.test - as.numeric(error.test),
          corr.test.diff = 
            out0$results$corr.test - as.numeric(corr.test)) %>% 
          
 mutate_at(vars(-model), function(x) round(as.numeric(x),3) )->
  all.out

all.out %>% arrange(desc(corr.test)) %>%  kable()
#-----------------------------------------------------
X <- dat$X

apply(out1$fit$fit$beta, 1, summary) |> as.data.frame() |>
  t() |>
  as.data.frame() |>
  mutate(prop_non_zero = apply(out1$fit$fit$beta, 1, function(x)
    sum(x != 0) / length(x))) |>
  `rownames<-` (colnames(X)) %>% 
  mutate(Model = out1$results$model) %>% 

    rbind(
    

apply(out2$fit$fit$beta, 1, summary) |> as.data.frame() |>
  t() |>
  as.data.frame() |>
  mutate(prop_non_zero = apply(out2$fit$fit$beta, 1, function(x)
    sum(x != 0) / length(x))) |>
  `rownames<-` (colnames(X)) %>% 
  mutate(Model = out2$results$model)
  ) %>% 

  rbind(
apply(out3$fit$fit$beta, 1, summary) |> as.data.frame() |>
  t() |>
  as.data.frame() |>
  mutate(prop_non_zero = apply(out3$fit$fit$beta, 1, function(x)
    sum(x != 0) / length(x))) |>
  `rownames<-` (colnames(X)) %>% 
  mutate(Model = out3$results$model)
  ) %>% 
  mutate_if(is.numeric, round, digits=3) %>% 
  kable("html", col.names = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "Prop Non Zero", "Model"),
        caption = "Covariate Summaries by Model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F, position = "center") %>%
  row_spec(0, bold = TRUE) %>%
  # pack_rows("Mao", 1, 3, hline_after = TRUE) %>%
  pack_rows("CASMC-Ridge", 1, 3, hline_after = TRUE) %>%
  pack_rows("CASMC-Nuclear", 4, 6, hline_after = TRUE) %>%
  pack_rows("CASMC-Lasso", 7, 9, hline_after = TRUE) 







#results <- rbind(results, all.out)
#print(results)
#}

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
