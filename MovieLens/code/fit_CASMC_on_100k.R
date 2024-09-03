setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")

source("./code_files/import_lib.R")
source("./MovieLens/code/load_data_100k.R")


dat <- load_movielens_100k("a", validp=0.1, seed = 2023)


out0 <- SImpute_Sim_Wrapper(dat) 

hpar_r <- CASMC_Ridge_hparams 
#hpar_r$beta$lambda.grid <- seq(10,7, length.out=20)
out1 <- CASMC_Ridge_Sim_Wrapper(dat, trace = T, 
                                hpar=hpar_r, return_fit = T, max_cores = 20)

hpar_n = CASMC_Nuclear_hparams
#hpar$beta$n.lambda = 80
out2 <- CASMC_Nuclear_Sim_Wrapper(dat, trace=T, hpar=hpar_n, return_fit = T)


hpar_l = CASMC_Lasso_hparams
#hpar_l$beta$n.lambda = 60
#hpar_l$beta$lambda.max = 1
out3 <- CASMC_Lasso_Sim_Wrapper(dat, trace=T, hpar = hpar_l, return_fit = T)

out4 <- Naive_Sim_Wrapper(dat)

out5 <- Mao_Sim_Wrapper(dat)


rbind(out0$results, out1$results, out2$results,
      out3$results, out4, out5) %>% 
 as.data.frame() %>% 
 mutate(total_miss = (sum(is.na(dat$Y)|dat$Y==0)/length(dat$Y) )%>% round(2),
        test_prop  = (sum(dat$W==0)/length(dat$Y)) %>% round(2)) %>% 
 dplyr::select(model, time, lambda.beta, lambda.M,
               error.test, corr.test, error.train,
               rank_M, rank_beta, sparse_all, total_miss, test_prop) %>% 
 # mutate(error.test.diff = 
 #        out0$results$error.test - as.numeric(error.test),
 #        corr.test.diff = 
 #         out0$results$corr.test - as.numeric(corr.test)) %>% 
 
 mutate_at(vars(-model), function(x) round(as.numeric(x),3) )->
 all.out
rownames(all.out) <- NULL
all.out %>% arrange(error.test, desc(corr.test)) %>%
 mutate(corr.test = paste0(  round(as.numeric(corr.test)*100,1),"%" )) %>% 
 select(-test_prop, -total_miss) %>% 
 kable("html", escape = F,
       caption = paste("total missing=",all.out$total_miss[1],
                       ", test proportion=", all.out$test_prop[1])) %>% 
 kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
               full_width = F,
               position = "center") #%>%
 column_spec(
  5:6,
  bold = T,
  background = grDevices::adjustcolor("pink", alpha.f = 0.1)
 ) %>% 
 row_spec(
  c(1,2,12,10),
  bold = T,
  background = grDevices::adjustcolor("yellow", alpha.f = 0.05)
 )
