setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")
select =  dplyr::select

file_list <- list.files(path = path_to_data, pattern = "Compare_MC_Models.*\\.csv$",
                        full.names = TRUE,recursive =F)

results = data.frame()
for(f in file_list)
   results = rbind(results, read.csv(f))

results %<>%
   mutate(simulation = case_when(missinginess == 0 ~ "Mao",
                                 missinginess == 0.8 & collinearity == TRUE ~ "80% Missing with Coll",
                                 missinginess == 0.8 & collinearity == FALSE ~ "80% Missing",
                                 missinginess == 0.9 & collinearity == TRUE ~ "90% Missing with Coll",
                                 missinginess == 0.9 & collinearity == FALSE ~ "90% Missing")) %>%
   filter(simulation != "90% Missing") 

# 1. plot time
results %>%
   select(time, dim, model, simulation) %>%
   ggplot(aes(x = dim, y = time, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "Time (seconds)", 
      color = "Model",  
      title = "Model Running Time by Dimension and Simulation Scenario",
      subtitle = "Comparing time used for fitting and computing predictions across different models"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()



# 2. plot A hat
results %>%
   select(error.test, dim, model, simulation) %>%
   ggplot(aes(x = dim, y = error.test, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "1 - R^2", 
      color = "Model",  
      title = "Model Performance on the Test Dataset (A) by Dimension and Simulation Scenario",
      subtitle = "Comparing 1 - goodness of fit  across different models"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()



# 2. plot A hat
results %>%
   select(error.beta, dim, model, simulation) %>%
   filter(! is.na(error.beta)) %>% 
   ggplot(aes(x = dim, y = error.beta, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "1 - R^2", 
      color = "Model",  
      title = "Model Performance on the Covariate Coefficients (Beta) by Dimension and Simulation Scenario",
      subtitle = "Comparing 1 - goodness of fit  across different models"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()


# 2. plot B hat
results %>%
   select(error.B, dim, model, simulation) %>%
   filter(! is.na(error.B)) %>% 
   ggplot(aes(x = dim, y = error.B, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "1 - R^2", 
      color = "Model",  
      title = "Model Performance on the Low Rank Matrix (B) by Dimension and Simulation Scenario",
      subtitle = "Comparing 1 - goodness of fit  across different models"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()


# 2. plot Rank
results %>%
   select(rank, true_rank, dim, model, simulation) %>%
   ggplot(aes(x = dim, y = rank, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_line(aes(y=true_rank), color="black", size=1, linetype="dashed") +
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "Rank", 
      color = "Model",  
      title = "Rank of Estimated Matrix by Dimension and Simulation Scenario",
      subtitle = "Dashed lines indicate the true rank"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()


# 2. plot Rank
results %>%
   select(lambda.2, dim, model, simulation) %>% 
   filter(!is.na(lambda.2), model != "Mao", model != "SoftImpute_Orig") %>% 
   ggplot(aes(x = dim, y = lambda.2, color = model, group = model)) +
   geom_line(size = 1.2) + 
   geom_point(size = 3, shape = 21, fill = "white") +  
   scale_color_brewer(palette = "Dark2") +  
   facet_wrap(~simulation, scales="free_y") +
   #theme_minimal(base_size = 14) + 
   labs(
      x = "Dimension", 
      y = "Parameter Value", 
      color = "Model",  
      title = " Rank Restriction Regularization Parameter (lambda_2) by Dimension and Simulation Scenario",
      subtitle = "Only the 4 soft-Impute with covariate models are considered"
   ) +
   theme(
      legend.position = "bottom", 
      panel.spacing = unit(1, "lines"),  
      strip.text.x = element_text(face = "bold")  
   ) +
   ggthemes::theme_stata()
