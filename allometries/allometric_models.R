# ==============================================================================
# Script Name:        Allometric Models
# Description:        workflow for fitting and comparing allometric models 
#                     (log-log linear and nonlinear) for tree height–DBH 
#                     relationships
#
# Author:             Hanna Komischke
# Date:               [2025-11-18]
#
# Input Data:         - BWI data (CSV)
#                     - iLand species database (SQLite)
# Output Data:        - Allometric model coefficients (CSV)
#                     - Plots (PNG)
#
# ==============================================================================

# directories -------------------------------------------------------------

setwd("/media/hanna/EXTERNAL_USB/Master/data/TreeSpecies/")

plotdir <- "/home/hanna/MA/Master/Orga/plots/"

# packages ----------------------------------------------------------------

# Load required packages
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(caret)
library(minpack.lm)
library(RSQLite)
library(stringr)

# functions ---------------------------------------------------------------

# calc performance metrics
model_performance <- function(df_pred) {
  df_perf <- df_pred %>%
    group_by(genus) %>%
    summarise(
      rmse = sqrt(mean((dbh - pred_dbh)^2)),
      bias = (mean(pred_dbh) - mean(dbh)) / mean(dbh) * 100,
      r2   = 1 - sum((dbh - pred_dbh)^2) / sum((dbh - mean(dbh))^2),
      .groups = "drop"
    ) 
  
  return(df_perf)
}

# Fun 1: Log-Log Linear Regression (deterministic) 
jucker_det <- function(df, coef=F) {
  # Fit allometric model per genus
  allom_params <- df %>%
    group_by(genus) %>%
    do({
      fit <- lm(log(height) ~ log(dbh), data = .)
      a <- exp(coef(fit)[1])  # intercept 
      b <- coef(fit)[2]       # slope
      data.frame(a = a, b = b)
    }) %>%
    ungroup()
  
  if(coef==T){
    return(allom_params)
  }else{
    # Join coefficients to data
    df_coef <- df %>% left_join(allom_params, by = "genus")
    
    # Predict DBH
    df_pred <- df_coef %>%
      mutate(pred_dbh = (height / a)^(1 / b))
    
    return(df_pred)
  }
}

# Fun 2: Nonlinear Regression (deterministic) 
parker_det <- function(df_all,df, coef=F) {
  
  # global start values
  start_vals <- df_all %>%
    group_by(genus) %>%
    summarise(
      b0_start = min(dbh),
      b1_start = (max(dbh) - min(dbh)) / 2,
      b2_start = 0.5
    )
  
  # Fit non-linear model per genus using nlsLM
  models <- df %>%
    group_by(genus) %>%
    group_map(~ {
      sv <- start_vals %>% filter(genus == unique(.x$genus))
      
      nls_fit <- tryCatch(
        nlsLM(dbh ~ b0 + b1 * (height)^b2,
              data = .x,
              start = list(
                b0 = sv$b0_start,
                b1 = sv$b1_start,
                b2 = sv$b2_start
              ),
              lower = c(-Inf, 0, -Inf),
              upper = c(Inf, Inf, Inf),
              control = nls.lm.control(maxiter = 500)),
        error = function(e) NULL
      )
      tibble(genus = unique(.x$genus), model = list(nls_fit))
    }, .keep = TRUE) %>%
    bind_rows()
  
  # Extract coefficients
  coef_df <- models %>%
    mutate(coef = map(model, ~ as.list(coef(.x)))) %>%
    select(genus, coef) %>%
    unnest_wider(coef)
  
  if(coef==T){
    return(coef_df)
  } else {
    # Join coefficients to data
    df_coef <- df %>% left_join(coef_df, by = "genus")
    
    # Predict DBH
    df_pred <- df_coef %>%
      mutate(pred_dbh = pmax(b0 + b1 * (height)^b2, 7))
    
    return(df_pred)
  }
}


# ==============================================================================
# data =========================================================================
# ==============================================================================

# species of interest -----------------------------------------------------

genus_list <- data.frame(
  genus = c("Betula", "Fagus", "Pseudotsuga", "Quercus", "Alnus", 
            "Picea", "Pinus", "Larix", "Abies",
            "Acer", "Fraxinus","Populus", "Salix","Caprinus"),
  species_eng = c("Birch","Beech","Douglas fir","Oak","Alder",
                  "Spruce","Pine","Larch","Fir","Maple","Ash","Poplar","Willow","Hornbeam"))

# BWI dataset  --------------------------------------------------

# load data
bwi_data <- read.csv("BWI/bwi2012_allometry.csv", sep = ";") %>% na.omit()
bwi_codes <- read.csv("BWI/bwi2012_species_codes.csv", sep = ";") %>% 
  rename(species = ICode)

# filter species with > 500 samples
bwi_species <- bwi_data %>% 
  select(species) %>%
  group_by(species) %>% mutate(n = n()) %>% 
  unique() %>% left_join(bwi_codes, by = "species") %>% 
  filter(n > 500)


# add genus info
bwi_species$genus <- c("Acer","Fagus","Picea","Alnus2","Fraxinus","Abies",
                       "Sorbus","Pinus","Salix","Ulmus","Larix","Betula",
                       "Betula2","Quercus","Prunus","Quercus2","Caprinus",
                       "Pseudotsuga","Robinia","Alnus","Tilia","Acer2",
                       "Quercus2","Acer2","Populus","Larix2","Populus2",
                       "Pinus2","Populus2","Castanea","Picea2")

# filter for species of interest and correct units
bwi_vals <- bwi_data %>%
  filter(species %in% bwi_species$species) %>%
  mutate(dbh = dbh/10, height = height/10) %>%
  inner_join(bwi_species %>% select(species, genus)) %>%
  inner_join(genus_list, by= join_by(genus))

# iLand -------------------------------------------------------------------

# connect
db <- dbConnect(RSQLite::SQLite(), "all_species_database.sqlite")

# list tables
db_tab <- dbListTables(db)

# read species table
species <- dbReadTable(db,"species")

# disconnect
dbDisconnect(db)

# add englisch name
iland_species <- c("Birch","Beech","Douglas fir","Oak","Alder","Spruce","Pine","Larch","Fir","Maple","Ash","Hornbeam","Willow","Poplar")
iland_genus <- c("Betula", "Fagus", "Pseudotsuga", "Quercus", "Alnus", 
                 "Picea", "Pinus", "Larix", "Abies",
                 "Acer", "Fraxinus","Caprinus", "Salix", "Populus")
iland_species_code <- c("bepe","fasy","psme","quro","algl","piab","pisy","lade","abal","acps","frex","cabe","saca","potr")
iland_code <- data.frame("genus"= iland_genus,"species_eng" = iland_species,"shortName"= iland_species_code)

# filter species
species_f <- species %>% filter(shortName %in% iland_species_code) %>% inner_join(iland_code, join_by(shortName))

# correct functions  
species_f$HDhigh[species_f$species_eng =="Ash"] <- "(467.781*1.001*(1-0.4839)*1.2)*d^-0.4839"
species_f$HDhigh[species_f$species_eng =="Spruce"] <- "195.547*1.004*(-0.2396+1)*d^-0.2396"
species_f$HDhigh[species_f$species_eng =="Maple"] <- "(350.654*1.000*(1-0.4012)*1.2)*d^-0.4012"

# ==============================================================================
# allometric models ============================================================
# ==============================================================================

data <- bwi_vals 

# Performance metrics ----------------------------------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae  <- function(obs, pred) mean(abs(obs - pred))
r2   <- function(obs, pred) cor(obs, pred)^2


# Stratified K-Fold CV ------------------------------------------------

set.seed(42)
K <- 10

# Stratify with genus
folds <- createFolds(data$genus, k = K, list = TRUE, returnTrain = FALSE)

# CV Loop 
results <- map_dfr(1:K, function(k){
  
  test_idx  <- folds[[k]]
  trainData <- data[-test_idx, ]
  testData  <- data[test_idx, ]
  
  
#--- Fun 1: Log-Log Linear Regression 
  # Fit 
  jucker_train <- jucker_det(trainData)
  # Predict
  test_jucker <- testData %>%
    left_join(
      jucker_train %>% select(genus, a, b) %>% distinct(),
      by = "genus"
    ) %>%
    mutate(pred_dbh = (height / a)^(1 / b))
  
#--- Fun 2: Nonlinear Regression
  # Fit 
  parker_train <- parker_det(data,trainData)
  # Predict
  test_parker <- testData %>%
    left_join(parker_train %>% select(genus, b0, b1,b2) %>% distinct(),
              by = "genus"
    ) %>%
    mutate(pred_dbh = pmax(b0 + b1 * (height)^b2, min(trainData$dbh)))

  
#--- Performance    
  bind_rows(
    test_jucker %>%
      group_by(genus) %>%
      summarise(
        Fold = k,
        RMSE = rmse(dbh, pred_dbh),
        MAE  = mae(dbh, pred_dbh),
        R2   = r2(dbh, pred_dbh)
      ) %>%
      mutate(Model = "Jucker"),
    
    test_parker %>%
      group_by(genus) %>%
      summarise(
        Fold = k,
        RMSE = rmse(dbh, pred_dbh),
        MAE  = mae(dbh, pred_dbh),
        R2   = r2(dbh, pred_dbh)
      ) %>%
      mutate(Model = "Parker")
  )
})


# Intern allometric model of iLand ------------------------------------------------------------

# get median dbh
median_dbh <- bwi_vals %>%
  filter(genus %in% iland_genus) %>%
  group_by(species,genus) %>%
  summarise(median_dbh=median(dbh,na.rm=T))

# get coeficents
iland_dbh_coef <- species_f %>% 
  select(genus, HDlow, HDhigh) %>%
  pivot_longer(-genus, names_to = "HDiland") %>%
  mutate(a = str_extract(value, ".*(?=\\*[ds]\\^)"),
         b = str_extract(value, "(?<=\\*[ds]\\^)-?\\d+\\.\\d+")
  )%>%
  mutate(
    a = sapply(a, function(x) eval(parse(text = x))),  
    b = as.numeric(b)                                 
  ) %>% select(-value) %>%
  inner_join(median_dbh,by = join_by(genus))

# predict for high stand density
iland_det_bwi_high <- bwi_vals %>%
  filter(genus %in% iland_genus) %>%
  left_join(
    iland_dbh_coef %>%
      filter(HDiland == "HDhigh") %>%
      select(genus, a, b,median_dbh),
    by = "genus"
  )%>%
  mutate(hd = a*median_dbh^(b),
         pred_dbh = (height/hd)*100,
         Model="iland high")

# predict for low stand density
iland_det_bwi_low <- bwi_vals %>%
  filter(genus %in% iland_genus) %>%
  left_join(
    iland_dbh_coef %>%
      filter(HDiland == "HDlow") %>%
      select(genus, a, b,median_dbh),
    by = "genus"
  ) %>%
  mutate(hd = a*median_dbh^(b),
         pred_dbh = (height/hd)*100,
         Model="iland low")

# combine
iland_det_bwi <- iland_det_bwi_high %>% rbind(iland_det_bwi_low)

#--- Performance

iland_perf <- iland_det_bwi %>%
  #filter(Model=="iland high")%>%
  group_by(genus,Model) %>%
  summarise(
    RMSE = rmse(dbh, pred_dbh),
    MAE  = mae(dbh, pred_dbh),
    R2   = r2(dbh, pred_dbh)
  )

iland_perf




# Compare Performance -----------------------------------------------------

# combine
results_all <- results %>%
  rbind(iland_perf %>% mutate(Fold=1)) 

# Performance per Genus 
results_summary_per_genus <-  results_all %>%
  group_by(genus, Model) %>%
  summarise(across(c(RMSE, MAE, R2),\(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  inner_join(genus_list, by= join_by(genus))

# Select best deterministic model 
best_model_name <- results_summary_per_genus %>% group_by(Model) %>% summarise(across(c(RMSE, MAE, R2), ~ mean(.x, na.rm = TRUE))) %>% arrange(RMSE) 


# Get mean coefficients and predict -------------------------------------------------------------

#--- Fun 1: Log-Log Linear Regression 

# get coefficients
coefs_jucker <- map_dfr(1:K, function(k){
  trainData <- data[-folds[[k]], ]
  
  coefs <- jucker_det(trainData, coef = TRUE) %>%
    mutate(Fold = k)
  
  return(coefs)
})

# mean and sd per genus
coefs_summary_jucker <- coefs_jucker %>%
  group_by(genus) %>%
  summarise(
    a_mean = mean(a),
    a_sd   = sd(a),
    b_mean = mean(b),
    b_sd   = sd(b),
    .groups = "drop"
  ) %>%
  rename(a = a_mean, b = b_mean)

# predict on all data
pred_jucker <- bwi_vals %>%
  left_join(coefs_summary_jucker, by = "genus") %>%
  mutate(pred_dbh = (height / a)^(1 / b))

#--- Fun 2: Nonlinear Regression

# get coefficients
coefs_parker <- map_dfr(1:K, function(k){
  trainData <- data[-folds[[k]], ]
  
  coefs <- parker_det(data,trainData, coef = TRUE) %>%
    mutate(Fold = k)
  
  return(coefs)
})

# mean and sd per genus
coefs_summary_parker <- coefs_parker %>%
  group_by(genus) %>%
  summarise(
    b0_mean = mean(b0),
    b0_sd   = sd(b0),
    b1_mean = mean(b1),
    b1_sd   = sd(b1),
    b2_mean = mean(b2),
    b2_sd   = sd(b2),
    .groups = "drop"
  ) %>% rename(b0=b0_mean,b1=b1_mean,b2=b2_mean)

# predict on all data
pred_parker <- bwi_vals %>%
  left_join(coefs_summary_parker, by = "genus") %>%
  mutate(pred_dbh = pmax(b0 + b1 * (height)^b2, 7))


# add engl name and save
coef_df <- coefs_summary_parker %>% inner_join(genus_list,by=join_by(genus)) %>%
  rename(species=species_eng) %>%
  select(genus,b0,b1,b2,species)

write.csv(coef_df,"Parker_dbh_coef.csv",row.names = F)


# ==============================================================================
# plots ========================================================================
# ==============================================================================

df_plot <- bind_rows(
  bwi_vals %>% select(species_eng,genus,dbh,height) %>% mutate(source = "BWI")
)

df_pred_lines <- bind_rows(
  pred_jucker %>% mutate(Model = "log-log linear"),
  pred_parker %>% mutate(Model = "nonlinear"),
  iland_det_bwi
) %>% filter( Model != "iland low")


#--- only main species
selected_genera <- c("Fagus", "Pinus", "Picea","Quercus")

df_plot_sub <- df_plot %>% filter(genus %in% selected_genera)
df_pred_lines_sub <- df_pred_lines %>% filter(genus %in% selected_genera)

# plot
p_dbh_height <-  ggplot() +
  geom_point(data = df_plot_sub,
             aes(x = height, y = dbh),
             color="grey80",shape = 21, fill = "grey80", alpha = 1, size = 2) +
  geom_point(data = df_plot_sub,
             aes(x = height, y = dbh, fill = source),
             color="grey50",shape = 21, alpha = 0.005, size = 2) +
  geom_line(data = df_pred_lines_sub,
            aes(x = height, y = pred_dbh, color = Model),
            size = 1) +
  scale_x_continuous(limits = c(0,52))+
  scale_y_continuous(limits = c(0,200))+
  facet_wrap(~species_eng) +
  scale_fill_manual(guide="none",values = c("BWI" = "grey50")) +
  scale_color_manual(values = c("log-log linear" = "#92C5DE",
                                "nonlinear" = "#CC79A7",
                                "iland high" = "#FDB863",  
                                "iland low" = "#A6D854"),
                     labels= c("log-log linear" = "log-log linear",
                               "nonlinear" = "nonlinear",
                               "iland high" = "iLand",  
                               "iland low" = "iland open")) +
  labs(
    x = "Height (m)",
    y = "DBH (cm)",
    fill = "Dataset",
    color = "Model"
    #title = "Height vs DBH for selected genera with predicted model lines"
  ) +
  theme_minimal() +
  theme( element_text(size = 18),legend.position = "bottom")

p_dbh_height

ggsave(
  filename = paste0(plotdir,"allo_height_dbh_v3.png"),  
  plot = p_dbh_height,                 
  width = 8,                
  height = 5,               
  dpi = 300                 
)


#--- all species
p_dbh_height_all <-  ggplot() +
  geom_point(data = df_plot,
             aes(x = height, y = dbh),
             color="grey80",shape = 21, fill = "grey80", alpha = 1, size = 2) +
  geom_point(data = df_plot,
             aes(x = height, y = dbh, fill = source),
             color="grey50",shape = 21, alpha = 0.005, size = 2) +
  geom_line(data = df_pred_lines,
            aes(x = height, y = pred_dbh, color = Model),
            size = 1) +
  scale_x_continuous(limits = c(0,52))+
  scale_y_continuous(limits = c(0,200))+
  facet_wrap(~species_eng, nrow=7) +
  scale_fill_manual(guide="none",values = c("BWI" = "grey50")) +
  scale_color_manual(values = c("log-log linear" = "#92C5DE",
                                "nonlinear" = "#CC79A7",
                                "iland high" = "#FDB863",  
                                "iland low" = "#A6D854"), 
                     labels= c("log-log linear" = "log-log linear",
                               "nonlinear" = "nonlinear",
                               "iland high" = "iLand",  
                               "iland low" = "iland open")) +
  labs(
    x = "Height (m)",
    y = "DBH (cm)",
    fill = "Dataset",
    color = "Model"
    #title = "Height vs DBH for selected genera with predicted model lines"
  ) +
  theme_minimal() +
  theme( element_text(size = 18),legend.position = "bottom")

p_dbh_height_all

ggsave(
  filename = paste0(plotdir,"allo_height_dbh_v3_all.png"),  
  plot = p_dbh_height_all,                 
  width = 8,                
  height = 12,               
  dpi = 300                 
)







