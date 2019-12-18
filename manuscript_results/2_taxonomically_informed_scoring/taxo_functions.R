#######################################################
######   Taxonomical Weighing Functions    ############
#######################################################
###############   Filtering Functions   ###############
#######################################################
#######  Plotting Functions   #########################
#######################################################

library(readr)
library(dplyr)
library(purrr)

TaxoWeighter <- function(
  score_kin,
  score_phy,
  score_cla, 
  score_ord, 
  score_fam,
  score_gen, 
  score_spe, 
  dfsel
  )
  {
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_kin <- paste0('MatchVal_Kin_', as.name(score_kin), collapse = NULL)
  new_colname_phy <- paste0('MatchVal_Phy_', as.name(score_phy), collapse = NULL)
  new_colname_cla <- paste0('MatchVal_Cla_', as.name(score_cla), collapse = NULL)
  new_colname_ord <- paste0('MatchVal_Ord_', as.name(score_ord), collapse = NULL)
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_spe <- paste0('MatchVal_Spe_', as.name(score_spe), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  df[new_colname_kin] <- ifelse(df[ , "Kingdom"] == df[ , "Kingdom_STDS"], score_kin, 0.0)
  df[new_colname_phy] <- ifelse(df[ , "Phylum"] == df[ , "Phylum_STDS"], score_phy, 0.0)
  df[new_colname_cla] <- ifelse(df[ , "Class"] == df[ , "Class_STDS"], score_cla, 0.0)
  df[new_colname_ord] <- ifelse(df[ , "Order"] == df[ , "Order_STDS"], score_ord, 0.0)
  df[new_colname_fam] <- ifelse(df[ , "Family"] == df[ , "Family_STDS"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus"] == df[ , "Genus_STDS"], score_gen, 0.0)
  df[new_colname_spe]<- ifelse(df[ , "Species"] == df[ , "Species_STDS"], score_spe, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_kin], df[new_colname_phy], df[new_colname_cla], df[new_colname_ord], df[new_colname_fam], df[new_colname_gen], df[new_colname_spe], na.rm = TRUE)
  df[new_colname_pmax][is.na(df[new_colname_pmax])] <- 0.0
  
  df <- df %>%
    mutate(Pondered_Score_query = (Taxonomical_Score_query + Normalized_Score_query)) %>%
    group_by(ID_query) %>%
    mutate(RANK_BEFORE = (dense_rank(-Normalized_Score_query))) %>%
    group_by(ID_query) %>%
    mutate(RANK_FINAL = (dense_rank(-Pondered_Score_query))) %>%
    group_by(ID_query) %>%
    arrange(RANK_FINAL,desc(-Pondered_Score_query))
  
  df <- df %>%
    filter(MatchVal == 'Match') %>%
    distinct(ID_query, .keep_all = TRUE)
  
  data.frame(df)
  }

TaxoWeighter_old <- function(
  score_fam,
  score_gen,
  score_spe,
  dfsel
  )
{
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_spe <- paste0('MatchVal_Spe_', as.name(score_spe), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  df[new_colname_fam] <- ifelse(df[ , "Family"] == df[ , "Family_STDS"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus"] == df[ , "Genus_STDS"], score_gen, 0.0)
  df[new_colname_spe]<- ifelse(df[ , "Species"] == df[ , "Species_STDS"], score_spe, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_fam], df[new_colname_gen], df[new_colname_spe], na.rm = TRUE)
  df[new_colname_pmax][is.na(df[new_colname_pmax])] <- 0.0
  
  df <- df %>%
    mutate(Pondered_Score_query = (Taxonomical_Score_query + Normalized_Score_query)) %>%
    group_by(ID_query) %>%
    mutate(RANK_BEFORE = (dense_rank(-Normalized_Score_query))) %>%
    group_by(ID_query) %>%
    mutate(RANK_FINAL = (dense_rank(-Pondered_Score_query))) %>%
    group_by(ID_query) %>%
    arrange(RANK_FINAL,desc(-Pondered_Score_query))
  
  df <- df %>%
    filter(MatchVal == 'Match') %>%
    distinct(ID_query, .keep_all = TRUE)
  
  data.frame(df)
}

TaxoWeighter_UseR <- function(
  score_kin,
  score_phy,
  score_cla,
  score_ord,
  score_fam,
  score_gen,
  score_spe,
  dfsel
  )
{
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_kin <- paste0('MatchVal_Kin_', as.name(score_kin), collapse = NULL)
  new_colname_phy <- paste0('MatchVal_Phy_', as.name(score_phy), collapse = NULL)
  new_colname_cla <- paste0('MatchVal_Cla_', as.name(score_cla), collapse = NULL)
  new_colname_ord <- paste0('MatchVal_Ord_', as.name(score_ord), collapse = NULL)
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_spe <- paste0('MatchVal_Spe_', as.name(score_spe), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  df[new_colname_kin] <- ifelse(df[ , "Kingdom_Bio"] == df[ , "Kingdom"], score_kin, 0.0)
  df[new_colname_phy] <- ifelse(df[ , "Phylum_Bio"] == df[ , "Phylum"], score_phy, 0.0)
  df[new_colname_cla] <- ifelse(df[ , "Class_Bio"] == df[ , "Class"], score_cla, 0.0)
  df[new_colname_ord] <- ifelse(df[ , "Order_Bio"] == df[ , "Order"], score_ord, 0.0)
  df[new_colname_fam] <- ifelse(df[ , "Family_Bio"] == df[ , "Family"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus_Bio"] == df[ , "Genus"], score_gen, 0.0)
  df[new_colname_spe]<- ifelse(df[ , "Species_Bio"] == df[ , "Species"], score_spe, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_kin], df[new_colname_phy], df[new_colname_cla], df[new_colname_ord], df[new_colname_fam], df[new_colname_gen], df[new_colname_spe], na.rm = TRUE)
  df[new_colname_pmax][is.na(df[new_colname_pmax])] <- 0.0
  
  df <- df %>%
    mutate(Pondered_Score_query = (Taxonomical_Score_query + Normalized_Score_query)) %>%
    group_by(ID_query) %>%
    arrange(-Pondered_Score_query) %>%
    distinct(ID_query, Short_IK_query, .keep_all = TRUE) %>%
    mutate(RANK_BEFORE = (dense_rank(-Normalized_Score_query))) %>%
    group_by(ID_query) %>%
    mutate(RANK_FINAL = (dense_rank(-Pondered_Score_query))) %>%
    group_by(ID_query) %>%
    arrange(RANK_FINAL,desc(-Pondered_Score_query))
  
  data.frame(df)
}





###ON GOING

TaxoWeighter_ncbi <- function(
  score_kin,
  score_phy,
  score_cla, 
  score_ord, 
  score_fam,
  score_gen, 
  score_spe, 
  dfsel
)
{
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_kin <- paste0('MatchVal_Kin_', as.name(score_kin), collapse = NULL)
  new_colname_phy <- paste0('MatchVal_Phy_', as.name(score_phy), collapse = NULL)
  new_colname_cla <- paste0('MatchVal_Cla_', as.name(score_cla), collapse = NULL)
  new_colname_ord <- paste0('MatchVal_Ord_', as.name(score_ord), collapse = NULL)
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_spe <- paste0('MatchVal_Spe_', as.name(score_spe), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  df[new_colname_kin] <- ifelse(df[ , "Kingdom"] == df[ , "Kingdom_STDS"], score_kin, 0.0)
  df[new_colname_phy] <- ifelse(df[ , "Phylum"] == df[ , "Phylum_STDS"], score_phy, 0.0)
  df[new_colname_cla] <- ifelse(df[ , "Class"] == df[ , "Class_STDS"], score_cla, 0.0)
  df[new_colname_ord] <- ifelse(df[ , "Order"] == df[ , "Order_STDS"], score_ord, 0.0)
  df[new_colname_fam] <- ifelse(df[ , "Family"] == df[ , "Family_STDS"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus"] == df[ , "Genus_STDS"], score_gen, 0.0)
  df[new_colname_spe]<- ifelse(df[ , "Species"] == df[ , "Species_STDS"], score_spe, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_kin], df[new_colname_phy], df[new_colname_cla], df[new_colname_ord], df[new_colname_fam], df[new_colname_gen], df[new_colname_spe], na.rm = TRUE)
  df[new_colname_pmax][is.na(df[new_colname_pmax])] <- 0.0
  
  df <- df %>%
    mutate(Pondered_Score_query = (Taxonomical_Score_query + Normalized_Score_query)) %>%
    group_by(ID_query) %>%
    mutate(RANK_BEFORE = (dense_rank(-Normalized_Score_query))) %>%
    group_by(ID_query) %>%
    mutate(RANK_FINAL = (dense_rank(-Pondered_Score_query))) %>%
    group_by(ID_query) %>%
    arrange(RANK_FINAL,desc(-Pondered_Score_query))
  
  df <- df %>%
    filter(MatchVal == 'Match') %>%
    distinct(ID_query, .keep_all = TRUE)
  
  data.frame(df)
}