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

TaxoWeighter <- function(score_fam, score_gen, score_sp, dfsel)
  {
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_sp <- paste0('MatchVal_Sp_', as.name(score_sp), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  
  df[new_colname_fam] <- ifelse(df[ , "Family"] == df[ , "Family_STDS"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus"] == df[ , "Genus_STDS"], score_gen, 0.0)
  df[new_colname_sp]<- ifelse(df[ , "Species"] == df[ , "Species_STDS"], score_sp, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_fam], df[new_colname_gen], df[new_colname_sp], na.rm = TRUE)
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

TaxoWeighter_UseR <- function(score_fam, score_gen, score_sp, dfsel)
  {
  df <- as.data.frame(dfsel)
  # appends a column name indicating the applied score
  new_colname_fam <- paste0('MatchVal_Fam_', as.name(score_fam), collapse = NULL)
  new_colname_gen <- paste0('MatchVal_Gen_', as.name(score_gen), collapse = NULL)
  new_colname_sp <- paste0('MatchVal_Sp_', as.name(score_sp), collapse = NULL)
  new_colname_pmax <- paste0('Taxonomical_Score_query', collapse = NULL)
  #Scoring Converting to numeric and normalizing
  df <- df %>%
    mutate_at(vars(matches("Score_query")), list(as.numeric)) %>%
    mutate(Normalized_Score_query = (Score_query - min(Score_query)) / (max(Score_query) - min(Score_query)))
  
  df[new_colname_fam] <- ifelse(df[ , "Family_Bio"] == df[ , "Family"], score_fam, 0.0)
  df[new_colname_gen] <- ifelse(df[ , "Genus_Bio"] == df[ , "Genus"], score_gen, 0.0)
  df[new_colname_sp]<- ifelse(df[ , "Species_Bio"] == df[ , "Species"], score_sp, 0.0)
  df[new_colname_pmax] <- pmax(df[new_colname_fam], df[new_colname_gen], df[new_colname_sp], na.rm = TRUE)
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