# By Dr. David V. Conti
library(MultiAssayExperiment)
library(Biobase)
library(tidyverse)

dat_dir = "~/Documents/Research/exposome_data_challenge/data/"
setwd(dat_dir)
load("exposome.RData")
load("proteome.RData")
load("genome.RData")
load("metabol_serum.Rdata")
load("metabol_urine.Rdata")


outdoor.exposures <- exposome[,c("ID", as.character(codebook$variable_name[codebook$domain=="Outdoor exposures"]))] %>% 
  column_to_rownames("ID") %>% 
  t() %>%
  DataFrame()
indoor.air <- exposome[,c("ID", as.character(codebook$variable_name[codebook$domain=="Indoor air"]))] %>% 
  column_to_rownames("ID") %>% 
  t() %>%
  DataFrame()
lifestyles <- exposome[,c("ID", as.character(codebook$variable_name[codebook$domain=="Lifestyles"]))] %>% 
  column_to_rownames("ID") %>% 
  t() %>%
  DataFrame()
chemicals <- exposome[,c("ID", as.character(codebook$variable_name[codebook$domain=="Chemicals"]))] %>% 
  column_to_rownames("ID") %>% 
  t() %>%
  DataFrame()
covariates <- covariates %>% 
  column_to_rownames("ID") %>% 
  t() %>%
  DataFrame()
phenotype <- phenotype %>% as.data.frame() # use as ColData for MultiAssayExperiment format
row.names(phenotype) <- paste0("X", phenotype$ID)

proteome.d <- proteome@assayData$exprs %>% DataFrame()
proteome.cov <- proteome@phenoData@data
proteome.cov <- proteome.cov[stats::complete.cases(proteome.cov),] %>% t() %>% DataFrame()

metabol_urine.d <- metabol_urine@assayData$exprs %>% DataFrame()
metabol_urine.cov <- metabol_urine@phenoData@data
metabol_urine.cov <- metabol_urine.cov[stats::complete.cases(metabol_urine.cov),] %>% t() %>% DataFrame()

metabol_serum.d <- metabol_serum@assayData$exprs %>% DataFrame()
metabol_serum.cov <- metabol_serum@phenoData@data
metabol_serum.cov <- metabol_serum.cov[stats::complete.cases(metabol_serum.cov),] %>% t() %>% DataFrame()


helix_ma <- MultiAssayExperiment(
  experiments= ExperimentList("outdoor.exposures"=outdoor.exposures,
                              "indoor.air"=indoor.air,
                              "lifestyles"=lifestyles,
                              "exposome"=chemicals,
                              "covariates"=covariates,
                              "proteome"=proteome.d,
                              "proteome.cov"=proteome.cov,
                              "metabol_urine"=metabol_urine.d,
                              "metabol_urine.cov"=metabol_urine.cov,
                              "metabol_serum"=metabol_serum.d,
                              "metabol_serum.cov"=metabol_serum.cov,
                              "genome"=G), 
  colData = phenotype)

upsetSamples(helix_ma, nintersects = 10, )