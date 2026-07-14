library(cmapR)
library(tidyverse)
library(biomaRt)
library(dplyr)
library(readxl)
library(ggpubr)
library(edgeR)
library(sva)
library(caret)
library(GSVA)
library(gridExtra)
library(pROC)
library(EnhancedVolcano)
library(gridExtra)
library(hrbrthemes)
library(gghalves)
library(gridExtra)
library(caret)
library(org.Hs.eg.db)
library(limma)
library(randomForest)
library(limma)


Expression_ackerman <- parse_gctx("Neuroblastoma_208Samples.gct")
metadata_ackerman <- read_delim("NBL_Ackerman_CompleteMeta.txt")
metadata_ackerman <- metadata_ackerman %>%
  mutate(RiskGroup = case_when(
    Risk == "YES" ~ "HR",
    TRUE ~ "Non.HR"
  ))


Expression_ackerman <- Expression_ackerman@mat
expr_ackerman <- log2(Expression_ackerman)
####################################################################################

# setting metadata order.
metadata_ackerman <- metadata_ackerman %>%
  mutate(Cohort = "Ackerman")

metadata_ackerman <- metadata_ackerman %>%
  mutate(
    ALT_Case = case_when(
      TMM_Category == "ALT" ~ "ALT",
      TRUE ~ "Non_ALT"
    )
)

expr_ackerman <- expr_ackerman[, colnames(expr_ackerman) %in% metadata_ackerman$SampleID]
metadata_ackerman <- metadata_ackerman[!duplicated(metadata_ackerman$SampleID), ]
expr_ackerman <- expr_ackerman[, match(metadata_ackerman$SampleID, colnames(expr_ackerman)), drop = FALSE]
