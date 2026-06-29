library(tidyverse)
library(dplyr)
library(ggpubr)
library(biomaRt)
library(cmapR)
library(GSVA)
library(edgeR)
library(uwot)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(readxl)
library(gghalves)

geneExpression <- readRDS("GECounts_0532Pts_02202026.RDS")
geneExpression <- round(geneExpression)

metadata_0532 <- read_xlsx("NEW COG_0532_TMM_classification_Rawdata_6.25.24 - Organized 2_13_25 - Updated 4_30_26 for Merit.xlsx")
metadata_0532_2 <- read_delim("TMM_052_MetaFinal_09162025.txt")
metadata_0532 <- left_join(metadata_0532_2[, c("SampleID", "RNAseq_SampleID")], metadata_0532,
                           by = c("SampleID" = "ID"))


colnames(metadata_0532) <- c("ID","SampleID", "C-CIRCLE CONTENT", "C-CIRCLE ASSAY CLASSIFICATION",
                             "TELOMERE CONTENT", "TERT Expression (AU)",
                             "RIN", "APBs", "HTF Assay", "TMM", "TRAP")
metadata_0532 <- metadata_0532[, c("ID", "SampleID","C-CIRCLE CONTENT", "C-CIRCLE ASSAY CLASSIFICATION",
                                   "TELOMERE CONTENT", "TERT Expression (AU)",
                                   "RIN", "APBs", "HTF Assay", "TMM", "TRAP")]

metadata_0532 <- metadata_0532 %>%
  mutate(TMM = case_when(
    TMM == "ALT+"  ~ "ALT",
    TMM == "TERT+" ~ "Telomerase",
    TMM == "TMM-"  ~ "NO_TMM"
  ))


# filtering for protein coding genes.
#mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 115)
# 
# target_biotypes <- c("protein_coding", "lncRNA", "snRNA", "snoRNA", "miRNA", "rRNA")
# 
# valid_gene_anno <- getBM(
#   attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
#   filters    = "biotype",
#   values     = target_biotypes,
#   mart       = mart
# )

geneExpression <- geneExpression[rownames(geneExpression) %in% valid_gene_anno$hgnc_symbol, , drop = FALSE]


# setting metadata order.
metadata_0532 <- metadata_0532 %>%
  mutate(TMM_Case = case_when(
    TMM == "ALT" ~ "TMM",
    TMM == "Telomerase" ~ "TMM",
    TMM == "NO_TMM" ~ "NO_TMM"
  ))

metadata_0532 <- metadata_0532 %>%
  mutate(ALT_Case = case_when(
    TMM == "ALT" ~ "ALT",
    TMM == "Telomerase" ~ "Non-ALT",
    TMM == "NO_TMM" ~ "Non-ALT"
  ))

# metadata_0532 <- metadata_0532 %>%
#   filter(TMM != "Telomerase" | `TELOMERE CONTENT` < 15)



metadata_0532 <- metadata_0532 %>%
  arrange(ALT_Case, TMM)

metadata_0532 <- metadata_0532 %>%
  mutate(Cohort = "ANBL0532")

metadata_0532 <- metadata_0532 %>%
  mutate(COG.Risk.Group = "High Risk")


colnames(geneExpression) <- sub("_.*", "", colnames(geneExpression))
geneExpression <- geneExpression[, match(metadata_0532$SampleID, colnames(geneExpression)), drop = FALSE]



##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.

metadata_0532$Group <- factor(metadata_0532$TMM, levels = c("ALT", "Telomerase", "NO_TMM"))


# creating differential gene expression object.
dge_TMM <- DGEList(counts=geneExpression,group=metadata_0532$Group)

keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]


# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")

tmm_cpm  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_0532 <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE)

## Filtering extra Telomerase samples.

tert_expr <- tmm_lcpm_0532["TERT", ]

plot_data <- data.frame(TERT_expression = as.numeric(tert_expr),
                        TMM = metadata_0532$TMM, SampleID = metadata_0532$SampleID)

ggplot(plot_data, aes(x = TMM, y = TERT_expression, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic() +  scale_x_discrete(labels = c("ALT" = "ALT+", 
                                                                                         "NO_TMM" = "TMM-ve", 
                                                                                         "Telomerase" = "Telomerase")) +
  theme(legend.position = "none") + scale_fill_manual(values = c("ALT" = "blue", 
                                                                 "NO_TMM" = "green", 
                                                                 "Telomerase" = "darkred")) +
  theme(legend.position = "none")

## filtering Telomerase samples.
# plot_data <- plot_data %>%
#  filter(TMM == "Telomerase" | TERT_expression < 0)

# metadata_0532 <- metadata_0532 %>%
#   filter(SampleID %in% plot_data$SampleID)
# metadata_0532 <- metadata_0532 %>%
#   arrange(ALT_Case, TMM)
# 
# 
# dge_TMM <- dge_TMM[, metadata_0532$SampleID, keep.lib.sizes = FALSE]
# 
# dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")
# 
# 
# tmm_cpm  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
# tmm_lcpm_0532 <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
########


### Using EXTEND to look at Telomerase activity.
extendScores <- RunEXTEND(tmm_lcpm_0532)
extendScores <- read_delim(file = "TelomeraseScores.txt")
extendScores$SampleID <- gsub(".", "-", extendScores$SampleID, fixed = TRUE)
extendScores <- left_join(extendScores, metadata_0532, by = "SampleID")

ggplot(extendScores, aes(x = TMM, y = NormEXTENDScores, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic()

#####################################################################################

design <- model.matrix(~ 0 + Group, data = metadata_0532)
colnames(design) <- levels(metadata_0532$Group)


### Performing differential expression between NO_TMM & TMM.

# building model matrix.
# Calculating dispersion and fitting the model.
d1 <- estimateDisp(dge_TMM, design, verbose=TRUE)
fit1 <- glmQLFit(d1, design, robust = TRUE)

# contrast parameter.
contrast_matrix <- makeContrasts(
  notmmVS_Tel    = NO_TMM - Telomerase,
  notmmVS_alt     = NO_TMM - ALT,
  altVS_Tel   = ALT - Telomerase,
  altVS_NOTMM     = ALT - NO_TMM,
  levels = design
)

# differential expression test.
fitNOTMM_vs_Tel    <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVS_Tel"])
fitNOTMM_vs_ALT     <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVS_alt"])
fitALT_vs_Tel    <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVS_Tel"])
fitALT_vs_NoTMM     <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVS_NOTMM"])

# results.
top_NOTMM_VS_Tel <- topTags(fitNOTMM_vs_Tel, n = Inf)
top_NOTMM_VS_ALT <- topTags(fitNOTMM_vs_ALT, n = Inf)

top_ALT_VS_Tel <- topTags(fitALT_vs_Tel, n = Inf)
top_ALT_VS_NOTMM <- topTags(fitALT_vs_NoTMM, n = Inf)



# filtering for candidate genes.
candidate_genes_NOTMM_vs_Tel <- subset(top_NOTMM_VS_Tel$table, FDR <= 0.01 & abs(logFC) >= 0.5)
candidate_genes_NOTMM_vs_ALT <- subset(top_NOTMM_VS_ALT$table, FDR <= 0.01 & abs(logFC) >= 0.5)


candidate_genes_ALT_vs_Tel <- subset(
  top_ALT_VS_Tel$table, FDR <= 0.01 & abs(logFC) >= 0.5
)
candidate_genes_ALT_vs_NOTMM <- subset(
  top_ALT_VS_NOTMM$table, FDR <= 0.01 & abs(logFC) >= 0.5
)


candidate_genes_tmm_0532 <- Reduce(union, list(
  rownames(candidate_genes_NOTMM_vs_Tel),
  rownames(candidate_genes_NOTMM_vs_ALT)
))

candidate_genes_alt_0532 <- Reduce(union, list(
  rownames(candidate_genes_ALT_vs_Tel),
  rownames(candidate_genes_ALT_vs_NOTMM)
))

#####
#####################################################################################
mycn_expr <- tmm_lcpm_0532["MYCN", colnames(tmm_lcpm_0532) %in% metadata_0532[metadata_0532$TMM == "Telomerase", ]$SampleID]
mycn_expr <- tmm_lcpm_0532["MYCN", ]


plot_data <- data.frame(
  Expression = as.numeric(mycn_expr),
  Phenotype = "Telomerase"
)

plot_data <- data.frame(
  Expression = as.numeric(mycn_expr),
  Phenotype = metadata_0532$TMM
)

ggplot(plot_data, aes(x = Expression)) +
  geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(aes(color = Phenotype), sides = "b", length = unit(0.05, "npc")) + 
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "MYCN Expression Density by Phenotype",
    x = "log2 CPM",
    color = "Phenotype"
  )

#####

metadata_0532$MYCN.status <- ifelse(mycn_expr >= 9, "Amplified",
                                    ifelse(mycn_expr <= 6.875, "Not-Amplified", NA))
metadata_0532 <- metadata_0532 %>%
  filter(TMM != "Telomerase" | (TMM == "Telomerase" & MYCN.status %in% c("Amplified", "Not-Amplified")))

metadata_0532 <- metadata_0532 %>%
  mutate(TMM = case_when(
    TMM == "Telomerase" & MYCN.status == "Amplified" ~ "Telomerase-Amplified",
    TMM == "Telomerase" & MYCN.status == "Not-Amplified" ~ "Telomerase-NotAmplified",
    TRUE ~ TMM 
  ))

metadata_0532 <- metadata_0532 %>%
  arrange(ALT_Case, TMM)
tmm_lcpm_0532 <- tmm_lcpm_0532[, colnames(tmm_lcpm_0532) %in% metadata_0532$SampleID]
tmm_lcpm_0532 <- tmm_lcpm_0532[, match(metadata_0532$SampleID, colnames(tmm_lcpm_0532))]


#############################

metadata0532 <- metadata_0532 ## metadata will be changed quite a bit in training set; metadata_TARGET is final.

###################################################################################
