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

## extend.
source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")

# Loading bulk-seq log counts RNA and metadata file.
Expression <- read_delim("TARGET-NBL.star_counts.tsv")
metadata <- read_excel("TARGETTMMclass_2025_Koneru.xlsx")

# renaming metadata column target.usi to sampleID.
colnames(metadata)[colnames(metadata) == "TARGET.USI"] <- "SampleID"

# matching metadata ID and Expression colnames.
Expression <- Expression %>%
  dplyr::select(Ensembl_ID, ends_with("-01A")) # only keeping sampleIDs ending with 01A.
colnames(Expression) <- sub("^TARGET-30-([^-]+)-.*$", "\\1", colnames(Expression))


# removing Ensembl ID data after . in the name in the counts data.
Expression$Ensembl_ID <- gsub("\\..*", "", Expression$Ensembl_ID)


#mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 115)
# target_biotypes <- c("protein_coding", "lncRNA", "snRNA", "snoRNA", "miRNA", "rRNA")
# valid_gene_anno <- getBM(
#   attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
#   filters    = "biotype",
#   values     = target_biotypes,
#   mart       = mart
# )
# 
# saveRDS(
#   valid_gene_anno,
#   file = "valid_gene_anno_ensembl115.rds"
# )

valid_gene_anno <- readRDS(
  "valid_gene_anno_ensembl115.rds"
)

Expression <- Expression %>%
  filter(Ensembl_ID %in% valid_gene_anno$ensembl_gene_id)


# creating annotation map.
annotation <- valid_gene_anno %>%
  dplyr::rename(ENSEMBL = ensembl_gene_id, SYMBOL = hgnc_symbol) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(ENSEMBL %in% Expression$Ensembl_ID) %>%
  distinct() %>%
  # Removing Ensembl IDs that map to multiple Symbols
  group_by(ENSEMBL) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  # Removing Symbols that map to multiple Ensembl IDs
  group_by(SYMBOL) %>%
  filter(n() == 1) %>%
  ungroup()


Expression <- Expression %>%
  inner_join(annotation, by = c("Ensembl_ID" = "ENSEMBL")) %>%
  mutate(avg_expr = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE)) %>%
  group_by(SYMBOL) %>%
  # Keeping the transcript with the highest average signal
  slice_max(order_by = avg_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-avg_expr, -gene_biotype)

Expression <- as.data.frame(Expression)
rownames(Expression) <- Expression$SYMBOL
Expression <- Expression[, !colnames(Expression) %in% c("Ensembl_ID", "SYMBOL", "avg_expr")]


Expression <- Expression[, colnames(Expression) %in% metadata$SampleID]
metadata <- as.data.frame(metadata[metadata$SampleID %in% colnames(Expression), ])


# adding TMM_Case made of TMM+ and TMM-.
metadata <- metadata %>%
  mutate(TMM_Case = case_when(
    TMM == "ALT+"  ~ "TMM",
    TMM == "TERT+" ~ "TMM",
    TMM == "TMM-"  ~ "NO_TMM",
    TRUE           ~ NA
  ))

metadata <- metadata %>%
  mutate(TMM = case_when(
    TMM == "ALT+"  ~ "ALT",
    TMM == "TERT+" ~ "Telomerase",
    TMM == "TMM-"  ~ "NO_TMM",
    TRUE           ~ NA
  ))

metadata <- metadata %>%
  mutate(ALT_Case = case_when(
    TMM == "ALT" ~ "ALT",
    TRUE ~ "Non-ALT"
  ))

##
## Filtering telomerase samples with TC > 15.
# metadata <- metadata %>%
#   filter(TMM != "Telomerase" | `Telomere Content` < 15)
##

###########################################################################################
##################################################################################################
gene_Expression <- round(2^Expression - 1)


# setting metadata order.
metadata <- metadata %>%
  filter(COG.Risk.Group == "High Risk")

metadata <- metadata %>%
  mutate(Cohort = "TARGET")

batch <- metadata$Cohort
group <- metadata$TMM

metadata <- metadata %>%
  arrange(ALT_Case, TMM)

gene_Expression <- gene_Expression[, colnames(gene_Expression) %in% metadata$SampleID]
gene_Expression <- gene_Expression[, match(metadata$SampleID, colnames(gene_Expression)), drop = FALSE]



### Performing differential expression between Non-ALT & ALT & NO_TMM vs. TMM.

# building model matrix.

# First, determining the factors of TMM.
metadata$Group <- factor(metadata$TMM, levels = c("ALT", "Telomerase", "NO_TMM"))

# creating differential gene expression object.
dge_TMM <- DGEList(counts=gene_Expression, group=metadata$Group)

keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# cpm and log cpm.
cpm_matrix <- cpm(dge_TMM, normalized.lib.sizes = TRUE)
tmm_lcpm_target <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE)




########
## Filtering extra Telomerase samples.

tert_expr <- tmm_lcpm_target["TERT", ]

plot_data <- data.frame(TERT_expression = as.numeric(tert_expr),
                        TMM = metadata$TMM, SampleID = metadata$SampleID)

ggplot(plot_data, aes(x = TMM, y = TERT_expression, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic() +
  scale_x_discrete(labels = c("ALT" = "ALT+", 
                              "NO_TMM" = "TMM-ve", 
                              "Telomerase" = "Telomerase")) +
  theme(legend.position = "none") + scale_fill_manual(values = c("ALT" = "blue", 
                                                                 "NO_TMM" = "green", 
                                                                 "Telomerase" = "darkred"))

## taking the top 90 percentile of Telomerase samples only.
# plot_data <- plot_data %>%
#   filter(TMM == "Telomerase") %>%
#   arrange(-TERT_expression) %>%
#   slice_max(order_by = TERT_expression, prop = 0.90)

# metadata <- metadata %>%
#   filter(TMM != "Telomerase" | SampleID %in% plot_data$SampleID)
metadata <- metadata %>%
  arrange(ALT_Case, TMM)



### Using EXTEND to look at Telomerase activity.
extendScores <- RunEXTEND(tmm_lcpm_target)
extendScores <- read_delim(file = "TelomeraseScores.txt")
extendScores <- left_join(extendScores, metadata, by = "SampleID")

ggplot(extendScores, aes(x = TMM, y = NormEXTENDScores, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic()


# ## Loop back with updated sample size.
# dge_TMM <- dge_TMM[, metadata$SampleID, keep.lib.sizes = FALSE]
# dge_TMM <- dge_TMM[, match(metadata$SampleID, colnames(dge_TMM))]
# 
# 
# # TMM normalization.
# dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")
# 
# # cpm and log cpm.
# cpm_matrix <- cpm(dge_TMM, normalized.lib.sizes = TRUE)
# tmm_lcpm_target <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE)
# 
# tmm_lcpm_target <- tmm_lcpm_target[, colnames(tmm_lcpm_target) %in% metadata$SampleID]
# tmm_lcpm_target <- tmm_lcpm_target[, match(metadata$SampleID, colnames(tmm_lcpm_target))]

design <- model.matrix(~ 0 + Group, data = metadata)
colnames(design) <- levels(metadata$Group)

### for tmm-ve vs. tmm+ve.
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


candidate_genes_tmm_target <- Reduce(union, list(
  rownames(candidate_genes_NOTMM_vs_Tel),
  rownames(candidate_genes_NOTMM_vs_ALT)
))

candidate_genes_alt_target <- Reduce(union, list(
  rownames(candidate_genes_ALT_vs_Tel),
  rownames(candidate_genes_ALT_vs_NOTMM)
))


## changing the TMM style of metadata.
metadata <- metadata%>%
  mutate(TMM = case_when(
    TMM == "Telomerase" & MYCN.status == "Amplified" ~ "Telomerase-Amplified",
    TMM == "Telomerase" & MYCN.status == "Not Amplified" ~ "Telomerase-NotAmplified",
    TRUE ~ TMM 
  ))
##########################################

## Looking at MYCN plots.

mycn_expr <- tmm_lcpm_target["MYCN", ]

plot_data <- data.frame(MYCN_Expression = as.numeric(mycn_expr),
                        Amplification = metadata$MYCN.status,
                        TMM = metadata$TMM,
                        SampleID = metadata$SampleID) 
plot_data <- plot_data %>%
  filter(TMM == "Telomerase-Amplified" | TMM == "Telomerase-NotAmplified")

ggplot(plot_data, aes(x = TMM, y = MYCN_Expression, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic() + scale_fill_manual(values = c("ALT" = "blue", 
                                                                                         "NO_TMM" = "green", 
                                                                                         "Telomerase-Amplified" = "darkred",
                                                                                         "Telomerase-NotAmplified" = "red")) + 
  theme(legend.position = "none") + theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"))

## density plot.
ggplot(plot_data, aes(x = MYCN_Expression)) +
  geom_density(fill = "steelblue", alpha = 0.3) +
  geom_rug(aes(color = TMM), sides = "b", length = unit(0.07, "npc")) + 
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "MYCN Expression Density by Phenotype",
    x = "log2 CPM",
    color = "Phenotype")

####################################################################################

metadata_TARGET <- metadata ## metadata will be changed quite a bit in training set; metadata_TARGET is final.

###################################################################################



