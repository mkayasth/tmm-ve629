library(dplyr)
library(tidyverse)
library(ggpubr)
library(edgeR)
library(ggfortify)
library(sva)
library(caret)
library(GSVA)
library(gridExtra)
library(pROC)
library(EnhancedVolcano)
library(readxl)
library(org.Hs.eg.db)
library(randomForest)
library(ggtext)
library(GSVA)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)


source("tmm-ve629/targetDataCleaning.R")
source("tmm-ve629/0532DataCleaning.R")


common_genes <- intersect(rownames(gene_Expression), rownames(geneExpression))


counts_combined <- cbind(geneExpression[common_genes, ],
                         gene_Expression[common_genes, ])


# making metadata.
metadata_target <- metadata_TARGET[, c("SampleID", "TMM", "TMM_Case", "ALT_Case", "COG.Risk.Group", "Cohort", "Telomere Content")]

metadata_0532_2 <- metadata0532[, c("SampleID", "TMM", "TMM_Case", "ALT_Case", "COG.Risk.Group", "Cohort", "TELOMERE CONTENT")]
colnames(metadata_0532_2)[grep("TELOMERE", colnames(metadata_0532_2))] <- "Telomere Content"

metadata_combined <- rbind(metadata_target, metadata_0532_2)
metadata_combined <- metadata_combined %>%
  arrange(ALT_Case, TMM, Cohort, COG.Risk.Group, `Telomere Content`)

batch <- metadata_combined$Cohort
group <- metadata_combined$TMM

counts_combined <- counts_combined[, match(metadata_combined$SampleID, colnames(counts_combined))]

## running combat-seq to merge the two raw datasets.
adjusted_counts <- ComBat_seq(
  counts = as.matrix(counts_combined),
  batch = batch,
  group = group
)

##################################################################################
##################################################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata_combined$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

# creating differential gene expression object.
dge_TMM <- DGEList(counts=adjusted_counts, samples = metadata_combined, group=group1)
keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes=FALSE]

dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")

# cpm and log cpm.
cpm_matrix <- cpm(dge_TMM, normalized.lib.sizes = TRUE)
lcpm_all <- cpm(dge_TMM, log = TRUE, prior.count = 2)


####################################################################################
####################################################################################

# Now, changing metadata and lcpm such that they are in metadata_TARGET or metadata0532.
metadata_combined <- metadata_combined %>%
  filter(SampleID %in% metadata_TARGET$SampleID | SampleID %in% metadata0532$SampleID)
metadata_combined <- metadata_combined %>%
  arrange(ALT_Case, TMM, Cohort, COG.Risk.Group)

lookup_mycn <- bind_rows(
  metadata_TARGET %>% dplyr::select(SampleID, MYCN.status),
  metadata0532   %>% dplyr::select(SampleID, MYCN.status)
) %>% distinct(SampleID, .keep_all = TRUE)

metadata_combined <- metadata_combined %>%
  left_join(lookup_mycn, by = "SampleID")

metadata_combined <- metadata_combined %>%
  mutate(MYCN.status = case_when(
    MYCN.status == "Not Amplified" ~ "Not-Amplified",
    TRUE ~ MYCN.status 
  ))

common_genes <- intersect(rownames(tmm_lcpm_target), rownames(tmm_lcpm_0532))
common_genes <- intersect(rownames(dge_TMM), common_genes)

dge_TMM <- dge_TMM[common_genes, metadata_combined$SampleID, keep.lib.sizes = FALSE]
keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes=FALSE]
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# cpm and log cpm.
cpm_matrix <- cpm(dge_TMM, normalized.lib.sizes = TRUE)
lcpm_all <- cpm(dge_TMM, log = TRUE, prior.count = 2)

lcpm_all <- lcpm_all[, match(metadata_combined$SampleID, colnames(lcpm_all))]

###
# checking if batch correction worked.
group1 <- as.factor(metadata_combined$TMM)

counts_combined <- counts_combined[common_genes, metadata_combined$SampleID]
dge_TMM_raw <- DGEList(counts=counts_combined, group=group1)

# TMM normalization.
dge_TMM_raw <- calcNormFactors(dge_TMM_raw, method = "TMM")

# cpm and log cpm.
cpm_matrix_raw <- cpm(dge_TMM_raw, normalized.lib.sizes = TRUE)
lcpm_raw <- cpm(dge_TMM_raw, log = TRUE, prior.count = 2)


pca_pre  <- prcomp(t(lcpm_raw))
pca_post <- prcomp(t(lcpm_all))

# Plot PCA colored by batch
autoplot(pca_pre,  data = metadata_combined, colour = 'Cohort') +
  ggtitle("Before Batch Correction") +
  theme_classic(base_size = 18)

autoplot(pca_post, data = metadata_combined, colour = 'Cohort') +
  ggtitle("After Batch Correction") +
  theme_classic(base_size = 18)

autoplot(pca_post, data = metadata_combined, colour = 'TMM') +
  ggtitle("After Batch Correction") +
  theme_classic(base_size = 18)

#####
pca_data <- as.data.frame(pca_pre$x)
pca_data$SampleID <- rownames(pca_data)

# Joining with metadata
pca_plot_df <- left_join(pca_data, metadata_combined, by = "SampleID")

# Variance Explained for axis labels
pc1_var <- round(summary(pca_pre)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_pre)$importance[2,2] * 100, 1)
# 
# Plot manually
ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Cohort, shape = TMM)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "Before Batch Correction") +
  theme_classic(base_size = 18)


## Post pca.
pca_data <- as.data.frame(pca_post$x)
pca_data$SampleID <- rownames(pca_data)

# Joining with metadata
pca_plot_df <- left_join(pca_data, metadata_combined, by = "SampleID")

# Variance Explained for axis labels
pc1_var <- round(summary(pca_pre)$importance[2,1] * 100, 1)
pc2_var <- round(summary(pca_pre)$importance[2,2] * 100, 1)
# 
# Plot manually
ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Cohort, shape = TMM)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "After Batch Correction") +
  theme_classic(base_size = 18)

##################################################################################
metadata_combined <- metadata_combined %>%
  mutate(TMM = case_when(
    SampleID %in% metadata_TARGET$SampleID & TMM == "Telomerase" ~ 
      metadata_TARGET$TMM[match(SampleID, metadata_TARGET$SampleID)],
    
    SampleID %in% metadata0532$SampleID & TMM == "Telomerase" ~ 
      metadata0532$TMM[match(SampleID, metadata0532$SampleID)],
    
    TRUE ~ TMM
  ))


###################################################################################
###################################################################################

### Data split.
Strata <- interaction(metadata_combined$Cohort, metadata_combined$TMM)

set.seed(17) # For reproducibility83, 41.
train_indices <- createDataPartition(Strata, p = 0.60, list = FALSE)

# training and testing set.
train_metadata <- metadata_combined[train_indices, ]
test_metadata  <- metadata_combined[-train_indices, ]

lcpm <- lcpm_all[common_genes, colnames(lcpm_all) %in% train_metadata$SampleID]
lcpm_test <- lcpm_all[common_genes, colnames(lcpm_all) %in% test_metadata$SampleID]


# edgeR -- for training set.
train_metadata <- train_metadata %>%
  arrange(ALT_Case, TMM)
adjusted_counts_train <- adjusted_counts[, colnames(adjusted_counts) %in% 
                                           train_metadata$SampleID]

adjusted_counts_train <- adjusted_counts_train[, match(train_metadata$SampleID, 
                                                       colnames(adjusted_counts_train))] 

train_metadata <- train_metadata %>%
  mutate(
    Group = case_when(
      TMM %in% c("Telomerase-Amplified", "Telomerase-NotAmplified") ~ "Telomerase",
      TRUE ~ TMM
    )
  )

train_metadata$Group <- factor(train_metadata$Group, levels = c("ALT", "Telomerase", "NO_TMM"))
train_metadata$TMM <- factor(train_metadata$TMM, levels = c("ALT", "Telomerase-Amplified", "Telomerase-NotAmplified", "NO_TMM"))


dge_TMM <- DGEList(counts= adjusted_counts_train, group= train_metadata$TMM)
keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


design <- model.matrix(~ 0 + TMM, data = train_metadata)
colnames(design) <- levels(train_metadata$TMM)
colnames(design) <- make.names(colnames(design))

### for tmm-ve vs. tmm+ve and ALT vs ALT-ve.
# Calculating dispersion and fitting the model.
d1 <- estimateDisp(dge_TMM, design, verbose=TRUE)
fit1 <- glmQLFit(d1, design, robust = TRUE)


# contrast parameter.
contrast_matrix <- makeContrasts(
  notmmVS_TelAmp    = NO_TMM - Telomerase.Amplified,
  notmmVS_TelNotAmp = NO_TMM - Telomerase.NotAmplified,
  notmmVS_ALT     = NO_TMM - ALT,
  
  altVS_TelAmp    = ALT - Telomerase.Amplified,
  altVS_TelNotAmp = ALT - Telomerase.NotAmplified,
  altVS_notmm    = ALT - NO_TMM,
  levels = design
)

# differential expression test.
fitnotmm_vs_TelAmp    <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVS_TelAmp"])
fitnotmm_vs_TelNotAmp <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVS_TelNotAmp"])
fitnotmm_vs_ALT     <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVS_ALT"])

fitALT_vs_TelAmp    <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVS_TelAmp"])
fitALT_vs_TelNotAmp <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVS_TelNotAmp"])
fitALT_vs_notmm     <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVS_notmm"])

# results.
top_notmm_VS_TelAmp_combined <- topTags(fitnotmm_vs_TelAmp, n = Inf)
top_notmm_VS_TelNotAmp_combined <- topTags(fitnotmm_vs_TelNotAmp, n = Inf)
top_notmm_VS_ALT_combined <- topTags(fitnotmm_vs_ALT, n = Inf)

top_ALT_VS_TelAmp_combined <- topTags(fitALT_vs_TelAmp, n = Inf)
top_ALT_VS_TelNotAmp_combined <- topTags(fitALT_vs_TelNotAmp, n = Inf)
top_ALT_VS_notmm_combined <- topTags(fitALT_vs_notmm, n = Inf)


# filtering for candidate genes.
candidate_genes_notmm_vs_TelAmp_combined <- subset(
  top_notmm_VS_TelAmp_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)
candidate_genes_notmm_vs_TelNotAmp_combined <- subset(
  top_notmm_VS_TelNotAmp_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)
candidate_genes_notmm_vs_ALT_combined <- subset(
  top_notmm_VS_ALT_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)

candidate_genes_ALT_vs_TelAmp_combined <- subset(
  top_ALT_VS_TelAmp_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)
candidate_genes_ALT_vs_TelNotAmp_combined <- subset(
  top_ALT_VS_TelNotAmp_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)
candidate_genes_ALT_vs_notmm_combined <- subset(
  top_ALT_VS_notmm_combined$table, FDR <= 0.01 & abs(logFC) >= 0.5
)

###############
# ALT pairwise - extract up/down per contrast

alt_up_TelAmp    <- rownames(candidate_genes_ALT_vs_TelAmp_combined[
  candidate_genes_ALT_vs_TelAmp_combined$logFC > 0 &
    rownames(candidate_genes_ALT_vs_TelAmp_combined) %in% rownames(lcpm), ])

alt_up_TelNotAmp <- rownames(candidate_genes_ALT_vs_TelNotAmp_combined[
  candidate_genes_ALT_vs_TelNotAmp_combined$logFC > 0 &
    rownames(candidate_genes_ALT_vs_TelNotAmp_combined) %in% rownames(lcpm), ])

alt_up_NoTMM     <- rownames(candidate_genes_ALT_vs_notmm_combined[
  candidate_genes_ALT_vs_notmm_combined$logFC > 0 &
    rownames(candidate_genes_ALT_vs_notmm_combined) %in% rownames(lcpm), ])

alt_down_TelAmp    <- rownames(candidate_genes_ALT_vs_TelAmp_combined[
  candidate_genes_ALT_vs_TelAmp_combined$logFC < 0 &
    rownames(candidate_genes_ALT_vs_TelAmp_combined) %in% rownames(lcpm), ])

alt_down_TelNotAmp <- rownames(candidate_genes_ALT_vs_TelNotAmp_combined[
  candidate_genes_ALT_vs_TelNotAmp_combined$logFC < 0 &
    rownames(candidate_genes_ALT_vs_TelNotAmp_combined) %in% rownames(lcpm), ])

alt_down_NoTMM     <- rownames(candidate_genes_ALT_vs_notmm_combined[
  candidate_genes_ALT_vs_notmm_combined$logFC < 0 &
    rownames(candidate_genes_ALT_vs_notmm_combined) %in% rownames(lcpm), ])

# Union then remove conflicting
alt_up_union   <- Reduce(union, list(alt_up_TelAmp,   alt_up_TelNotAmp,   alt_up_NoTMM))
alt_down_union <- Reduce(union, list(alt_down_TelAmp, alt_down_TelNotAmp, alt_down_NoTMM))

alt_conflicting <- intersect(alt_up_union, alt_down_union)

candidate_genes_alt_combined_up   <- setdiff(alt_up_union,   alt_conflicting)
candidate_genes_alt_combined_down <- setdiff(alt_down_union, alt_conflicting)


###################################################################################
# NO_TMM pairwise - extract up/down per contrast

notmm_up_TelAmp    <- rownames(candidate_genes_notmm_vs_TelAmp_combined[
  candidate_genes_notmm_vs_TelAmp_combined$logFC > 0 &
    rownames(candidate_genes_notmm_vs_TelAmp_combined) %in% rownames(lcpm), ])

notmm_up_TelNotAmp <- rownames(candidate_genes_notmm_vs_TelNotAmp_combined[
  candidate_genes_notmm_vs_TelNotAmp_combined$logFC > 0 &
    rownames(candidate_genes_notmm_vs_TelNotAmp_combined) %in% rownames(lcpm), ])

notmm_up_ALT       <- rownames(candidate_genes_notmm_vs_ALT_combined[
  candidate_genes_notmm_vs_ALT_combined$logFC > 0 &
    rownames(candidate_genes_notmm_vs_ALT_combined) %in% rownames(lcpm), ])

notmm_down_TelAmp    <- rownames(candidate_genes_notmm_vs_TelAmp_combined[
  candidate_genes_notmm_vs_TelAmp_combined$logFC < 0 &
    rownames(candidate_genes_notmm_vs_TelAmp_combined) %in% rownames(lcpm), ])

notmm_down_TelNotAmp <- rownames(candidate_genes_notmm_vs_TelNotAmp_combined[
  candidate_genes_notmm_vs_TelNotAmp_combined$logFC < 0 &
    rownames(candidate_genes_notmm_vs_TelNotAmp_combined) %in% rownames(lcpm), ])

notmm_down_ALT       <- rownames(candidate_genes_notmm_vs_ALT_combined[
  candidate_genes_notmm_vs_ALT_combined$logFC < 0 &
    rownames(candidate_genes_notmm_vs_ALT_combined) %in% rownames(lcpm), ])

# Union then remove conflicting
notmm_up_union   <- Reduce(union, list(notmm_up_TelAmp,   notmm_up_TelNotAmp,   notmm_up_ALT))
notmm_down_union <- Reduce(union, list(notmm_down_TelAmp, notmm_down_TelNotAmp, notmm_down_ALT))

notmm_conflicting <- intersect(notmm_up_union, notmm_down_union)

candidate_genes_tmm_combined_up   <- setdiff(notmm_up_union,   notmm_conflicting)
candidate_genes_tmm_combined_down <- setdiff(notmm_down_union, notmm_conflicting)


###############




###########
gene_contrast_counts_alt <- table(c(
  rownames(candidate_genes_ALT_vs_TelAmp_combined),
  rownames(candidate_genes_ALT_vs_TelNotAmp_combined),
  rownames(candidate_genes_ALT_vs_notmm_combined)
)) %>% as.data.frame() %>%
  dplyr::rename(Gene = Var1) %>%
  filter(Gene %in% c(candidate_genes_alt_combined_up, candidate_genes_alt_combined_down)) %>%
  arrange(desc(Freq))

gene_contrast_counts_tmm <- table(c(
  rownames(candidate_genes_notmm_vs_TelAmp_combined),
  rownames(candidate_genes_notmm_vs_TelNotAmp_combined),
  rownames(candidate_genes_notmm_vs_ALT_combined)
)) %>%
  as.data.frame() %>%
  dplyr::rename(Gene = Var1) %>%
  filter(Gene %in% c(candidate_genes_tmm_combined_up, candidate_genes_tmm_combined_down)) %>%
  arrange(desc(Freq))


#################################################
candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[
  candidate_genes_alt_combined_up %in% 
    gene_contrast_counts_alt$Gene[gene_contrast_counts_alt$Freq >= 2]
]

candidate_genes_alt_combined_down <- candidate_genes_alt_combined_down[
  candidate_genes_alt_combined_down %in% 
    gene_contrast_counts_alt$Gene[gene_contrast_counts_alt$Freq >= 2]
]

# same for TMM
candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[
  candidate_genes_tmm_combined_up %in% 
    gene_contrast_counts_tmm$Gene[gene_contrast_counts_tmm$Freq >= 2]
]

candidate_genes_tmm_combined_down <- candidate_genes_tmm_combined_down[
  candidate_genes_tmm_combined_down %in% 
    gene_contrast_counts_tmm$Gene[gene_contrast_counts_tmm$Freq >= 2]
]

###################################################################################
# Seeing where the genes are coming from.

tmm_up_membership <- data.frame(
  Gene = candidate_genes_tmm_combined_up,
  
  notmm_vs_TelAmp =
    candidate_genes_tmm_combined_up %in% notmm_up_TelAmp,
  
  notmm_vs_TelNotAmp =
    candidate_genes_tmm_combined_up %in% notmm_up_TelNotAmp,
  
  notmm_vs_ALT =
    candidate_genes_tmm_combined_up %in% notmm_up_ALT
)

tmm_up_membership$Total <- rowSums(tmm_up_membership[, -1])

tmm_up_membership <- tmm_up_membership %>%
  arrange(desc(Total))

##########################################

tmm_down_membership <- data.frame(
  Gene = candidate_genes_tmm_combined_down,
  
  notmm_vs_TelAmp =
    candidate_genes_tmm_combined_down %in% notmm_down_TelAmp,
  
  notmm_vs_TelNotAmp =
    candidate_genes_tmm_combined_down %in% notmm_down_TelNotAmp,
  
  notmm_vs_ALT =
    candidate_genes_tmm_combined_down %in% notmm_down_ALT
)

tmm_down_membership$Total <- rowSums(tmm_down_membership[, -1])

tmm_down_membership <- tmm_down_membership %>%
  arrange(desc(Total))


#####
alt_up_membership <- data.frame(
  Gene = candidate_genes_alt_combined_up,
  
  alt_vs_TelAmp =
    candidate_genes_alt_combined_up %in% alt_up_TelAmp,
  
  alt_vs_TelNotAmp =
    candidate_genes_alt_combined_up %in% alt_up_TelNotAmp,
  
  alt_vs_NoTMM =
    candidate_genes_alt_combined_up %in% alt_up_NoTMM
)

alt_up_membership$Total <- rowSums(alt_up_membership[, -1])

alt_up_membership <- alt_up_membership %>%
  arrange(desc(Total))



alt_down_membership <- data.frame(
  Gene = candidate_genes_alt_combined_down,
  
  alt_vs_TelAmp =
    candidate_genes_alt_combined_down %in% alt_down_TelAmp,
  
  alt_vs_TelNotAmp =
    candidate_genes_alt_combined_down %in% alt_down_TelNotAmp,
  
  alt_vs_NoTMM =
    candidate_genes_alt_combined_down %in% alt_down_NoTMM
)

alt_down_membership$Total <- rowSums(alt_down_membership[, -1])

alt_down_membership <- alt_down_membership %>%
  arrange(desc(Total))




###################################################################################
candidate_genes_alt_combined <- c(candidate_genes_alt_combined_up, candidate_genes_alt_combined_down)

candidate_genes_tmm_combined <- c(candidate_genes_tmm_combined_up, candidate_genes_tmm_combined_down)

###############################################################################################

## MYCN check.
mycn_expr <- lcpm["MYCN", ]

plot_data <- data.frame(MYCN_Expression = as.numeric(mycn_expr),
                        TMM = train_metadata$TMM,
                        SampleID = train_metadata$SampleID)
plot_data <- plot_data %>%
  filter(TMM == "Telomerase-Amplified" | TMM == "Telomerase-NotAmplified")

ggplot(plot_data, aes(x = TMM, y = MYCN_Expression, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic()

## TERT check.

tert_expr <- lcpm["TERT", ]

plot_data <- data.frame(TERT_Expression = as.numeric(tert_expr),
                        TMM = train_metadata$TMM,
                        SampleID = train_metadata$SampleID)


ggplot(plot_data, aes(x = TMM, y = TERT_Expression, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic()

## EXTEND check.
extendScores <- RunEXTEND(lcpm)
extendScores <- read_delim(file = "TelomeraseScores.txt")

extendScores$SampleID <- gsub(".", "-", extendScores$SampleID, fixed = TRUE)
extendScores <- left_join(extendScores, train_metadata, by = "SampleID")


ggplot(extendScores, aes(x = TMM, y = NormEXTENDScores, fill = TMM)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + theme_classic()


####################################################################################
# ## including RRHO2 filtered genes.
# # source("tmm-ve/rrhoTrainingSet.R")
# 
# candidate_genes_tmm_combined_up <- rownames(candidate_genes_TMM_combined[candidate_genes_TMM_combined$logFC > 0 & 
#                                                                            rownames(candidate_genes_TMM_combined) %in% rownames(lcpm), ])
# candidate_genes_tmm_combined_down <- rownames(candidate_genes_TMM_combined[candidate_genes_TMM_combined$logFC < 0 & 
#                                                                              rownames(candidate_genes_TMM_combined) %in% rownames(lcpm), ])
# 
# candidate_genes_alt_combined_up <- rownames(candidate_genes_ALT_combined[candidate_genes_ALT_combined$logFC > 0 & 
#                                                                            rownames(candidate_genes_ALT_combined) %in% rownames(lcpm), ])
# candidate_genes_alt_combined_down <- rownames(candidate_genes_ALT_combined[candidate_genes_ALT_combined$logFC < 0 & 
#                                                                            rownames(candidate_genes_ALT_combined) %in% rownames(lcpm), ])
# 
# 
# # candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[candidate_genes_tmm_combined_up %in% up_up_genes_target0532_tmm]
# # candidate_genes_tmm_combined_down <- candidate_genes_tmm_combined_down[candidate_genes_tmm_combined_down %in% down_down_genes_target0532_tmm]
# # 
# # candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[candidate_genes_alt_combined_up %in% up_up_genes_target0532_alt]
# # candidate_genes_alt_combined_down <- candidate_genes_alt_combined_down[candidate_genes_alt_combined_down %in% down_down_genes_target0532_alt]
# 
# ###
# candidate_genes_tmm_combined <- c(candidate_genes_tmm_combined_up, candidate_genes_tmm_combined_down)
# candidate_genes_alt_combined <- c(candidate_genes_alt_combined_up,  candidate_genes_alt_combined_down)
#####################################################################################

###################################################################################


# Step 2: Random Forest: first for TMM-ve.


# Running RF multiple times for stability.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 12,
  "ALT"    = 6,
  "Telomerase-Amplified" = 6,
  "Telomerase-NotAmplified" = 6
)



for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  rf <- randomForest(
    x = t(lcpm[candidate_genes_tmm_combined, ]),
    y = as.factor(train_metadata$TMM_Case),
    sampsize   = sample_sizes,
    strata     = train_metadata$TMM,
    ntree = 1500,
    importance = TRUE
  )
  imp <- as.data.frame(importance(rf))
  imp$gene <- rownames(imp)
  imp_seed_id <- i
  imp$rank_in_accuracy <- rank(-imp$MeanDecreaseAccuracy)
  imp$rank_in_Gini <- rank(-imp$MeanDecreaseGini)
  importance_list[[i]] <- imp
}


# Data Fusion & Consistency Calculation
all_seed_data <- bind_rows(importance_list)

# Summarize metrics across all seeds
stability_metrics <- all_seed_data %>%
  group_by(gene) %>%
  summarise(
    # Core Averages
    mean_accuracy = mean(MeanDecreaseAccuracy),
    mean_gini     = mean(MeanDecreaseGini),
    mean_no_tmm   = mean(NO_TMM),
    mean_tmm      = mean(TMM),
    
    # Stability/Consistency Metrics.
    seeds_present_in_top = sum(rank_in_accuracy <= top_n_threshold & rank_in_Gini <= top_n_threshold ),
    avg_accuracy_rank_in_seeds    = mean(rank_in_accuracy),
    avg_rank_in_Gini_in_seeds = mean(rank_in_Gini)
  ) %>%
  # Quality filter: must have positive importance and valid contribution.
  filter(mean_no_tmm > 0, mean_tmm > 0, mean_accuracy > quantile(mean_accuracy, 0.50)) %>%
  arrange(desc(mean_accuracy)) %>%
  mutate(gene_index = row_number())


# Visualization: Stability vs. Magnitude
p2 <- ggplot(stability_metrics, aes(x = seeds_present_in_top, y = mean_accuracy)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#2E86AB") +
  geom_text_repel(
    data = subset(stability_metrics, seeds_present_in_top >= 8),
    aes(label = gene),
    size = 4.5,
    max.overlaps = 20
  ) +
  labs(
    title = "Consistency Selection (Top 20 across Seeds)",
    x = paste0(
      "Times in Top ", top_n_threshold,
      " (out of ", length(seeds), " seeds)"
    ),
    y = "Average Mean Decrease Accuracy"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

# Top 20 Lollipop Plot
top_20_tmm <- stability_metrics %>% slice_head(n = 20)

p3 <- ggplot(top_20_tmm, aes(x = reorder(gene, mean_accuracy), y = mean_accuracy)) +
  geom_segment(aes(xend = gene, y = 0, yend = mean_accuracy), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Mean Decrease Accuracy") +
  theme_classic()

print(p2)
print(p3)







########################################################################################

####################################################################################
candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[candidate_genes_tmm_combined_up %in% stability_metrics$gene]
stability_metrics[stability_metrics$gene %in% candidate_genes_tmm_combined_up, ]

candidate_genes_tmm_combined_down <- candidate_genes_tmm_combined_down[candidate_genes_tmm_combined_down %in% stability_metrics$gene]
stability_metrics[stability_metrics$gene %in% candidate_genes_tmm_combined_down, ]

tmm_up_membership <- tmm_up_membership[tmm_up_membership$Gene %in% candidate_genes_tmm_combined_up, ]
tmm_down_membership <- tmm_down_membership[tmm_down_membership$Gene %in% candidate_genes_tmm_combined_down, ]
###################################################################################

# Step 2: Random Forest: first for ALT.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 6,
  "ALT"    = 12,
  "Telomerase-Amplified" = 6,
  "Telomerase-NotAmplified" = 6
)



for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  rf <- randomForest(
    x = t(lcpm[candidate_genes_alt_combined, ]),
    y = as.factor(train_metadata$ALT_Case),
    sampsize   = sample_sizes,
    strata     = train_metadata$TMM,
    ntree = 1500,
    importance = TRUE
  )
  imp <- as.data.frame(importance(rf))
  imp$gene <- rownames(imp)
  imp_seed_id <- i
  imp$rank_in_accuracy <- rank(-imp$MeanDecreaseAccuracy)
  imp$rank_in_Gini <- rank(-imp$MeanDecreaseGini)
  importance_list[[i]] <- imp
}


# Data Fusion & Consistency Calculation
all_seed_data2 <- bind_rows(importance_list)

# Summarize metrics across all seeds
stability_metrics2 <- all_seed_data2 %>%
  group_by(gene) %>%
  summarise(
    # Core Averages
    mean_accuracy = mean(MeanDecreaseAccuracy),
    mean_gini     = mean(MeanDecreaseGini),
    mean_alt   = mean(ALT),
    mean_nonalt      = mean(`Non-ALT`),
    
    # Stability/Consistency Metrics.
    seeds_present_in_top = sum(rank_in_accuracy <= top_n_threshold & rank_in_Gini <= top_n_threshold),
    avg_accuracy_rank_in_seeds = mean(rank_in_accuracy),
    avg_Gini_rank_in_seeds = mean(rank_in_Gini)
  ) %>%
  # Quality filter: must have positive importance and valid contribution.
  filter(mean_alt > 0, mean_nonalt > 0, mean_accuracy > quantile(mean_accuracy, 0.50)) %>%
  arrange(desc(mean_accuracy)) %>%
  mutate(gene_index = row_number())


# Visualization: Stability vs. Magnitude
p5 <- ggplot(stability_metrics2, aes(x = seeds_present_in_top, y = mean_accuracy)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#2E86AB") +
  geom_text_repel(
    data = subset(stability_metrics2, seeds_present_in_top >= 10),
    aes(label = gene),
    size = 4.5,
    max.overlaps = 20
  ) +
  labs(
    title = "Consistency Selection (Top 20 across Seeds)",
    x = paste0(
      "Times in Top ", top_n_threshold,
      " (out of ", length(seeds), " seeds)"
    ),
    y = "Average Mean Decrease Accuracy"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

# Top 20 Lollipop Plot
top_20_alt <- stability_metrics2 %>% slice_head(n = 20)

p6 <- ggplot(top_20_alt, aes(x = reorder(gene, mean_accuracy), y = mean_accuracy)) +
  geom_segment(aes(xend = gene, y = 0, yend = mean_accuracy), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Mean Decrease Accuracy") +
  theme_classic()

print(p5)
print(p6)

########################################################################################

candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[candidate_genes_alt_combined_up %in% stability_metrics2$gene]
stability_metrics2[stability_metrics2$gene %in% candidate_genes_alt_combined_up, ]

candidate_genes_alt_combined_down <- candidate_genes_alt_combined_down[candidate_genes_alt_combined_down %in% stability_metrics2$gene]
stability_metrics2[stability_metrics2$gene %in% candidate_genes_alt_combined_down, ]

###################################################################################

alt_up_membership <- alt_up_membership[alt_up_membership$Gene %in% candidate_genes_alt_combined_up, ]
alt_down_membership <- alt_down_membership[alt_down_membership$Gene %in% candidate_genes_alt_combined_down, ]
###################################################################################
### for heatmap, making extra strats.
train_metadata$Plot_Group <- as.character(train_metadata$TMM)


# separating TMM patients with better survival.
train_metadata$Vital_Status <- metadata_TARGET$Vital.Status[
  match(train_metadata$SampleID, metadata_TARGET$SampleID)]
train_metadata$OS_days <- metadata_TARGET$Overall.Survival.Time.in.Days[
  match(train_metadata$SampleID, metadata_TARGET$SampleID)]
train_metadata$OS_event <- ifelse(train_metadata$Vital_Status == "Dead", 1, 0)

# Split each TMM group by optimal survival cutpoint.
tmm_groups_with_surv <- c("ALT", "Telomerase-Amplified", "Telomerase-NotAmplified")

surv_group_labels <- rep(NA_character_, nrow(train_metadata))

for (grp in tmm_groups_with_surv) {
  
  idx <- which(
    train_metadata$Cohort == "TARGET" &
      train_metadata$TMM == grp &
      !is.na(train_metadata$OS_days) &
      !is.na(train_metadata$OS_event)
  )
  
  if (length(idx) < 10) {
    surv_group_labels[idx] <- grp
    next
  }
  
  sub_meta <- train_metadata[idx, ]
  sub_meta$OS_days_var <- sub_meta$OS_days
  
  cut_res <- surv_cutpoint(
    sub_meta,
    time      = "OS_days",
    event     = "OS_event",
    variables = "OS_days_var"
  )
  cut_cat <- surv_categorize(cut_res)
  
  # cutpoint value
  cp_value <- cut_res$cutpoint$cutpoint
  
  labels_grp <- ifelse(
    cut_cat$OS_days_var == "high",
    paste0(grp, "_LongSurv"),
    paste0(grp, "_ShortSurv")
  )
  
  # Alive patients below the cutpoint are late-entry/short follow-up,
  # not true short survivors — reassign to LongSurv
  early_censored <- which(
    sub_meta$OS_event == 0 &        # alive (censored)
      sub_meta$OS_days < cp_value   # below the cutpoint
  )
  labels_grp[early_censored] <- paste0(grp, "_LongSurv")
  
  surv_group_labels[idx] <- labels_grp
}

# Applying survival labels for TARGET TMM groups
target_tmm_idx <- which(
  train_metadata$Cohort == "TARGET" &
    train_metadata$TMM %in% tmm_groups_with_surv &
    !is.na(surv_group_labels)
)
train_metadata$Plot_Group[target_tmm_idx] <- surv_group_labels[target_tmm_idx]

# NO_TMM high telomere content.
train_metadata$Plot_Group[
  train_metadata$TMM == "NO_TMM" & train_metadata$`Telomere Content` > 15
] <- "NO_TMM_HighTC"

train_metadata$Plot_Group[
  train_metadata$TMM == "Telomerase-Amplified" & train_metadata$`Telomere Content` > 15
] <- "TA_HighTC"

train_metadata$Plot_Group[
  train_metadata$TMM == "Telomerase-NotAmplified" & train_metadata$`Telomere Content` > 15
] <- "TNA_HighTC"


###################


## Heatmap of TMM-ve upregulated genes in training set.

selected_genes <- c(candidate_genes_tmm_combined_up, candidate_genes_tmm_combined_down)
row_split_vector <- c(rep("Up-regulated", length(candidate_genes_tmm_combined_up)), 
                      rep("Down-regulated", length(candidate_genes_tmm_combined_down)))


plot_matrix        <- lcpm[selected_genes, ]
plot_matrix_scaled <- t(scale(t(plot_matrix)))
col_fun            <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))

group_colors <- c(
  # TARGET — with survival splits
  "ALT_LongSurv"                        = "#1A5C9E",
  "ALT_ShortSurv"                       = "#7EB6E8",
  "Telomerase-Amplified_LongSurv"       = "#9C1A1A",
  "Telomerase-Amplified_ShortSurv"      = "#E07070",
  "Telomerase-NotAmplified_LongSurv"    = "#B05000",
  "Telomerase-NotAmplified_ShortSurv"   = "#E0A870",
  "TNA_HighTC" = "#E0A999",
  
  # Non-TARGET — no survival data, plain TMM label
  "ALT"                                 = "#A8C8E8",   # lighter shade of ALT blue
  "Telomerase-Amplified"                = "#F0B0B0",   # lighter shade of Tel-Amp red
  "Telomerase-NotAmplified"             = "#F0D0A0",   # lighter shade of Tel-NotAmp orange
  
  # NO_TMM groups (no survival data available regardless of cohort)
  "NO_TMM"                              = "#418561",
  "NO_TMM_HighTC"                       = "#2A6B50"
)

present_colors <- group_colors[names(group_colors) %in% unique(train_metadata$Plot_Group)]

ht <- Heatmap(
  plot_matrix_scaled,
  name   = "Z-score",
  border = "black",
  
  cluster_rows     = TRUE,
  show_row_names   = FALSE,
  row_split        = row_split_vector,
  row_title_gp     = gpar(fontsize = 8, fontface = "bold"),
  
  column_split          = factor(train_metadata$Plot_Group),
  column_title_rot = 45, 
  column_title_gp       = gpar(fontsize = 5, fontface = "bold"),
  show_column_names     = TRUE,
  column_names_gp       = gpar(fontsize = 4),
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  
  col = col_fun,
  
  top_annotation = HeatmapAnnotation(
    TMM_Status           = train_metadata$Plot_Group,
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    col                  = list(TMM_Status = present_colors)
  )
)

pdf("Heatmap-tmm-training-equalsplit.pdf", width = 15, height = 15)
draw(ht)
dev.off()

#####

## Heatmap of TMM-ve upregulated genes in testing set.

test_metadata$Plot_Group <- as.character(test_metadata$TMM)


test_metadata$Vital_Status <- metadata_TARGET$Vital.Status[
  match(test_metadata$SampleID, metadata_TARGET$SampleID)]
test_metadata$OS_days <- metadata_TARGET$Overall.Survival.Time.in.Days[
  match(test_metadata$SampleID, metadata_TARGET$SampleID)]
test_metadata$OS_event <- ifelse(test_metadata$Vital_Status == "Dead", 1, 0)

tmm_groups_with_surv <- c("ALT", "Telomerase-Amplified", "Telomerase-NotAmplified")

surv_group_labels <- rep(NA_character_, nrow(test_metadata))

for (grp in tmm_groups_with_surv) {
  
  idx <- which(
    test_metadata$Cohort == "TARGET" &
      test_metadata$TMM == grp &
      !is.na(test_metadata$OS_days) &
      !is.na(test_metadata$OS_event)
  )
  
  if (length(idx) < 5) {
    surv_group_labels[idx] <- grp
    next
  }
  
  sub_meta <- test_metadata[idx, ]
  sub_meta$OS_days_var <- sub_meta$OS_days
  
  cut_res <- surv_cutpoint(
    sub_meta,
    time      = "OS_days",
    event     = "OS_event",
    variables = "OS_days_var"
  )
  cut_cat <- surv_categorize(cut_res)
  
  surv_group_labels[idx] <- ifelse(
    cut_cat$OS_days_var == "high",
    paste0(grp, "_LongSurv"),
    paste0(grp, "_ShortSurv")
  )
}



target_tmm_idx <- which(
  test_metadata$Cohort == "TARGET" &
    test_metadata$TMM %in% tmm_groups_with_surv &
    !is.na(surv_group_labels)
)
test_metadata$Plot_Group[target_tmm_idx] <- surv_group_labels[target_tmm_idx]

# NO_TMM high telomere content.
test_metadata$Plot_Group[
  test_metadata$TMM == "NO_TMM" & test_metadata$`Telomere Content` > 15
] <- "NO_TMM_HighTC"

test_metadata$Plot_Group[
  test_metadata$TMM == "Telomerase-Amplified" & test_metadata$`Telomere Content` > 15
] <- "TA_HighTC"

test_metadata$Plot_Group[
  test_metadata$TMM == "Telomerase-NotAmplified" & test_metadata$`Telomere Content` > 15
] <- "TNA_HighTC"

#####

selected_genes   <- c(candidate_genes_tmm_combined_up, candidate_genes_tmm_combined_down)
row_split_vector <- c(
  rep("Up-regulated",   length(candidate_genes_tmm_combined_up)),
  rep("Down-regulated", length(candidate_genes_tmm_combined_down))
)

plot_matrix        <- lcpm_test[selected_genes, ]
plot_matrix_scaled <- t(scale(t(plot_matrix)))
col_fun            <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))

group_colors <- c(
  # TARGET — with survival splits
  "ALT_LongSurv"                        = "#1A5C9E",
  "ALT_ShortSurv"                       = "#7EB6E8",
  "Telomerase-Amplified_LongSurv"       = "#9C1A1A",
  "Telomerase-Amplified_ShortSurv"      = "#E07070",
  "Telomerase-NotAmplified_LongSurv"    = "#B05000",
  "Telomerase-NotAmplified_ShortSurv"   = "#E0A870",
  "TNA_HighTC" = "#E0A999",
  
  # Non-TARGET — no survival data, plain TMM label
  "ALT"                                 = "#A8C8E8",   # lighter shade of ALT blue
  "Telomerase-Amplified"                = "#F0B0B0",   # lighter shade of Tel-Amp red
  "Telomerase-NotAmplified"             = "#F0D0A0",   # lighter shade of Tel-NotAmp orange
  
  # NO_TMM groups (no survival data available regardless of cohort)
  "NO_TMM"                              = "#418561",
  "NO_TMM_HighTC"                       = "#2A6B50"
)

present_colors <- group_colors[names(group_colors) %in% unique(test_metadata$Plot_Group)]

ht <- Heatmap(
  plot_matrix_scaled,
  name   = "Z-score",
  border = "black",
  
  cluster_rows     = TRUE,
  show_row_names   = FALSE,
  row_split        = row_split_vector,
  row_title_gp     = gpar(fontsize = 8, fontface = "bold"),
  
  column_split          = factor(test_metadata$Plot_Group),
  column_title_rot      = 90,
  column_title_gp       = gpar(fontsize = 6, fontface = "bold"),
  show_column_names     = TRUE,
  column_names_gp       = gpar(fontsize = 4),
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  
  col = col_fun,
  
  top_annotation = HeatmapAnnotation(
    TMM_Status           = test_metadata$Plot_Group,
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    col                  = list(TMM_Status = present_colors)
  )
)

pdf("Heatmap-tmm-testing-equalSplit.pdf", width = 10, height = 15)
draw(ht)
dev.off()

##################################################################################
## same thing for ALT.

## Heatmap of TMM-ve upregulated genes in training set.

selected_genes <- c(candidate_genes_alt_combined_up, candidate_genes_alt_combined_down)
row_split_vector <- c(rep("Up-regulated", length(candidate_genes_alt_combined_up)), 
                      rep("Down-regulated", length(candidate_genes_alt_combined_down)))


plot_matrix        <- lcpm[selected_genes, ]
plot_matrix_scaled <- t(scale(t(plot_matrix)))
col_fun            <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))

group_colors <- c(
  # TARGET — with survival splits
  "ALT_LongSurv"                        = "#1A5C9E",
  "ALT_ShortSurv"                       = "#7EB6E8",
  "Telomerase-Amplified_LongSurv"       = "#9C1A1A",
  "Telomerase-Amplified_ShortSurv"      = "#E07070",
  "Telomerase-NotAmplified_LongSurv"    = "#B05000",
  "Telomerase-NotAmplified_ShortSurv"   = "#E0A870",
  "TNA_HighTC" = "#E0A999",
  
  # Non-TARGET — no survival data, plain TMM label
  "ALT"                                 = "#A8C8E8",   # lighter shade of ALT blue
  "Telomerase-Amplified"                = "#F0B0B0",   # lighter shade of Tel-Amp red
  "Telomerase-NotAmplified"             = "#F0D0A0",   # lighter shade of Tel-NotAmp orange
  
  # NO_TMM groups (no survival data available regardless of cohort)
  "NO_TMM"                              = "#418561",
  "NO_TMM_HighTC"                       = "#2A6B50"
)


present_colors <- group_colors[names(group_colors) %in% unique(train_metadata$Plot_Group)]

ht <- Heatmap(
  plot_matrix_scaled,
  name   = "Z-score",
  border = "black",
  
  cluster_rows     = TRUE,
  show_row_names   = FALSE,
  row_split        = row_split_vector,
  row_title_gp     = gpar(fontsize = 8, fontface = "bold"),
  
  column_split          = factor(train_metadata$Plot_Group),
  column_title_rot = 45, 
  column_title_gp       = gpar(fontsize = 5, fontface = "bold"),
  show_column_names     = TRUE,
  column_names_gp       = gpar(fontsize = 4),
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  
  col = col_fun,
  
  top_annotation = HeatmapAnnotation(
    TMM_Status           = train_metadata$Plot_Group,
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    col                  = list(TMM_Status = present_colors)
  )
)

pdf("Heatmap-alt-training-equalSplit.pdf", width = 15, height = 15)
draw(ht)
dev.off()

#####

## Heatmap of TMM-ve upregulated genes in testing set.

selected_genes   <- c(candidate_genes_alt_combined_up, candidate_genes_alt_combined_down)
row_split_vector <- c(
  rep("Up-regulated",   length(candidate_genes_alt_combined_up)),
  rep("Down-regulated", length(candidate_genes_alt_combined_down))
)

plot_matrix        <- lcpm_test[selected_genes, ]
plot_matrix_scaled <- t(scale(t(plot_matrix)))
col_fun            <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))

group_colors <- c(
  # TARGET — with survival splits
  "ALT_LongSurv"                        = "#1A5C9E",
  "ALT_ShortSurv"                       = "#7EB6E8",
  "Telomerase-Amplified_LongSurv"       = "#9C1A1A",
  "Telomerase-Amplified_ShortSurv"      = "#E07070",
  "Telomerase-NotAmplified_LongSurv"    = "#B05000",
  "Telomerase-NotAmplified_ShortSurv"   = "#E0A870",
  "TNA_HighTC" = "#E0A999",
  
  # Non-TARGET — no survival data, plain TMM label
  "ALT"                                 = "#A8C8E8",   # lighter shade of ALT blue
  "Telomerase-Amplified"                = "#F0B0B0",   # lighter shade of Tel-Amp red
  "Telomerase-NotAmplified"             = "#F0D0A0",   # lighter shade of Tel-NotAmp orange
  
  # NO_TMM groups (no survival data available regardless of cohort)
  "NO_TMM"                              = "#418561",
  "NO_TMM_HighTC"                       = "#2A6B50"
)


present_colors <- group_colors[names(group_colors) %in% unique(test_metadata$Plot_Group)]

ht <- Heatmap(
  plot_matrix_scaled,
  name   = "Z-score",
  border = "black",
  
  cluster_rows     = TRUE,
  show_row_names   = FALSE,
  row_split        = row_split_vector,
  row_title_gp     = gpar(fontsize = 8, fontface = "bold"),
  
  column_split          = factor(test_metadata$Plot_Group),
  column_title_rot      = 90,
  column_title_gp       = gpar(fontsize = 6, fontface = "bold"),
  show_column_names     = TRUE,
  column_names_gp       = gpar(fontsize = 4),
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  
  col = col_fun,
  
  top_annotation = HeatmapAnnotation(
    TMM_Status           = test_metadata$Plot_Group,
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    col                  = list(TMM_Status = present_colors)
  )
)

pdf("Heatmap-alt-testing-equalSplit.pdf", width = 10, height = 15)
draw(ht)
dev.off()

##################################################################################

# # SIGNATURE BUILDING.
# 
# tmm_up_zscore <- run_harmonic_cv_selection_TMMve(expr_total = lcpm, meta_total = train_metadata,
#                                                       expr_test = lcpm_test, meta_test = test_metadata,
#                                                       candidate_genes = candidate_genes_tmm_combined_up,
#                                                       phenotype_col = "TMM", batch_col = "Cohort",
#                                                       label_neg = "NO_TMM",
#                                                       min_per_subgroup= 2,
#                                                       n_folds = 3, n_repeats = 5,
#                                                       n_cores = 8,
#                                                       max_genes = 15,
#                                                       seed_genes = NULL,
#                                                       lcb_conf = 0.95,
#                                                       lcb_boot_R = 500)
# 
# 
# tmm_down_zscore <- run_harmonic_cv_selection_TMMve(expr_total = lcpm, meta_total = train_metadata,
#                                                    expr_test = lcpm_test, meta_test = test_metadata,
#                                                    candidate_genes = candidate_genes_tmm_combined_down,
#                                                    phenotype_col = "TMM", batch_col = "Cohort",
#                                                    label_neg = "NO_TMM",
#                                                    min_per_subgroup= 2,
#                                                    n_folds = 3, n_repeats = 5,
#                                                    n_cores = 8,
#                                                    max_genes = 15,
#                                                    seed_genes = NULL,
#                                                    lcb_conf = 0.95,
#                                                    lcb_boot_R = 500)
# 
# alt_up_zscore <- run_harmonic_cv_selection_ALT(expr_total = lcpm, meta_total = train_metadata,
#                                                  expr_test = lcpm_test, meta_test = test_metadata,
#                                                  candidate_genes = candidate_genes_tmm_combined_up,
#                                                  phenotype_col = "TMM", batch_col = "Cohort",
#                                                  label_neg = "ALT",
#                                                  min_per_subgroup= 2,
#                                                  n_folds = 3, n_repeats = 5,
#                                                  n_cores = 8,
#                                                  max_genes = 15,
#                                                  seed_genes = NULL,
#                                                  lcb_conf = 0.95,
#                                                  lcb_boot_R = 500)
# 
# 
# 
# alt_down_zscore <- run_harmonic_cv_selection_ALT(expr_total = lcpm, meta_total = train_metadata,
#                                                expr_test = lcpm_test, meta_test = test_metadata,
#                                                candidate_genes = candidate_genes_tmm_combined_down,
#                                                phenotype_col = "TMM", batch_col = "Cohort",
#                                                label_neg = "ALT",
#                                                min_per_subgroup= 2,
#                                                n_folds = 3, n_repeats = 5,
#                                                n_cores = 8,
#                                                max_genes = 15,
#                                                seed_genes = NULL,
#                                                lcb_conf = 0.95,
#                                                lcb_boot_R = 500)
# 
# 
# #######################################################################################
# 
# tmm_singscore <- run_harmonic_cv_selection_TMMve_singscore(expr_total = lcpm , meta_total = train_metadata,
#                                                            expr_test = lcpm_test, meta_test = test_metadata,
#                                                            candidate_genes_up = candidate_genes_tmm_combined_up,
#                                                            candidate_genes_down = candidate_genes_tmm_combined_down,
#                                                            phenotype_col = "TMM", batch_col = "Cohort",
#                                                            label_neg = "NO_TMM",
#                                                            min_per_subgroup = 2,
#                                                            n_folds = 3, n_repeats = 5,
#                                                            n_cores = 8,
#                                                            max_genes = 20,
#                                                            seed_genes_up   = c("MYO9A", "FAXDC2", "ZSWIM6"),
#                                                            seed_genes_down = c("MRPL58", "SUV39H1", "SELENOH"),
#                                                            lcb_conf   = 0.95,
#                                                            lcb_boot_R = 500,
#                                                            perm_R     = 500)
# 
# 


###################################################################################
###################################################################################

# TMM_COLORS <- c(
#   "ALT"                     = "#2196F3", # Blue
#   "NO_TMM"                  = "#43A047", # Green
#   "TMM-ve"                  = "#43A047", # Green (Alias)
#   "Telomerase"              = "#E53935", # Red
#   "Telomerase-Amplified"    = "#FB8C00", # Orange
#   "Telomerase-NotAmplified" = "#8E24AA"  # Purple
# )
# 
# # This function handles data prep, plotting, and pagination for any dataset
# generate_gene_report <- function(expr_matrix, metadata, gene_list, output_file) {
# 
#   # A. Identify present genes
#   genes_present <- gene_list[gene_list %in% rownames(expr_matrix)]
#   if(length(genes_present) == 0) {
#     message("No genes from the list found in the provided matrix.")
#     return(NULL)
#   }
# 
#   # B. Prepare long-format data
#   lcpm_sub <- expr_matrix[genes_present, , drop = FALSE]
#   lcpm_long <- as.data.frame(t(lcpm_sub)) %>%
#     tibble::rownames_to_column("SampleID") %>%
#     pivot_longer(cols = -SampleID, names_to = "gene", values_to = "logCPM") %>%
#     left_join(metadata[, c("SampleID", "TMM")], by = "SampleID") %>%
#     mutate(
#       TMM = factor(TMM),
#       gene = factor(gene, levels = genes_present)
#     )
# 
#   # C. Internal Boxplot Generator
#   make_boxplot <- function(gene_name) {
#     df_gene <- lcpm_long %>% filter(gene == gene_name)
# 
#     # Define comparisons against ALT/NO_TMM
#     # Only create comparisons for levels that actually exist in this specific dataset
#     existing_levels <- levels(droplevels(df_gene$TMM))
#     other_levels <- setdiff(existing_levels, "NO_TMM")
#     comparisons <- lapply(other_levels, function(x) c("NO_TMM", x))
# 
#     ggplot(df_gene, aes(x = TMM, y = logCPM, fill = TMM)) +
#       geom_boxplot(outlier.shape = NA, width = 0.40, linewidth = 0.4) +
#       geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
#     
#       scale_fill_manual(values = TMM_COLORS) +
#       labs(title = gene_name, x = NULL, y = "log CPM") +
#       stat_compare_means(
#         comparisons = comparisons,
#         method = "t.test",
#         method.args = list(var.equal = FALSE),
#         label = "p.signif",
#         tip.length = 0.01,
#         step.increase = 0.05,
#         size = 3.5
#       ) +
#       theme_classic() +
#       theme(
#         plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#         legend.position = "none",
#         axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
#         axis.text.y = element_text(size = 10),
#         plot.margin = unit(c(5, 5, 5, 5), "pt")
#       )
#   }
# 
#   # D. Generate and Paginate Plots
#   all_plots <- lapply(genes_present, make_boxplot)
#   plots_per_page <- 9
#   n_pages <- ceiling(length(all_plots) / plots_per_page)
# 
#   pdf(output_file, width = 10, height = 14)
#   for (p in seq_len(n_pages)) {
#     idx_start <- (p - 1) * plots_per_page + 1
#     idx_end   <- min(p * plots_per_page, length(all_plots))
#     page_plots <- all_plots[idx_start:idx_end]
# 
#     # Pad with blank panels if needed
#     while (length(page_plots) < plots_per_page) {
#       page_plots[[length(page_plots) + 1]] <- ggplot() + theme_void()
#     }
# 
#     grid.arrange(
#       grobs = page_plots, ncol = 3, nrow = 3,
#       top = grid::textGrob(
#         label = paste0(output_file, " | Page ", p, " of ", n_pages),
#         gp = grid::gpar(fontsize = 10, fontface = "bold")
#       )
#     )
#   }
#   dev.off()
#   message("Report saved to: ", output_file)
# }
# 
# 
# target_genes <- c("LINC01783", "CCNB3", "LMNTD2", "DNAH14", "RNU5F-1")
# 
# # Train Report
# generate_gene_report(
#   expr_matrix = lcpm,
#   metadata = train_metadata,
#   gene_list = target_genes,
#   output_file = "pivot-ALTtraining.pdf"
# )
# 
# # Test Report
# generate_gene_report(
#   expr_matrix = lcpm_test,
#   metadata = test_metadata,
#   gene_list = target_genes,
#   output_file = "pivot_ALT_test.pdf"
# )

#######################################################################################

# target_genes <- c("MYO9A", "MRPL58", "FAXDC2", "SUV39H1", "SELENOH", "ZSWIM6", "NCAM2", "DUS1L",
#                   "TSEN54", "NUP85", "RNU6-722P", "VPS13C", "TEDC1", "SECISBP2L")
#
# # Generate Train Report
# generate_gene_report(
#   expr_matrix = lcpm,
#   metadata = train_metadata,
#   gene_list = target_genes,
#   output_file = "pivot-tmmtraining.pdf"
# )
#
# # Generate Test Report
# generate_gene_report(
#   expr_matrix = lcpm_test,
#   metadata = test_metadata,
#   gene_list = target_genes,
#   output_file = "pivot_tmm_test.pdf"
# )

####################################################################################
####################################################################################

## heatmap using pivot genes -- NO_TMM.

# Using the genes identified in your finalized clean_gene_list.txt
# new_gene_list <- c(
#   "MYO9A", "MRPL58", "FAXDC2", "SUV39H1", "SELENOH", "ZSWIM6", "NCAM2", "DUS1L",
#   "TSEN54", "NUP85", "RNU6-722P", "VPS13C", "TEDC1", "SECISBP2L"
# )
#
# col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))
#
# group_colors <- c(
#   "ALT_LongSurv"                      = "#1A5C9E",
#   "ALT_ShortSurv"                     = "#7EB6E8",
#   "Telomerase-Amplified_LongSurv"      = "#9C1A1A",
#   "Telomerase-Amplified_ShortSurv"     = "#E07070",
#   "Telomerase-NotAmplified_LongSurv"   = "#B05000",
#   "Telomerase-NotAmplified_ShortSurv"  = "#E0A870",
#   "ALT"                               = "#A8C8E8",
#   "Telomerase-Amplified"              = "#F0B0B0",
#   "Telomerase-NotAmplified"           = "#F0D0A0",
#   "NO_TMM"                            = "#418561",
#   "NO_TMM_HighTC"                     = "#2A6B50"
# )
#
# generate_split_heatmap <- function(expr_mat, metadata, genes, output_name, dataset_label) {
#
#   # A. Filter genes present in matrix
#   valid_genes <- intersect(genes, rownames(expr_mat))
#
#   # B. Calculate Directionality (Up vs Down in NO_TMM)
#   # We compare NO_TMM vs the average of everyone else
#   is_target <- metadata$TMM_Case == "NO_TMM"
#   target_means <- rowMeans(expr_mat[valid_genes, is_target, drop=FALSE])
#   other_means  <- rowMeans(expr_mat[valid_genes, !is_target, drop=FALSE])
#
#   up_genes   <- valid_genes[target_means > other_means]
#   down_genes <- valid_genes[target_means <= other_means]
#
#   ordered_genes <- c(up_genes, down_genes)
#   row_split_vec <- factor(c(rep("Up-regulated", length(up_genes)),
#                             rep("Down-regulated", length(down_genes))),
#                           levels = c("Up-regulated", "Down-regulated"))
#
#   # C. Prepare Matrix
#   plot_matrix <- t(scale(t(expr_mat[ordered_genes, metadata$SampleID])))
#
#   # D. Colors
#   present_colors <- group_colors[names(group_colors) %in% unique(metadata$Plot_Group)]
#
#   # E. Build Heatmap
#   ht <- Heatmap(
#     plot_matrix,
#     name   = "Z-score",
#     border = "black",
#     cluster_rows   = TRUE,
#     show_row_names = TRUE, # Enabled to see which new genes are where
#     row_names_gp   = gpar(fontsize = 6),
#     row_split      = row_split_vec,
#     row_title_gp   = gpar(fontsize = 10, fontface = "bold"),
#
#     column_split   = factor(metadata$Plot_Group),
#     column_title_rot = 45,
#     column_title_gp  = gpar(fontsize = 7, fontface = "bold"),
#     show_column_names = FALSE,
#     cluster_columns   = TRUE,
#     cluster_column_slices = FALSE,
#
#     col = col_fun,
#
#     top_annotation = HeatmapAnnotation(
#       TMM_Status = metadata$Plot_Group,
#       show_legend = FALSE,
#       show_annotation_name = FALSE,
#       col = list(TMM_Status = present_colors),
#       simple_anno_size = unit(0.5, "cm")
#     )
#   )
#
#   # F. Export
#   pdf(output_name, width = 16, height = 6)
#   draw(ht, column_title = paste("TMM Pivot Heatmap  -", dataset_label),
#        column_title_gp = gpar(fontsize = 16, fontface = "bold"))
#   dev.off()
#
#   message("Generated: ", output_name)
# }
# # Training Set
# generate_split_heatmap(
#   expr_mat = lcpm,
#   metadata = train_metadata,
#   genes = new_gene_list,
#   output_name = "Heatmap-pivot-Training.pdf",
#   dataset_label = "Training Set"
# )
#
# # Testing Set
# generate_split_heatmap(
#   expr_mat = lcpm_test,
#   metadata = test_metadata,
#   genes = new_gene_list,
#   output_name = "Heatmap-pivot-Testing.pdf",
#   dataset_label = "Testing Set"
# )
#
#
# #########################################################################################
#
# ## heatmap using pivot genes -- ALT
#
# new_gene_list <- c("ZNF285", "BBOX1-AS1", "LINC01783", "CCNB3", "TCAIM", "TADA3", "PTPRG",
#                    "OLFM2", "STRIP2", "MAGEA9B", "LINC01916", "XK", "MYZAP", "CYFIP1",
#                    "COL22A1", "OTULINL"
#
# )
#
# # --- 2. GLOBAL VISUAL CONFIGURATION ---
# col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#9C2424"))
#
# group_colors <- c(
#   "ALT_LongSurv"                      = "#1A5C9E",
#   "ALT_ShortSurv"                     = "#7EB6E8",
#   "Telomerase-Amplified_LongSurv"      = "#9C1A1A",
#   "Telomerase-Amplified_ShortSurv"     = "#E07070",
#   "Telomerase-NotAmplified_LongSurv"   = "#B05000",
#   "Telomerase-NotAmplified_ShortSurv"  = "#E0A870",
#   "ALT"                               = "#A8C8E8",
#   "Telomerase-Amplified"              = "#F0B0B0",
#   "Telomerase-NotAmplified"           = "#F0D0A0",
#   "NO_TMM"                            = "#418561",
#   "NO_TMM_HighTC"                     = "#2A6B50"
# )
#
# # --- 3. HELPER FUNCTION: PREPARE AND DRAW HEATMAP (ALT vs Non-ALT) ---
# generate_split_heatmap <- function(expr_mat, metadata, genes, output_name, dataset_label) {
#
#   # A. Filter genes present in matrix
#   valid_genes <- intersect(genes, rownames(expr_mat))
#
#   # B. Calculate Directionality (ALT vs Non-ALT)
#   # We compare ALT vs the average of everyone else (Telomerase + NO_TMM)
#   # Assuming 'TMM' column contains the strings 'ALT', 'Telomerase', and 'NO_TMM'
#   is_alt <- metadata$TMM == "ALT"
#
#   alt_means   <- rowMeans(expr_mat[valid_genes, is_alt, drop=FALSE])
#   other_means <- rowMeans(expr_mat[valid_genes, !is_alt, drop=FALSE])
#
#   up_genes   <- valid_genes[alt_means > other_means]
#   down_genes <- valid_genes[alt_means <= other_means]
#
#   ordered_genes <- c(up_genes, down_genes)
#   row_split_vec <- factor(c(rep("ALT Upregulated", length(up_genes)),
#                             rep("ALT Downregulated", length(down_genes))),
#                           levels = c("ALT Upregulated", "ALT Downregulated"))
#
#   # C. Prepare Matrix
#   plot_matrix <- t(scale(t(expr_mat[ordered_genes, metadata$SampleID])))
#
#   # D. Colors
#   present_colors <- group_colors[names(group_colors) %in% unique(metadata$Plot_Group)]
#
#   # E. Build Heatmap
#   ht <- Heatmap(
#     plot_matrix,
#     name   = "Z-score",
#     border = "black",
#     cluster_rows   = TRUE,
#     show_row_names = TRUE,
#     row_names_gp   = gpar(fontsize = 8),
#     row_split      = row_split_vec,
#     row_title_gp   = gpar(fontsize = 10, fontface = "bold"),
#
#     column_split   = factor(metadata$Plot_Group),
#     column_title_rot = 45,
#     column_title_gp  = gpar(fontsize = 7, fontface = "bold"),
#     show_column_names = FALSE,
#     cluster_columns   = TRUE,
#     cluster_column_slices = FALSE,
#
#     col = col_fun,
#
#     top_annotation = HeatmapAnnotation(
#       TMM_Status = metadata$Plot_Group,
#       show_legend = FALSE,
#       show_annotation_name = FALSE,
#       col = list(TMM_Status = present_colors),
#       simple_anno_size = unit(0.5, "cm")
#     )
#   )
#
#   # F. Export
#   # Width and height adjusted for a focused pivot gene set
#   pdf(output_name, width = 16, height = 6)
#   draw(ht, column_title = paste("ALT vs Non-ALT Pivot Heatmap -", dataset_label),
#        column_title_gp = gpar(fontsize = 16, fontface = "bold"))
#   dev.off()
#
#   message("Generated: ", output_name)
# }
#
# # --- 4. EXECUTION ---
#
# # Training Set
# generate_split_heatmap(
#   expr_mat = lcpm,
#   metadata = train_metadata,
#   genes = new_gene_list,
#   output_name = "Heatmap-ALT-Pivot-Training.pdf",
#   dataset_label = "Training Set"
# )
#
# # Testing Set
# generate_split_heatmap(
#   expr_mat = lcpm_test,
#   metadata = test_metadata,
#   genes = new_gene_list,
#   output_name = "Heatmap-ALT-Pivot-Testing.pdf",
#   dataset_label = "Testing Set"
# )

