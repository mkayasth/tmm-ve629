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
lcpm_all <- cpm(dge_TMM, log = TRUE)

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

set.seed(17) # For reproducibility83, 41, 17.
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
  notmmVStmm = NO_TMM - (ALT + Telomerase.Amplified + Telomerase.NotAmplified)/3,
  altVSnonalt = ALT - (Telomerase.Amplified + Telomerase.NotAmplified + NO_TMM)/3,
  levels = design
)

# differential expression test.
fitTMM_combined <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVStmm"])
fitALT_combined <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVSnonalt"])


# results.
top_tmm_combined <- topTags(fitTMM_combined, n = Inf)
top_ALT_combined <- topTags(fitALT_combined, n = Inf)


# filtering for candidate genes.
candidate_genes_TMM_combined <- subset(top_tmm_combined$table, FDR <= 0.05 & abs(logFC) >= 0.5)
candidate_genes_tmm_combined <- rownames(candidate_genes_TMM_combined)

candidate_genes_ALT_combined <- subset(top_ALT_combined$table, FDR <= 0.05 & abs(logFC) >= 0.5)
candidate_genes_alt_combined <- rownames(candidate_genes_ALT_combined)

######################################################################################

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
## including RRHO2 filtered genes.
source("tmm-ve629/rrhoTrainingSet.R")

candidate_genes_tmm_combined_up <- rownames(candidate_genes_TMM_combined[candidate_genes_TMM_combined$logFC > 0 &
                                                                           rownames(candidate_genes_TMM_combined) %in% rownames(lcpm), ])
candidate_genes_tmm_combined_down <- rownames(candidate_genes_TMM_combined[candidate_genes_TMM_combined$logFC < 0 &
                                                                             rownames(candidate_genes_TMM_combined) %in% rownames(lcpm), ])

candidate_genes_alt_combined_up <- rownames(candidate_genes_ALT_combined[candidate_genes_ALT_combined$logFC > 0 &
                                                                           rownames(candidate_genes_ALT_combined) %in% rownames(lcpm), ])
candidate_genes_alt_combined_down <- rownames(candidate_genes_ALT_combined[candidate_genes_ALT_combined$logFC < 0 &
                                                                             rownames(candidate_genes_ALT_combined) %in% rownames(lcpm), ])


candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[candidate_genes_tmm_combined_up %in% up_up_genes_target0532_tmm]
candidate_genes_tmm_combined_down <- candidate_genes_tmm_combined_down[candidate_genes_tmm_combined_down %in% down_down_genes_target0532_tmm]

candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[candidate_genes_alt_combined_up %in% up_up_genes_target0532_alt]
candidate_genes_alt_combined_down <- candidate_genes_alt_combined_down[candidate_genes_alt_combined_down %in% down_down_genes_target0532_alt]

###################################################################################

# Step 2: Random Forest: first for TMM-ve.


# Running RF multiple times for stability.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 10,
  "ALT"    = 10,
  "Telomerase-Amplified" = 5,
  "Telomerase-NotAmplified" = 5
)

candidate_genes_tmm_combined <- candidate_genes_tmm_combined[candidate_genes_tmm_combined %in%
                                                               candidate_genes_tmm_combined_up |
                                                               candidate_genes_tmm_combined %in%
                                                               candidate_genes_tmm_combined_down]

candidate_genes_tmm_combined <- candidate_genes_tmm_combined[candidate_genes_tmm_combined %in% rownames(lcpm)]

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
  filter(mean_no_tmm > 0, mean_tmm > 0, mean_accuracy > 0, mean_gini > 0) %>%
  arrange(avg_accuracy_rank_in_seeds)
  #mutate(gene_index = row_number())


# Visualization: Stability vs. Magnitude
p2 <- ggplot(stability_metrics, aes(x = seeds_present_in_top, y = mean_accuracy)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#2E86AB") +
  geom_text_repel(
    data = subset(stability_metrics, seeds_present_in_top >= 6),
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
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

# Top 20 Lollipop Plot
top_20_tmm <- stability_metrics %>% slice_head(n = 20)

p3 <- ggplot(top_20_tmm, aes(x = reorder(gene, avg_accuracy_rank_in_seeds), y = avg_accuracy_rank_in_seeds)) +
  geom_segment(aes(xend = gene, y = 0, yend = avg_accuracy_rank_in_seeds), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Average Gene Ranks") +
  theme_classic()

print(p2)
print(p3)

# Step 2: Random Forest: first for ALT.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 10,
  "ALT"    = 10,
  "Telomerase-Amplified" = 5,
  "Telomerase-NotAmplified" = 5
)

candidate_genes_alt_combined <- candidate_genes_alt_combined[candidate_genes_alt_combined %in%
                                                               candidate_genes_alt_combined_up |
                                                               candidate_genes_alt_combined %in%
                                                               candidate_genes_alt_combined_down]
candidate_genes_alt_combined <- candidate_genes_alt_combined[candidate_genes_alt_combined %in% rownames(lcpm)]

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
  filter(mean_alt > 0, mean_nonalt > 0, mean_accuracy > 0, mean_gini > 0) %>%
  arrange(avg_accuracy_rank_in_seeds)


# Visualization: Stability vs. Magnitude
p5 <- ggplot(stability_metrics2, aes(x = seeds_present_in_top, y = mean_accuracy)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#2E86AB") +
  geom_text_repel(
    data = subset(stability_metrics2, seeds_present_in_top >= 7),
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
    plot.title   = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

# Top 20 Lollipop Plot
top_20_alt <- stability_metrics2 %>% slice_head(n = 20)

p6 <- ggplot(top_20_alt, aes(x = reorder(gene, avg_accuracy_rank_in_seeds), y = avg_accuracy_rank_in_seeds)) +
  geom_segment(aes(xend = gene, y = 0, yend = avg_accuracy_rank_in_seeds), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Mean Decrease Accuracy") +
  theme_classic()

print(p5)
print(p6)


#####
candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[candidate_genes_alt_combined_up %in% stability_metrics2$gene]
stability_metrics2[stability_metrics2$gene %in% candidate_genes_alt_combined_up, ]

candidate_genes_alt_combined_down <- candidate_genes_alt_combined_down[candidate_genes_alt_combined_down %in% stability_metrics2$gene]
stability_metrics2[stability_metrics2$gene %in% candidate_genes_alt_combined_down, ]

candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[candidate_genes_tmm_combined_up %in% stability_metrics$gene]
stability_metrics[stability_metrics$gene %in% candidate_genes_tmm_combined_up, ]

candidate_genes_tmm_combined_down <- candidate_genes_tmm_combined_down[candidate_genes_tmm_combined_down %in% stability_metrics$gene]
stability_metrics[stability_metrics$gene %in% candidate_genes_tmm_combined_down, ]

#####################

stability_metrics$combined_rank <- (
  2 * stability_metrics$avg_accuracy_rank_in_seeds +
    stability_metrics$avg_rank_in_Gini_in_seeds
) / 3
stability_metrics <- stability_metrics %>%
  arrange((combined_rank))

stability_metrics2$combined_rank <- (
  2 * stability_metrics2$avg_accuracy_rank_in_seeds +
    stability_metrics2$avg_Gini_rank_in_seeds
) / 3
stability_metrics2 <- stability_metrics2 %>%
  arrange((combined_rank))



#########

tmm_signature <- run_harmonic_cv_selection_singscore(expr_total = lcpm, meta_total = train_metadata,
                                                    expr_test_list = list("Testing Set" = lcpm_test),  # named list of expression matrices
                                                    meta_test_list = list("Testing Set" = test_metadata),   # named list of metadata data frames
                                                    candidate_genes_up = candidate_genes_tmm_combined_up,
                                                    candidate_genes_down = candidate_genes_tmm_combined_down,
                                                    phenotype_col = "TMM", batch_col = "Cohort",
                                                    label_neg = "NO_TMM",
                                                    min_per_subgroup = 2,
                                                    n_folds = 3 , n_repeats = 5,
                                                    n_cores = 8,
                                                    max_genes = 15,
                                                    pivot_genes = c(SUV39H1 = "DOWN", FMO5 = "UP",
                                                                    KCTD21 = "UP", CAB39L = "UP"),
                                                    n_pivots    = 4L,
                                                    lcb_conf   = 0.95,
                                                    lcb_boot_R = 500,
                                                    perm_R  = 0)
save.image()

alt_signature <- run_harmonic_cv_selection_singscore(expr_total = lcpm, meta_total = train_metadata,
                                                     expr_test_list = list("Testing Set" = lcpm_test),  # named list of expression matrices
                                                     meta_test_list = list("Testing Set" = test_metadata),   # named list of metadata data frames
                                                     candidate_genes_up = candidate_genes_alt_combined_up,
                                                     candidate_genes_down = candidate_genes_alt_combined_down,
                                                     phenotype_col = "TMM", batch_col = "Cohort",
                                                     label_neg = "ALT",
                                                     min_per_subgroup = 2,
                                                     n_folds = 3 , n_repeats = 5,
                                                     n_cores = 8,
                                                     max_genes = 15,
                                                     pivot_genes = c(LMNTD2 = "UP", LINC01783 = "UP",
                                                                     CCNB3 = "UP", OR56A3 = "UP"),
                                                     n_pivots    = 4L,
                                                     lcb_conf   = 0.95,
                                                     lcb_boot_R = 500,
                                                     perm_R  = 00)

save.image()

#####################################################################################
#####################################################################################
#####################################################################################

# Taking only upregulated genes -- starting with random forest.

# Running RF multiple times for stability.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 10,
  "ALT"    = 10,
  "Telomerase-Amplified" = 5,
  "Telomerase-NotAmplified" = 5
)

candidate_genes_tmm_combined_up <- candidate_genes_tmm_combined_up[candidate_genes_tmm_combined_up %in% rownames(lcpm)]

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  rf <- randomForest(
    x = t(lcpm[candidate_genes_tmm_combined_up, ]),
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
stability_metrics3 <- all_seed_data %>%
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
  filter(mean_no_tmm > 0, mean_tmm > 0, mean_accuracy > 0, mean_gini > 0) %>%
  arrange(avg_accuracy_rank_in_seeds)
#mutate(gene_index = row_number())


# Visualization: Stability vs. Magnitude
p7 <- ggplot(stability_metrics3, aes(x = seeds_present_in_top, y = mean_accuracy)) +
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
    y = "Average Ranks"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14)
  )

# Top 20 Lollipop Plot
top_20_tmm <- stability_metrics3 %>% slice_head(n = 20)

p8 <- ggplot(top_20_tmm, aes(x = reorder(gene, avg_accuracy_rank_in_seeds), y = avg_accuracy_rank_in_seeds)) +
  geom_segment(aes(xend = gene, y = 0, yend = avg_accuracy_rank_in_seeds), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Mean Decrease Accuracy") +
  theme_classic()

print(p7)
print(p8)

# Step 2: Random Forest: first for ALT.

seeds <- c(49, 77, 123, 134, 531, 424, 636, 4562, 46464, 55443)
top_n_threshold <- 20
importance_list <- list()


sample_sizes <- c(
  "NO_TMM" = 10,
  "ALT"    = 10,
  "Telomerase-Amplified" = 5,
  "Telomerase-NotAmplified" = 5
)

candidate_genes_alt_combined_up <- candidate_genes_alt_combined_up[candidate_genes_alt_combined_up %in% rownames(lcpm)]

for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  rf <- randomForest(
    x = t(lcpm[candidate_genes_alt_combined_up, ]),
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
stability_metrics4 <- all_seed_data2 %>%
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
  filter(mean_alt > 0, mean_nonalt > 0, mean_accuracy > 0, mean_gini > 0) %>%
  arrange(avg_accuracy_rank_in_seeds)

# Visualization: Stability vs. Magnitude
p9 <- ggplot(stability_metrics4, aes(x = seeds_present_in_top, y = mean_accuracy)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#2E86AB") +
  geom_text_repel(
    data = subset(stability_metrics2, seeds_present_in_top >= 9),
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
top_20_alt <- stability_metrics4 %>% slice_head(n = 20)

p10 <- ggplot(top_20_alt, aes(x = reorder(gene, avg_accuracy_rank_in_seeds), y = avg_accuracy_rank_in_seeds)) +
  geom_segment(aes(xend = gene, y = 0, yend = avg_accuracy_rank_in_seeds), linewidth = 0.6) +
  geom_point(size = 4, color = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 20 Pivot Genes (Averaged)", x = NULL, y = "Averaged Gene Ranks") +
  theme_classic()

print(p9)
print(p10)

tmm_signature2 <- run_harmonic_cv_selection_singscore(expr_total = lcpm, meta_total = train_metadata,
                                                     expr_test_list = list("Testing Set" = lcpm_test),  # named list of expression matrices
                                                     meta_test_list = list("Testing Set" = test_metadata),   # named list of metadata data frames
                                                     candidate_genes_up = candidate_genes_tmm_combined_up,
                                                     candidate_genes_down = NULL,
                                                     phenotype_col = "TMM", batch_col = "Cohort",
                                                     label_neg = "NO_TMM",
                                                     min_per_subgroup = 2,
                                                     n_folds = 3 , n_repeats = 5,
                                                     n_cores = 8,
                                                     max_genes = 15,
                                                     pivot_genes = c(SUV39H1 = "DOWN", FMO5 = "UP",
                                                                     KCTD21 = "UP", CAB39L = "UP"),
                                                     n_pivots    = 4L,
                                                     lcb_conf   = 0.95,
                                                     lcb_boot_R = 500,
                                                     perm_R  = 0)

save.image()

alt_signature2 <- run_harmonic_cv_selection_singscore(expr_total = lcpm, meta_total = train_metadata,
                                                     expr_test_list = list("Testing Set" = lcpm_test),  # named list of expression matrices
                                                     meta_test_list = list("Testing Set" = test_metadata),   # named list of metadata data frames
                                                     candidate_genes_up = candidate_genes_alt_combined_up,
                                                     candidate_genes_down = NULL,
                                                     phenotype_col = "TMM", batch_col = "Cohort",
                                                     label_neg = "ALT",
                                                     min_per_subgroup = 2,
                                                     n_folds = 3 , n_repeats = 5,
                                                     n_cores = 8,
                                                     max_genes = 15,
                                                     pivot_genes = c(LINC01783 = "UP", LMNTD2 = "UP",
                                                                     CCNB3 = "UP"),
                                                     n_pivots    = 3L,
                                                     lcb_conf   = 0.95,
                                                     lcb_boot_R = 500,
                                                     perm_R  = 200)

save.image()
