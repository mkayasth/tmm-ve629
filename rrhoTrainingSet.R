library(RRHO2)

#####

# TARGET first.

# Loading bulk-seq log counts RNA and metadata file.
Expression <- read_delim("TARGET-NBL.star_counts.tsv")


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

# only taking training samples.
Expression <- Expression[, colnames(Expression) %in% train_metadata$SampleID]


###########################################################################################
##################################################################################################
gene_Expression <- round(2^Expression - 1)

train_metadata_target <- train_metadata %>%
  filter(Cohort == "TARGET")
train_metadata_target <- train_metadata_target %>%
  arrange(ALT_Case, TMM)
gene_Expression <- gene_Expression[, match(train_metadata_target$SampleID, colnames(gene_Expression))]

batch <- train_metadata_target$Cohort
group <- train_metadata_target$TMM


### Performing differential expression between Non-ALT & ALT & NO_TMM vs. TMM.

# building model matrix.

# creating differential gene expression object.
dge_TMM <- DGEList(counts=gene_Expression, group=train_metadata_target$Group)

keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# cpm and log cpm.
cpm_matrix <- cpm(dge_TMM, normalized.lib.sizes = TRUE)
tmm_lcpm_target_train <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE)


########

### Performing differential expression.
design <- model.matrix(~ 0 + TMM, data = train_metadata_target)
colnames(design) <- levels(train_metadata_target$TMM)
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
fitTMM_train_target <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVStmm"])
fitALT_train_target <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVSnonalt"])


# results.
top_tmm_train_target <- topTags(fitTMM_train_target, n = Inf)
top_ALT_train_target <- topTags(fitALT_train_target, n = Inf)


####################################################################################

# 0532 now.

geneExpression <- readRDS("GECounts_0532Pts_02202026.RDS")
geneExpression <- round(geneExpression)


geneExpression <- geneExpression[rownames(geneExpression) %in% valid_gene_anno$hgnc_symbol, , drop = FALSE]


# setting metadata order.
train_metadata_0532 <- train_metadata %>%
  filter(Cohort == "ANBL0532")

train_metadata_0532 <- train_metadata_0532 %>%
  arrange(ALT_Case, TMM)

colnames(geneExpression) <- sub("_.*", "", colnames(geneExpression))

geneExpression <- geneExpression[, colnames(geneExpression) %in% train_metadata_0532$SampleID]
geneExpression <- geneExpression[, match(train_metadata_0532$SampleID, colnames(geneExpression)), drop = FALSE]



##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.

train_metadata_0532$Group <- factor(train_metadata_0532$TMM, levels = c("ALT", "Telomerase", "NO_TMM"))


# creating differential gene expression object.
dge_TMM <- DGEList(counts=geneExpression,group=train_metadata_0532$Group)

keep <- filterByExpr(dge_TMM)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]


# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")

tmm_cpm  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_0532_train <- cpm(dge_TMM, normalized.lib.sizes = TRUE, log = TRUE)

########

### Performing differential expression.
design <- model.matrix(~ 0 + TMM, data = train_metadata_0532)
colnames(design) <- levels(train_metadata_0532$TMM)
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
fitTMM_train_0532 <- glmQLFTest(fit1, contrast = contrast_matrix[, "notmmVStmm"])
fitALT_train_0532 <- glmQLFTest(fit1, contrast = contrast_matrix[, "altVSnonalt"])


# results.
top_tmm_train_0532 <- topTags(fitTMM_train_0532, n = Inf)
top_ALT_train_0532 <- topTags(fitALT_train_0532, n = Inf)

##################################################################################
###################################################################################
##################################################################################

## Now, RRHO -- first TMM-ve vs. TMM+ve.
df_target <- as.data.frame(top_tmm_train_target)
df_0532 <- as.data.frame(top_tmm_train_0532)


# Prevent Inf values by replacing P-value 0 with the smallest double.
df_target$FDR[df_target$FDR == 0] <- .Machine$double.xmin
df_0532$FDR[df_0532$FDR == 0] <- .Machine$double.xmin

# Creating the ranking vectors: -log10(P) * sign(logFC)
# Higher score = Stronger TMM-ve signal
# Lower score  = Stronger TMM+ signal
rank_0532 <- df_0532$logFC * sqrt(pmax(df_0532$F, 0))
names(rank_0532) <- rownames(df_0532)

rank_target <- df_target$logFC * sqrt(pmax(df_target$F, 0))
names(rank_target) <- rownames(df_target)

# Aligning the gene sets (RRHO requires identical gene lists)
list_of_names <- list(names(rank_0532), names(rank_target))
common_genes <- Reduce(intersect, list_of_names)


list_0532 <- rank_0532[common_genes]
list_0532 <- as.data.frame(list_0532)
colnames(list_0532)[1] <- "Score"
list_0532$Gene <- rownames(list_0532)
list_0532 <- list_0532[, c("Gene", "Score")]


list_target <- rank_target[common_genes]
list_target <- as.data.frame(list_target)
colnames(list_target)[1] <- "Score"
list_target$Gene <- rownames(list_target)
list_target <- list_target[, c("Gene", "Score")]

# Running RRHO2
# log10.p = TRUE calculates the significance of the overlap
rrho_obj <- RRHO2_initialize(
  list_target,
  list_0532,
  labels = c("TMM_TARGET", "TMM_0532"),
  log10.ind = TRUE
)

RRHO2_heatmap(rrho_obj)


up_up_genes_target0532_tmm <- rrho_obj$genelist_uu$gene_list_overlap_uu
down_down_genes_target0532_tmm <- rrho_obj$genelist_dd$gene_list_overlap_dd

#################################################################################
#################################################################################

## Now, RRHO -- alt vs. alt-ve
df_target <- as.data.frame(top_ALT_train_target)

df_0532 <- as.data.frame(top_ALT_train_0532)


# Prevent Inf values by replacing P-value 0 with the smallest double.
#df_target$FDR[df_target$FDR == 0] <- .Machine$double.xmin
#df_0532$FDR[df_0532$FDR == 0] <- .Machine$double.xmin


# Creating the ranking vectors: -log10(P) * sign(logFC)
# Higher score = Stronger TMM-ve signal
# Lower score  = Stronger TMM+ signal
rank_0532 <- df_0532$logFC * sqrt(pmax(df_0532$F, 0))
names(rank_0532) <- rownames(df_0532)

rank_target <- df_target$logFC * sqrt(pmax(df_target$F, 0))
names(rank_target) <- rownames(df_target)

# Aligning the gene sets (RRHO requires identical gene lists)
list_of_names <- list(names(rank_0532), names(rank_target))
common_genes <- Reduce(intersect, list_of_names)


list_0532 <- rank_0532[common_genes]
list_0532 <- as.data.frame(list_0532)
colnames(list_0532)[1] <- "Score"
list_0532$Gene <- rownames(list_0532)
list_0532 <- list_0532[, c("Gene", "Score")]


list_target <- rank_target[common_genes]
list_target <- as.data.frame(list_target)
colnames(list_target)[1] <- "Score"
list_target$Gene <- rownames(list_target)
list_target <- list_target[, c("Gene", "Score")]

# Running RRHO2
# log10.p = TRUE calculates the significance of the overlap
rrho_obj2 <- RRHO2_initialize(
  list_target,
  list_0532,
  labels = c("ALT_TARGET", "ALT_0532"),
  log10.ind = TRUE
)

RRHO2_heatmap(rrho_obj2)


up_up_genes_target0532_alt <- rrho_obj2$genelist_uu$gene_list_overlap_uu
down_down_genes_target0532_alt <- rrho_obj2$genelist_dd$gene_list_overlap_dd

