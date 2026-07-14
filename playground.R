### Content: pivot gene visualization, auc calculator for signatures, line graph for random forest.

## seeing if pivot genes work across training and testing.

library(tidyverse)
library(patchwork)

# Order of groups
tmm_levels <- c("NO_TMM", "ALT", "Telomerase-Amplified", "Telomerase-NotAmplified")

plot_gene <- function(gene){
  
  ## Training
  train_df <- data.frame(
    Sample = colnames(lcpm),
    Expression = as.numeric(lcpm[gene, ]),
    TMM = factor(train_metadata$TMM, levels = tmm_levels),
    Dataset = "Training"
  )
  
  ## Testing
  test_df <- data.frame(
    Sample = colnames(lcpm_test),
    Expression = as.numeric(lcpm_test[gene, ]),
    TMM = factor(test_metadata$TMM, levels = tmm_levels),
    Dataset = "Testing"
  )
  
  plot_df <- bind_rows(train_df, test_df)
  
  ggplot(plot_df,
         aes(x = TMM,
             y = Expression,
             fill = TMM)) +
    geom_boxplot(outlier.shape = NA,
                 alpha = 0.7,
                 width = 0.65) +
    geom_jitter(width = 0.15,
                size = 1.8,
                alpha = 0.7) +
    facet_wrap(~Dataset, nrow = 1) +
    labs(
      title = gene,
      x = "",
      y = "logCPM"
    ) +
    scale_fill_manual(values = c(
      "NO_TMM" = "#4E79A7",
      "ALT" = "#F28E2B",
      "Telomerase" = "#59A14F"
    )) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
    ))
}

# Producing one plot per gene
plots <- lapply(c("SUV39H1", "FMO5", "PLCXD1", "KCTD21", "CAB39L", "ADRA1A", "OR56A3"), plot_gene)


cairo_pdf(
  "Pivot_Genes_Boxplots.pdf",
  width = 15,
  height = ceiling(length(plots)) * 3
)

wrap_plots(plots, ncol = 2)

dev.off()

##################################################################################

# AUC calculator.

library(singscore)
library(pROC)

calc_singscore_auc <- function(expr_mat,
                               metadata,
                               group_col,
                               positive_group,
                               negative_group,
                               up_genes,
                               down_genes = NULL,
                               centerScore = TRUE) {
  
  ## Keep only the two groups
  keep <- metadata[[group_col]] %in% c(positive_group, negative_group)
  
  expr_sub <- expr_mat[, keep, drop = FALSE]
  meta_sub <- metadata[keep, , drop = FALSE]
  
  ## Rank genes
  ranked <- rankGenes(as.matrix(expr_sub))
  
  ## Calculate singscore
  if (is.null(down_genes)) {
    scores <- simpleScore(
      rankData = ranked,
      upSet = up_genes,
      centerScore = centerScore
    )$TotalScore
  } else {
    scores <- simpleScore(
      rankData = ranked,
      upSet = up_genes,
      downSet = down_genes,
      centerScore = centerScore
    )$TotalScore
  }
  
  ## Binary outcome
  response <- factor(
    meta_sub[[group_col]],
    levels = c(negative_group, positive_group)
  )
  
  ## ROC/AUC
  roc_obj <- roc(
    response = response,
    predictor = scores,
    levels = c(negative_group, positive_group),
    direction = "<",
    quiet = TRUE
  )
  
  ## Cohen's d
  score_df <- data.frame(
    Score = scores,
    Group = response
  )
  
  d <- effsize::cohen.d(
    Score ~ Group,
    data = score_df,
    hedges.correction = TRUE
  )$estimate
  
  ## Distribution overlap
  overlap <- 2 * pnorm(-abs(d) / 2)
  
  list(
    auc = as.numeric(auc(roc_obj)),
    roc = roc_obj,
    scores = scores,
    cohens_d = unname(d),
    overlap = overlap
  )
}

############## Testing set.
res <- calc_singscore_auc(
  expr_mat = lcpm_test,
  metadata = test_metadata,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "CAB39L", "PLXNA4", "WSPAR", "HECW2", "PRSS51", "BLCAP", "ARHGAP21", "GRHL2", "RAPGEF5", "LINC00987"),
  down_genes = c("SUV39H1", "PLCXD1", "GALK1", "LINC01163")
)

res$auc
res$overlap

res_alt <- calc_singscore_auc(
  expr_mat = lcpm_test,
  metadata = test_metadata,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "DNAJB13"),
  down_genes = c("ADRA1A", "TERT"))

res_alt$auc
res_alt$overlap

### Upregulated signatures -- TMM-ve.

res_up <- calc_singscore_auc(
  expr_mat = lcpm_test,
  metadata = test_metadata,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "FAXDC2", "CAB39L", "MYORG", "FLRT3", "SLA")
)

res_up$auc
res_up$overlap


res_up_alt <- calc_singscore_auc(
  expr_mat = lcpm_test,
  metadata = test_metadata,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "SLCO1A2", "SCARA5", "PLK3", "DCAF8L1", "F8"))

res_up_alt$auc
res_up_alt$overlap
############## Training set.
res2 <- calc_singscore_auc(
  expr_mat = lcpm,
  metadata = train_metadata,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "CAB39L", "PLXNA4", "WSPAR", "HECW2", "PRSS51", "BLCAP", "ARHGAP21", "GRHL2", "RAPGEF5", "LINC00987"),
  down_genes = c("SUV39H1", "PLCXD1", "GALK1", "LINC01163")
)

res2$auc
res2$overlap

res_alt2 <- calc_singscore_auc(
  expr_mat = lcpm,
  metadata = train_metadata,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "DNAJB13"),
  down_genes = c("ADRA1A", "TERT"))

res_alt2$auc
res_alt2$overlap

### Upregulated signatures.

res_up2 <- calc_singscore_auc(
  expr_mat = lcpm,
  metadata = train_metadata,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "FAXDC2", "CAB39L", "MYORG", "FLRT3", "SLA")
)

res_up2$auc
res_up2$overlap

res_up_alt2 <- calc_singscore_auc(
  expr_mat = lcpm,
  metadata = train_metadata,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "SLCO1A2", "SCARA5", "PLK3", "DCAF8L1", "F8"))

res_up_alt2$auc
res_up_alt2$overlap

################### Ackerman dataset.
res3 <- calc_singscore_auc(
  expr_mat = expr_ackerman,
  metadata = metadata_ackerman,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "CAB39L", "PLXNA4", "HECW2",  "BLCAP", "ARHGAP21", "GRHL2", "RAPGEF5"),
  down_genes = c("SUV39H1", "PLCXD1", "GALK1")
)

res3$auc
res3$overlap

res_alt3 <- calc_singscore_auc(
  expr_mat = expr_ackerman,
  metadata = metadata_ackerman,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "DNAJB13"),
  down_genes = c("ADRA1A", "TERT"))

res_alt3$auc
res_alt3$overlap

### Upregulated signatures.

res_up3 <- calc_singscore_auc(
  expr_mat = expr_ackerman,
  metadata = metadata_ackerman,
  group_col = "TMM_Case",
  positive_group = "NO_TMM",
  negative_group = "TMM",
  up_genes = c("FMO5", "KCTD21", "FAXDC2", "CAB39L", "MYORG", "FLRT3", "SLA")
)

res_up3$auc
res_up3$overlap

res_up_alt3 <- calc_singscore_auc(
  expr_mat = expr_ackerman,
  metadata = metadata_ackerman,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "SLCO1A2", "SCARA5", "PLK3", "DCAF8L1", "F8"))

res_up_alt3$auc
res_up_alt3$overlap
###################################################################################
###################################################################################
###################################################################################

#z Line graph -- genes TMM-ve.
library(readxl)
library(tidyverse)

df <- read_excel("tmmSignature.xlsx")

# removing pivot genes.
df_plot <- df %>%
  filter(method == "SINGSCORE")

# preserving order.
df_plot$Genes <- factor(df_plot$Genes)


# auc plot.
pivot_n <- 5
auc_long <- df_plot %>%
  dplyr::select(Genes, `auc_training`, `auc_testing`) %>%
  pivot_longer(cols = -Genes, names_to = "Dataset", values_to = "AUC") %>%
  mutate(AUC = as.numeric(AUC))

auc_long <- auc_long %>%
  mutate(
    AUC = as.numeric(AUC)
  ) %>%
  filter(!is.na(AUC))

auc_long <- auc_long %>%
  group_by(Dataset) %>%
  mutate(
    Added_Genes = row_number(),
    Signature_Size = pivot_n + Added_Genes
  ) %>%
  ungroup()

pivot_labels <- c(
  "SUV39H1", "FMO5", "KCTD21", "CAB39L", "PLCXD1"
)

selection_labels <- auc_long %>%
  filter(Dataset == unique(Dataset)[1]) %>%
  arrange(Signature_Size) %>%
  pull(Genes) %>%
  as.character()

all_labels <- c(pivot_labels, selection_labels)

all_breaks <- c(
  1:pivot_n,
  sort(unique(auc_long$Signature_Size))
)

all_breaks <- seq_along(all_labels)

highlight_genes <- c("GRHL2")

vline_df <- auc_long %>%
  filter(Genes %in% highlight_genes) %>%
  distinct(Genes, Signature_Size)


ggplot(
  auc_long,
  aes(
    x = Signature_Size,
    y = AUC,
    color = Dataset,
    group = Dataset
  )
) +
  
  # Pivot region
  geom_rect(
    xmin = 0,
    xmax = pivot_n,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey90",
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Training/Test curves
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  
  # Labels
  annotate(
    "text",
    x = pivot_n/2,
    y = 0.92,
    label = "Pivot Genes",
    fontface = "bold",
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = all_breaks,
    labels = all_labels,
    limits = c(1, max(all_breaks))
  )+
  
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.6, 1.00, 0.1)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = 0,   # removes line
        shape = 16,
        linewidth = 0
      )
    )
  ) +
  
  labs(
    x = "Total Signature Size",
    y = "Harmonic Mean AUC",
    color = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "right",
    axis.title = element_text(
      face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
    
  ) +
  
  geom_vline(
    data = vline_df,
    aes(xintercept = Signature_Size),
    linetype = "dashed"
  )

###################################################################################
###################################################################################
###################################################################################

#z Line graph -- genes ALT.
library(readxl)
library(tidyverse)

df <- read_excel("altSignature.xlsx")

# removing pivot genes.
df_plot <- df %>%
  filter(method == "SINGSCORE")

# preserving order.
df_plot$Genes <- factor(df_plot$Genes)


# auc plot.
pivot_n <- 5
auc_long <- df_plot %>%
  dplyr::select(Genes, `auc_training`, `auc_testing`) %>%
  pivot_longer(cols = -Genes, names_to = "Dataset", values_to = "AUC") %>%
  mutate(AUC = as.numeric(AUC))

auc_long <- auc_long %>%
  mutate(
    AUC = as.numeric(AUC)
  ) %>%
  filter(!is.na(AUC))

auc_long <- auc_long %>%
  group_by(Dataset) %>%
  mutate(
    Added_Genes = row_number(),
    Signature_Size = pivot_n + Added_Genes
  ) %>%
  ungroup()

pivot_labels <- c(
  "LMNTD2", "LINC01783", "CCNB3", "OR56A3", "ADRA1A"
)

selection_labels <- auc_long %>%
  filter(Dataset == unique(Dataset)[1]) %>%
  arrange(Signature_Size) %>%
  pull(Genes) %>%
  as.character()

all_labels <- c(pivot_labels, selection_labels)

all_breaks <- c(
  1:pivot_n,
  sort(unique(auc_long$Signature_Size))
)

all_breaks <- seq_along(all_labels)

highlight_genes <- c("DNAJB13")

vline_df <- auc_long %>%
  filter(Genes %in% highlight_genes) %>%
  distinct(Genes, Signature_Size)


ggplot(
  auc_long,
  aes(
    x = Signature_Size,
    y = AUC,
    color = Dataset,
    group = Dataset
  )
) +
  
  # Pivot region
  geom_rect(
    xmin = 0,
    xmax = pivot_n,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey90",
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Training/Test curves
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  
  # Labels
  annotate(
    "text",
    x = pivot_n/2,
    y = 0.92,
    label = "Pivot Genes",
    fontface = "bold",
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = all_breaks,
    labels = all_labels,
    limits = c(1, max(all_breaks))
  )+
  
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.6, 1.00, 0.1)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = 0,   # removes line
        shape = 16,
        linewidth = 0
      )
    )
  ) +
  
  labs(
    x = "Total Signature Size",
    y = "Harmonic Mean AUC",
    color = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "right",
    axis.title = element_text(
      face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
    
  ) +
  
  geom_vline(
    data = vline_df,
    aes(xintercept = Signature_Size),
    linetype = "dashed"
  )

###################################################################################
###################################################################################
###################################################################################

#z Line graph -- upregulated genes TMM-ve.
library(readxl)
library(tidyverse)

df <- read_excel("tmmSignatureUpregulated.xlsx")

# removing pivot genes.
df_plot <- df %>%
  filter(method == "SINGSCORE")

# preserving order.
df_plot$Genes <- factor(df_plot$Genes)


# auc plot.
pivot_n <- 5
auc_long <- df_plot %>%
  dplyr::select(Genes, `auc_training`, `auc_testing`) %>%
  pivot_longer(cols = -Genes, names_to = "Dataset", values_to = "AUC") %>%
  mutate(AUC = as.numeric(AUC))

auc_long <- auc_long %>%
  mutate(
    AUC = as.numeric(AUC)
  ) %>%
  filter(!is.na(AUC))

auc_long <- auc_long %>%
  group_by(Dataset) %>%
  mutate(
    Added_Genes = row_number(),
    Signature_Size = pivot_n + Added_Genes
  ) %>%
  ungroup()

pivot_labels <- c(
  "FMO5", "KCTD21", "FAXDC2", "CAB39L", "MYORG"
)

selection_labels <- auc_long %>%
  filter(Dataset == unique(Dataset)[1]) %>%
  arrange(Signature_Size) %>%
  pull(Genes) %>%
  as.character()

all_labels <- c(pivot_labels, selection_labels)

all_breaks <- c(
  1:pivot_n,
  sort(unique(auc_long$Signature_Size))
)

all_breaks <- seq_along(all_labels)

highlight_genes <- c("GRHL2")

vline_df <- auc_long %>%
  filter(Genes %in% highlight_genes) %>%
  distinct(Genes, Signature_Size)


ggplot(
  auc_long,
  aes(
    x = Signature_Size,
    y = AUC,
    color = Dataset,
    group = Dataset
  )
) +
  
  # Pivot region
  geom_rect(
    xmin = 0,
    xmax = pivot_n,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey90",
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Training/Test curves
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  
  # Labels
  annotate(
    "text",
    x = pivot_n/2,
    y = 0.92,
    label = "Pivot Genes",
    fontface = "bold",
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = all_breaks,
    labels = all_labels,
    limits = c(1, max(all_breaks))
  )+
  
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.6, 1.00, 0.1)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = 0,   # removes line
        shape = 16,
        linewidth = 0
      )
    )
  ) +
  
  labs(
    x = "Total Signature Size",
    y = "Harmonic Mean AUC",
    color = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "right",
    axis.title = element_text(
      face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
    
  ) +
  
  geom_vline(
    data = vline_df,
    aes(xintercept = Signature_Size),
    linetype = "dashed"
  )

###################################################################################
###################################################################################
###################################################################################

#z Line graph -- upregulated genes ALT.
library(readxl)
library(tidyverse)

df <- read_excel("altSignatureUpregulated.xlsx")

# removing pivot genes.
df_plot <- df %>%
  filter(method == "SINGSCORE")

# preserving order.
df_plot$Genes <- factor(df_plot$Genes)


# auc plot.
pivot_n <- 5
auc_long <- df_plot %>%
  dplyr::select(Genes, `auc_training`, `auc_testing`) %>%
  pivot_longer(cols = -Genes, names_to = "Dataset", values_to = "AUC") %>%
  mutate(AUC = as.numeric(AUC))

auc_long <- auc_long %>%
  mutate(
    AUC = as.numeric(AUC)
  ) %>%
  filter(!is.na(AUC))

auc_long <- auc_long %>%
  group_by(Dataset) %>%
  mutate(
    Added_Genes = row_number(),
    Signature_Size = pivot_n + Added_Genes
  ) %>%
  ungroup()

pivot_labels <- c(
  "LINC01783", "LMNTD2", "CCNB3", "OR56A3", "SLCO1A2"
)

selection_labels <- auc_long %>%
  filter(Dataset == unique(Dataset)[1]) %>%
  arrange(Signature_Size) %>%
  pull(Genes) %>%
  as.character()

all_labels <- c(pivot_labels, selection_labels)

all_breaks <- c(
  1:pivot_n,
  sort(unique(auc_long$Signature_Size))
)

all_breaks <- seq_along(all_labels)

highlight_genes <- c("F8")

vline_df <- auc_long %>%
  filter(Genes %in% highlight_genes) %>%
  distinct(Genes, Signature_Size)


ggplot(
  auc_long,
  aes(
    x = Signature_Size,
    y = AUC,
    color = Dataset,
    group = Dataset
  )
) +
  
  # Pivot region
  geom_rect(
    xmin = 0,
    xmax = pivot_n,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey90",
    color = "grey50",
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Training/Test curves
  geom_line(linewidth = 1.3) +
  geom_point(size = 2) +
  
  # Labels
  annotate(
    "text",
    x = pivot_n/2,
    y = 0.92,
    label = "Pivot Genes",
    fontface = "bold",
    size = 4
  ) +
  
  scale_x_continuous(
    breaks = all_breaks,
    labels = all_labels,
    limits = c(1, max(all_breaks))
  )+
  
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.6, 1.00, 0.1)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = 0,   # removes line
        shape = 16,
        linewidth = 0
      )
    )
  ) +
  
  labs(
    x = "Total Signature Size",
    y = "Harmonic Mean AUC",
    color = NULL
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "right",
    axis.title = element_text(
      face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
    
  ) +
  
  geom_vline(
    data = vline_df,
    aes(xintercept = Signature_Size),
    linetype = "dashed"
  )

###################################################################################
###################################################################################
###################################################################################