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
  
  list(
    auc = as.numeric(auc(roc_obj)),
    roc = roc_obj,
    scores = scores
  )
}

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

res_alt2 <- calc_singscore_auc(
  expr_mat = lcpm_test,
  metadata = test_metadata,
  group_col = "ALT_Case",
  positive_group = "ALT",
  negative_group = "Non-ALT",
  up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "DNAJB13"),
  down_genes = c("ADRA1A", "TERT"))

res_alt2$auc
