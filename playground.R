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
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# Producing one plot per gene
plots <- lapply(c("MFSD6", "ZNF511", "XRCC3", "WDR47", "SUV39H1", "SLCO1A2", "BBOX1-AS1", "SLC5A12", "CCNB3", "BBOX1"), plot_gene)


cairo_pdf(
  "Pivot_Genes_Boxplots.pdf",
  width = 12,
  height = ceiling(length(plots) / 2) * 3
)

wrap_plots(plots, ncol = 2)

dev.off()

