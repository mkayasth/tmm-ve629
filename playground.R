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
plots <- lapply(c("SUV39H1", "FMO5", "KCTD21", "TRPM3", "CAB39L", "RFC2", "LMNTD2", "LINC01783", "CCNB3", "OR56A3", "BBOX1-AS1"), plot_gene)


cairo_pdf(
  "Pivot_Genes_Boxplots.pdf",
  width = 15,
  height = ceiling(length(plots)) * 3
)

wrap_plots(plots, ncol = 2)

dev.off()

