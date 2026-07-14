### Contents: boxplots for signature.

# TMM-ve first -- training set.

up_genes = c("FMO5", "KCTD21", "CAB39L", "PLXNA4", "WSPAR", "HECW2", "PRSS51", "BLCAP", "ARHGAP21", "GRHL2", "RAPGEF5", "LINC00987")
down_genes = c("SUV39H1", "PLCXD1", "GALK1", "LINC01163")

ranked_lcpm <- rankGenes(as.matrix(lcpm))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  train_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "ALT",
  "NO_TMM",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c("#0a75ad", "#666666", "#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
  linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore Scores (Training Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## TMM-ve first -- testing set.
ranked_lcpm <- rankGenes(as.matrix(lcpm_test))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  test_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "ALT",
  "NO_TMM",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c("#0a75ad", "#666666", "#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "errorbar",
    width = 0.8,
    linewidth = 1
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Testing Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## TMM-ve first -- ackerman set.
ranked_lcpm <- rankGenes(as.matrix(expr_ackerman))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  metadata_ackerman[, c("SampleID", "TMM_Category")],
  by = "SampleID"
)

# Plot
groups <- c(
  "ALT",
  "NO_TMM",
  "Telomerase"

)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c("#0a75ad", "#666666", "#f08080")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Ackerman)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

####################################################################################
####################################################################################
####################################################################################

## Now, ALT -- training set.

up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "DNAJB13")
down_genes = c("ADRA1A", "TERT")

ranked_lcpm <- rankGenes(as.matrix(lcpm))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  train_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c("#666666", "#0a75ad","#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore Scores (Training Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## aLT -- testing set.
ranked_lcpm <- rankGenes(as.matrix(lcpm_test))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  test_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c( "#666666", "#0a75ad", "#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "errorbar",
    width = 0.8,
    linewidth = 1
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Testing Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## ALT -- ackerman set.
ranked_lcpm <- rankGenes(as.matrix(expr_ackerman))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  metadata_ackerman[, c("SampleID", "TMM_Category")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase"
  
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c( "#666666", "#0a75ad", "#f08080")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Ackerman)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

##################################################################################
##################################################################################
##################################################################################

## Now, ALT upregulated -- training set.

up_genes = c("LMNTD2", "LINC01783", "CCNB3", "OR56A3", "SLCO1A2", "SCARA5", "PLK3", "DCAF8L1", "F8")

ranked_lcpm <- rankGenes(as.matrix(lcpm))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  train_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c("#666666", "#0a75ad","#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore Scores (Training Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## aLT -- testing set.
ranked_lcpm <- rankGenes(as.matrix(lcpm_test))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  test_metadata[, c("SampleID", "TMM")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase-Amplified",
  "Telomerase-NotAmplified"
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c( "#666666", "#0a75ad", "#f08080", "#ffc3a0")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "errorbar",
    width = 0.8,
    linewidth = 1
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[4], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Testing Set)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

## ALT -- ackerman set.
ranked_lcpm <- rankGenes(as.matrix(expr_ackerman))

# singscore
score_res <- simpleScore(
  rankData = ranked_lcpm,
  upSet = up_genes,
  centerScore = TRUE
)

# plotting dataframe
plot_df <- data.frame(
  SampleID = rownames(score_res),
  Score = score_res$TotalScore
)

# Joining metadata
plot_df <- left_join(
  plot_df,
  metadata_ackerman[, c("SampleID", "TMM_Category")],
  by = "SampleID"
)

# Plot
groups <- c(
  "NO_TMM",
  "ALT",
  "Telomerase"
  
)

plot_df$TMM <- factor(plot_df$TMM, levels = groups)
group_colors <- c( "#666666", "#0a75ad", "#f08080")
names(group_colors) <- groups

# Darker versions for points
darken <- function(color, factor = 0.7) {
  rgb_vals <- col2rgb(color) / 255
  rgb(
    rgb_vals[1] * factor,
    rgb_vals[2] * factor,
    rgb_vals[3] * factor
  )
}

point_colors <- sapply(group_colors, darken)

ggplot(
  plot_df,
  aes(
    x = TMM,
    y = Score,
    fill = TMM
  )
) +
  geom_violin(
    aes(color = TMM),
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.8
  ) +
  stat_summary(
    aes(color = TMM),
    fun = median,
    geom = "crossbar",
    width = 0.80,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = TMM),
    width = 0.15,
    size = 3,
    alpha = 0.8
  ) +
  stat_compare_means(
    comparisons = list(
      c(groups[1], groups[2]),
      c(groups[3], groups[2])
    ),
    method = "wilcox.test",
    step.increase = 0.10,
    size = 6
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = point_colors) +
  labs(
    x = "",
    y = "Singscore (Ackerman)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = "none"
  )

#################################################################################
#################################################################################
#################################################################################


cohen_res <- cohen.d(
  Score ~ RiskGroup,
  data = plot_df,
  hedges.correction = TRUE
)


d <- abs(cohen_res$estimate)

overlap_target_up3 <- 2 * pnom(-d/2)
