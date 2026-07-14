library(ape)
library(tidyverse)
library(reshape2)
library(vegan)
library(philentropy)

#input data
ecology <- read.table('21_ecology_factor.txt', header = TRUE, row.names = NULL)
geography <- read.table('21_geography_factor.txt', header = TRUE, row.names = NULL)
phy_spe <- read.table('21_phylo.spe_factor.txt', header = TRUE, row.names = NULL)
strategy <- read.table('21_metabolism_rDsr.txt', header = TRUE, row.names = NULL)
combine <- merge(ecology,geography,by = 'accession')
combine <- merge(combine,phy_spe,by = 'accession')
combine <- merge(combine,strategy,by = 'accession')

RF_result <- read.csv('21_RF_result.csv')
RF_result$X <- factor(RF_result$X,levels = c('spe_PC1','spe_PC2','spe_PC3',
                                             'geo_PC1','geo_PC2','geo_PC3',
                                             'Oxygen','Temperature','pH','Salinity'))
create_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""   
  )
}

RF_result$stars <- create_stars(RF_result$adjusted_p_BH)
library(ggplot2)
ggplot(RF_result, aes(x = X, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#bebdbd", width = 0.5) +
  # 添加显著性符号（调整 nudge_x 控制位置�?
  geom_text(
    aes(label = stars), 
    nudge_x = 0,          # 横向位置微调（根据数据范围调整）
    nudge_y = max(RF_result$MeanDecreaseGini) * 0.05,  # 纵向位置
    hjust = -0.2,         # 文字左对齐（在条形右侧）
    size = 5,             # 符号大小
    color = "black"
  ) +
  scale_y_continuous(limit = c(0,400))+
  coord_flip() +          # 翻转坐标�?
  labs(title = NULL, x = NULL, y = 'MeanDecreaseGini') +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black")
  ) +
  ggh4x::force_panelsizes(rows = unit(18, "cm"), cols = unit(8, "cm")) +

  scale_x_discrete(expand = expansion(add = c(0.5, 0.5)))

df <- combine %>%
  group_by(Oxygen, strategy) %>%
  summarise(count = n(), .groups = "drop")

df$Oxygen <- factor(df$Oxygen, levels = c('aerobe','facultative','anaerobe'))
df$strategy <- factor(df$strategy, levels = c('others','rDsr','soxCD','both'))

# 绘制堆积柱状�?
ggplot(df, aes(x = Oxygen, y = count, fill = strategy)) +
  scale_fill_manual(values = c('#bebebe','#4f9ace','#d24a5a','#5866a2'))+
  geom_bar(stat = "identity", position = "fill", width = 0.5) +
  labs(x = NULL, y = "Proportion", fill = "Variable") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 14), 
    axis.text.y = element_text(color = "black", size = 14) 
  )+
  ggh4x::force_panelsizes(rows = unit(8, "cm"), cols = unit(8, "cm"))



analyze_genotype_environment <- function(data) {
  contingency_table <- table(data$genotype, data$environment)
  chi_test <- chisq.test(contingency_table)
  std_residuals <- chi_test$stdres
  p_values <- 2 * pnorm(abs(std_residuals), lower.tail = FALSE)
  adj_p_values <- p.adjust(p_values, method = "BH")
  residuals_df <- as.data.frame(as.table(std_residuals))
  names(residuals_df) <- c("Genotype", "Environment", "Residual")
  residuals_df$Significance <- ifelse(
    adj_p_values < 0.001, "***",
    ifelse(adj_p_values < 0.01, "**",
           ifelse(adj_p_values < 0.05, "*", ""))
  )
  heatmap_plot <- ggplot(residuals_df, 
                         aes(x = Environment, y = Genotype, fill = Residual)) +
    geom_tile(color = "gray90") +
    scale_fill_gradient2(
      low = "#4f9ace", 
      mid = "white",
      high = "#d24a5a", 
      midpoint = 0,
      limits = c(-max(abs(residuals_df$Residual)), max(abs(residuals_df$Residual)))
    ) +
    geom_text(aes(label = sprintf("%.2f", Residual)), 
              color = "black", size = 4, vjust = 1) +
    geom_text(aes(label = Significance), 
              color = "black", size = 6, vjust = -0.5) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(fill = "Standardized\nResidual",
         x = NULL, y = NULL) +
    coord_fixed()
  results <- list(
    contingency_table = contingency_table,
    chi_test = chi_test,
    std_residuals = std_residuals,
    adj_p_values = matrix(adj_p_values, 
                          nrow = nrow(std_residuals),
                          dimnames = dimnames(std_residuals)),
    plot = heatmap_plot
  )
  
  cat("\nChi-square test results:\n")
  print(chi_test)
  
  cat("\nContingency table:\n")
  print(contingency_table)
  
  cat("\nStandardized residuals with BH-adjusted p-values:\n")
  print(cbind(
    as.data.frame.table(std_residuals),
    Adj.P.value = matrix(adj_p_values, ncol = 1)
  ))
  
  print(heatmap_plot)
  
  return(invisible(results))
}
df <- combine[combine$strategy != 'others', ]
df <- df[, c(2, 12)]
colnames(df) <- c("environment", "genotype")
df$environment <- factor(df$environment, levels = c('aerobe', 'facultative', 'anaerobe'))
df$genotype <- factor(df$genotype, levels = c('soxCD', 'rDsr', 'both'))

results <- analyze_genotype_environment(df)


tree <- read.tree('bac120_r214.tree')
all <- tree[["tip.label"]]
id <- all[all %in% combine$accession]
tree2 <- keep.tip(tree, tip=id)
rm(tree)
phylo_dist <- cophenetic(tree2)
nmds <- metaMDS(phylo_dist, k = 2, trymax = 100)

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$accession <- rownames(nmds_scores)
nmds_scores <- merge(nmds_scores,combine,by = 'accession')

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = strategy)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x = "NMDS1",y = "NMDS2") +
  scale_x_continuous(limit = c(-1.5,2))+
  scale_y_continuous(limit = c(-2,1.2))+
  scale_color_manual(values = c('#5866a2','#bebebe', '#4f9ace', '#d24a5a')) +  
  annotate("text", x = Inf, y = Inf, 
           label = paste("Stress =", round(nmds$stress, 3)),
           hjust = 1.1, vjust = 1.1, size = 4)+
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),  
    axis.text.x = element_text(color = "black", size = 14), 
    axis.text.y = element_text(color = "black", size = 14)  
  )+
  ggh4x::force_panelsizes(rows = unit(7, "cm"), cols = unit(8, "cm"))

ggplot(nmds_scores[nmds_scores$strategy %in% c('others','both') == F, ], aes(x = NMDS1))+
  geom_density(
    aes(fill = strategy),
    alpha = 0.7,
    color = 'white',    
    linewidth = 0.1
  )+
  scale_x_continuous(limit = c(-1.5,2))+
  scale_fill_manual(values = c('#4f9ace', '#d24a5a')) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),  
    axis.text.x = element_text(color = "black", size = 14), 
    axis.text.y = element_text(color = "black", size = 14)  
  ) +
  ggh4x::force_panelsizes(rows = unit(1, "cm"), cols = unit(8, "cm"))

ggplot(nmds_scores[nmds_scores$strategy %in% c('others','both') == F, ], aes(x = NMDS2))+
  geom_density(
    aes(fill = strategy),
    alpha = 0.7,
    color = 'white',    
    linewidth = 0.1
  )+
  scale_x_continuous(limit = c(-2,1.2))+
  scale_fill_manual(values = c('#4f9ace', '#d24a5a')) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),  
    axis.text.x = element_text(color = "black", size = 14), 
    axis.text.y = element_text(color = "black", size = 14)  
  ) +
  ggh4x::force_panelsizes(rows = unit(1, "cm"), cols = unit(7, "cm"))

result <- adonis2(phylo_dist ~ strategy, data = combine, permutations = 999)
print(result)
library(pairwiseAdonis)
samples_order <- rownames(phylo_dist)
combine_ordered <- combine[match(samples_order, combine$accession), ]
group <- factor(combine$strategy)
pairwise_result <- pairwise.adonis(phylo_dist, group, p.adjust.m = "BH")
print(pairwise_result)




#dsr
tree <- read.tree('30_con_dsrAB_trim.best.nwk')
all <- tree[["tip.label"]]
id <- all[all %in% combine$accession]
tree2 <- keep.tip(tree, tip=id)
rm(tree)
rdsr_phylo_dist <- cophenetic(tree2)
#soxCD
tree <- read.tree('31_con_soxCD_trim.best.nwk')
all <- tree[["tip.label"]]
id <- all[all %in% combine$accession]
tree2 <- keep.tip(tree, tip=id)
rm(tree)
soxCD_phylo_dist <- cophenetic(tree2)
#soxBYZ
tree <- read.tree('32_con_soxBYZ_trim.best.nwk')
all <- tree[["tip.label"]]
id <- all[all %in% combine$accession]
tree2 <- keep.tip(tree, tip=id)
rm(tree)
soxB_phylo_dist <- cophenetic(tree2)
#rdsr soxB
common_rows <- intersect(rownames(rdsr_phylo_dist), rownames(soxB_phylo_dist)) 
common_cols <- intersect(colnames(rdsr_phylo_dist), colnames(soxB_phylo_dist)) 
soxB_phylo_dist_filtered <- soxB_phylo_dist[common_rows, common_cols, drop = FALSE]
df <- data.frame(x = as.vector(soxB_phylo_dist_filtered),y = as.vector(rdsr_phylo_dist))
model <- lm(y ~ x + 0, data = df)  
summary(model)
model_summary <- summary(model)
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients["x", "Pr(>|t|)"]  
label <- sprintf(
  "R² = %.2f\nP %s",
  r_squared,
  ifelse(p_value < 0.001, "< 0.001", format.pval(p_value, eps = 0.001))
)

ggplot(df, aes(x = log10(x+1), y = log10(y+1))) +
  geom_hex(bins = 100) +
  geom_smooth(method = "lm", formula = y ~ x + 0, se = T, size = 1, 
              color = '#c1272d', fullrange = T, alpha = 0.15) +
  scale_fill_gradientn(
    colours = c("#bebebe40", "#0071bcFF"),  
    values = scales::rescale(c(0,1))) +
  labs(x = 'Truncated Sox system', y = 'rDsrAB') +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label,
    hjust = -0.1,
    vjust = 1.5,
    size = 5,
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14)
  ) +
  ggh4x::force_panelsizes(rows = unit(8, "cm"), cols = unit(8, "cm"))

#soxCD soxBYZ
common_rows <- intersect(rownames(soxCD_phylo_dist), rownames(soxB_phylo_dist)) 
common_cols <- intersect(colnames(soxCD_phylo_dist), colnames(soxB_phylo_dist)) 
soxB_phylo_dist_filtered <- soxB_phylo_dist[common_rows, common_cols, drop = FALSE]
df <- data.frame(x = as.vector(soxB_phylo_dist_filtered),y = as.vector(soxCD_phylo_dist))
model <- lm(y ~ x + 0, data = df)  
summary(model)
model_summary <- summary(model)
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients["x", "Pr(>|t|)"]  
label <- sprintf(
  "R² = %.2f\nP %s",
  r_squared,
  ifelse(p_value < 0.001, "< 0.001", format.pval(p_value, eps = 0.001))
)

ggplot(df, aes(x = log10(x+1), y = log10(y+1))) +
  geom_hex(bins = 100) +
  geom_smooth(method = "lm", formula = y ~ x + 0, se = T, size = 1, 
              color = '#c1272d', fullrange = T, alpha = 0.15) +
  scale_fill_gradientn(
    colours = c("#bebebe40", "#0071bcFF"),  
    values = scales::rescale(c(0, 0.1, 0.3, 1))) +
  labs(x = 'Truncated Sox system', y = 'SoxCD') +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label,
    hjust = -0.1,
    vjust = 1.5,
    size = 5,
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14)
  ) +
  ggh4x::force_panelsizes(rows = unit(8, "cm"), cols = unit(8, "cm"))


