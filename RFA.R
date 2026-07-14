
library(tidyverse) 
library(Matching)  
library(RFA)

ecology <- read.table('21_ecology_factor.txt', header = TRUE, row.names = NULL)
geography <- read.table('21_geography_factor.txt', header = TRUE, row.names = NULL)
phy_spe <- read.table('21_phylo.spe_factor.txt', header = TRUE, row.names = NULL)
strategy <- read.table('21_metabolism_rDsr.txt', header = TRUE, row.names = NULL)
combine <- merge(ecology,geography,by = 'accession')
combine <- merge(combine,phy_spe,by = 'accession')
combine <- merge(combine,strategy,by = 'accession')
combine <- combine[,-1]
combine[, c(1,11)] <- lapply(combine[, c(1,11)], as.factor)


rfa_fit <- rfa(
  formula = strategy ~ Oxygen, # treatment variable
  covariates = ~ spe_PC1 + spe_PC2 + spe_PC3, # confounders
  data = combine,
  classification = T
)
rfa_fit$fit
rfa_fit$yrf
rfa_fit$xrf
head(rfa_fit$data)
summary_rfa(rfa_fit) %>%
  mutate_if(is.numeric, function(x) round(x, 3)) 
plot_rfa(rfa_fit)
plot_rfa(rfa_fit, varname = "Phone Call") +
  labs(
    x = "Estimate",
    y = NULL,
    title = "Effect of GOT phone calls on voter turnout"
  ) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_test()

