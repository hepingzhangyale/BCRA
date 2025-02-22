rm(list = ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd('/home/wd278/BCRA')
path_to_data <- "./demo/reproduce_figure_tables/data_to_replicate_figure_tables/"

## Figure 1: UKB P-values
load(path_to_data, 'pval_UKB.RData')
all_pvals_sp %>% group_by(chr, blk, iter, super_var) %>%
  summarise(p = min(min(constant, probability), chisquared)) %>%
  mutate(p = ifelse(p < 1e-6, 1e-6, p)) %>%
  ggplot(aes(x = iter, y = -log10(p), color = super_var)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05/2723), linetype = 'dashed', color = 'grey40') +
  facet_wrap(~super_var, nrow = 5) +
  theme_bw()


