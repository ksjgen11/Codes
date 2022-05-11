setwd('d:/Projects/UCPH/Work/siRNA/')
library(tidyverse)
library(ggpubr)

# data load
data <- read_tsv('2021-11-03 LINC_MET_RNF43_PANC1 SiRNA_2nd_data.txt') %>%
  unique() %>%
  na.omit(.)

colnames(data) <- c('sample', 'target', 'ct_mean')

samples <- str_split_fixed(data$sample, '[0-9]$', 2)
data$sample <- samples[,1]

head(data)

# calculate delta ct
hk_gene <- data %>%
  filter(target == 'beta-actin')

deltact <- left_join(data, hk_gene, by = 'sample', suffix = c('', '.hk')) %>%
  select(-target.hk) %>%
  mutate(delta_ct = ct_mean - ct_mean.hk)

# calculate delta delta ct and  foldchange, mean and sd
nt_data <- deltact %>%
  filter(sample == 'NT') %>%
  group_by(target) %>%
  summarize(mean_of_dct_nt = mean(delta_ct))

ddct <- left_join(deltact, nt_data, by = 'target')  %>%
  mutate(ddct = delta_ct - mean_of_dct_nt, Foldchange = 2^(-ddct)) %>%
  na.omit(.) %>%
  filter(target == 'LINC' | target == 'MET')

mean_sd <- ddct %>%
  group_by(target, sample) %>%
  summarize(mean_of_Fc = mean(Foldchange), sd_of_Fc = sd(Foldchange)) %>%
  filter(target == 'LINC' | target == 'MET')

my_comparisons = list(c('LINC', 'LIP'),c('LINC','MET'),c('LINC', 'NT'),c('LINC', 'WT'), 
                      c('LIP', 'MET'),c('LIP', 'NT'), c('LIP', 'WT'), 
                      c('MET', 'NT'), c('MET', 'WT'), c('NT', 'WT'))


# draw a plot



ggplot(ddct, aes(x = sample, y = Foldchange, fill = target)) +
  geom_boxplot() +
  #facet_wrap(~target) +
  #scale_y_continuous(limits = c(-1, 10)) +
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif', method = 't.test')

qpcr_mean_plot <- ggplot(mean_sd, aes(x = sample, y = mean_of_Fc, fill = target)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~target) +
  geom_errorbar(aes(ymin = mean_of_Fc - sd_of_Fc, ymax = mean_of_Fc + sd_of_Fc)) +
  theme_bw()

ggsave(qpcr_mean_plot, file = 'qpcr_mean_plot.png')
