library(tidyverse)
library(ggpubr)

growth <- read_csv('growth_2.csv') %>%
  group_by(Treatment, Elapsed, Concentration) %>%
  summarize(mean = mean(growth), sd = sd(growth))
growth$Treatment <-factor(growth$Treatment, levels = c('Panc1', 'Lip', 'NT', 'LINC', 'MET'))



confluency <- read_tsv('cellconfluency_92hr.txt')
confluency$Treatment <-factor(confluency$Treatment, levels = c('Panc1', 'Lip', 'NT', 'LINC', 'MET'))

my_comparisons = list(c('Panc1', 'Lip'),c('Panc1', 'NT'), c('Panc1', 'LINC'), c('Panc1','MET'),
                      c('Lip', 'NT'), c('Lip', 'LINC'), c('Lip', 'MET'), 
                      c('NT', 'LINC'), c('NT', 'MET'), 
                      c('LINC', 'MET') )

confluent <- ggplot(confluency, aes(x = Treatment, y = Confluency)) +
  geom_boxplot()  + facet_wrap(~Concentration) +
  theme_bw()  #+ stat_compare_means(comparisons = my_comparisons, label = 'p.signif', method = 'wilcox')

confluent

confluent2 <- ggplot(confluency, aes(x = as.factor(Concentration), y = Confluency)) +
  geom_boxplot()  + facet_wrap(~Treatment) +
  theme_bw()

confluent2

ggsave('Confluency_test.pdf', plot = confluent)
ggsave('Confluency_dedendent_concentration.pdf', plot = confluent2)

growth_plot <- ggplot(growth, aes(x = Elapsed, y = mean, col = Treatment)) +
  geom_line() + facet_wrap(~Concentration) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3) +
  theme_bw()

growth_plot2 <- ggplot(growth, aes(x = Elapsed, y = mean, col = Treatment)) +
  geom_line() + facet_wrap(~Concentration) +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  theme_bw()


growth_plot
growth_plot2
ggsave('Growth.pdf', plot = growth_plot)
ggsave('Growth_without_errorbar.pdf', plot = growth_plot2)
