library(dplyr)
library(forcats)
library(ggplot2)
library(readr)
library(writexl)

spatstat <- read_tsv('data/0.025.tsv', col_types = cols_only(
    ID = col_character(),
    cell_type = col_factor(),
    prop_close = col_double(),
    n_close_per_panCK = col_double()))

spatstat <- spatstat %>%
    dplyr::filter(cell_type != 'immune_cell')

clin <- read_tsv('data/clin_DBL_v_str_dens.tsv', col_types = cols_only(
    ID = col_character(),
    Cascon = col_factor(levels=c('0', '1'))))

clinst <- left_join(spatstat, clin, by = 'ID') %>%
    mutate(Cascon=fct_recode(Cascon, control='0', case='1'))

pdf('plots/boxplot_spatstat.pdf')

ggplot(clinst, aes_string(x='Cascon', y='n_close_per_panCK')) +
    geom_boxplot(outlier.alpha=0.0) +
    geom_jitter(width = 0.2, height = 0, shape=1, size=0.7) +
    facet_wrap(vars(cell_type), scales = "free_y", ncol=4) +
    theme_classic()

dev.off()


##############################
# Compute summary statistics #
##############################

summary_stats <- clinst %>%
  group_by(cell_type) %>%
  summarize(
    median = median(n_close_per_panCK, na.rm = T),
    p25 = quantile(n_close_per_panCK, probs = .25, na.rm = T),
    p75 = quantile(n_close_per_panCK, probs = .75, na.rm = T),
    n_obs = n()
  )

write_xlsx(summary_stats, 'results/summary_stats_spatstat.xlsx')
