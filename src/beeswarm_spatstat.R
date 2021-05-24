library(dplyr)
library(forcats)
library(ggplot2)
library(readr)

spatstat <- read_tsv('data/0.025.tsv', col_types = cols_only(
    sample = col_character(),
    cell_type = col_factor(),
    prop_close = col_double(),
    n_close_per_panCK = col_double()))

clin <- read_tsv('data/clin_DBL_v_str_dens.tsv', col_types = cols_only(
    t_number = col_character(),
    Cascon = col_factor(levels=c('0', '1'))))

clinst <- left_join(spatstat, clin, by = c(sample = 't_number')) %>%
    mutate(Cascon=fct_recode(Cascon, control='0', case='1'))

pdf('plots/boxplot_spatstat.pdf')

ggplot(clinst, aes_string(x='Cascon', y='n_close_per_panCK')) +
    geom_boxplot(outlier.alpha=0.0) +
    geom_jitter(width = 0.2, height = 0, shape=1, size=0.7) +
    facet_wrap(vars(cell_type), scales = "free_y", ncol=4)

dev.off()
