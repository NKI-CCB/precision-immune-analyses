library(conflicted)

library(broom)
library(coin, quietly = T)
library(dplyr)
library(forcats)
library(ggplot2)
library(glue)
library(purrr)
library(readr)
library(tidyr)
library(writexl)

options(dplyr.summarise.inform = F)

tidy.ScalarIndependenceTest <- function(x) {
    tibble::tibble(
        p.value = as.numeric(coin::pvalue(x)),
        statistic = coin::statistic(x),
        method = x@method,
        alternative = x@statistic@alternative)
}

tidy.QuadTypeIndependenceTest <- function(x) {
    tibble::tibble(
        p.value = as.numeric(coin::pvalue(x)),
        statistic = coin::statistic(x),
        method = x@method)
}

non_par_test <- function(x, y) {
  stopifnot(length(x) == length(y))
  excluded <- is.na(x) | is.na(y)
  x <- x[!excluded]
  y <- y[!excluded]
  if (is.double(x) & is.double(y)) {
    test <- tidy(coin::spearman_test(y ~ x))
    test$estimate <- cor(x, y, method = 'spearman')
  } else if (is.factor(x) & is.double(y)) {
    x <- fct_drop(x)
    if (length(levels(x)) == 2) {
      test <- tidy(coin::wilcox_test(y ~ x))
    } else {
      test <- tidy(coin::kruskal_test(y ~ x))
    }
  } else if (is.double(x) & is.factor(y)) {
    y <- fct_drop(y)
    if (length(levels(y)) == 2) {
      test <- tidy(coin::wilcox_test(x ~ y))
    } else {
      test <- tidy(coin::kruskal_test(x ~ y))
    }
  } else if (is.factor(x) & is.factor(y)) {
    x <- fct_drop(x)
    y <- fct_drop(y)
    test <- tidy(coin::chisq_test(x ~ y))
  } else {
    str(x)
    str(y)
    stop("Combination of variable types not supported.")
  }
  test$n_obs <- length(x)
  rename(test, nominal_p = p.value)
}

density_vars <- c("lymphocytes", "all_t_cells", "CD20+", "CD3+_CD8-",
                  "CD3+_CD8+", "CD3+_FOXP3+_NA", "CD68+", "density_Tissue_CD8._Ki67.")
clin_vars_cat <- c('Cascon', 'ER_status', 'grade', 'Her2status', 'COX2_status', 'fibrosis_yn',
                   'age_cat1', 'PR_status', 'clin_pres', 'Yearsto_iIBC_cat_cases', 'diameter_cat',
                   'margin', 'dom_growthpat', 'necrosis', 'calcs')
clin_vars_cont <- c('ki67perc_t', 'ag_zone_t', 'TLS_GC_t')

clin_vars <- c(clin_vars_cat, clin_vars_cont)

col_spec <- c(
    structure(map(density_vars, ~ col_double()), names=density_vars),
    structure(map(clin_vars_cat, ~ col_factor()), names=clin_vars_cat),
    structure(map(clin_vars_cont, ~ col_double()), names=clin_vars_cont))
col_spec[['ID']] <- col_character()
col_spec[['Cascon']] <- col_factor(levels=c('0', '1'))
col_spec[['fibrosis_yn']] <- col_factor(levels=c('0', '1'))
col_spec[['grade']] <- col_factor(levels=c('1', '2', '3'))

clin <- read_tsv(
        'data/clin_DBL_v_str_dens.tsv',
        col_types = do.call(cols_only, col_spec)) %>%
    mutate(
        Cascon=fct_recode(Cascon, control='0', case='1'),
        fibrosis_yn=fct_recode(fibrosis_yn, absent='0', present='1'),
        ki67_cat = factor(ki67perc_t >= 14, levels = c(F, T), labels = c('lt14', 'gte14')))

# Calculate per area
area_ki67 <- read_tsv('data/ki67_counts.tsv', col_types = cols_only(
        `Classifier Label` = col_character(),
        ID = col_character(),
        area = col_double())) %>%
    dplyr::filter(`Classifier Label` == 'Tissue') %>%
    rename(area_ki67 = area) %>%
    distinct()

area_vectra <- read_tsv('data/area_Stroma.tsv', col_types = cols_only(
        ID = col_character(),
        area = col_double())) %>%
    rename(area_vectra = area)

clin <- clin %>%
    left_join(area_ki67, by='ID') %>%
    left_join(area_vectra, by='ID') %>%
    mutate(
        ag_zone_t = ag_zone_t / area_vectra,
        TLS_GC_t = TLS_GC_t / area_ki67)

# Ki67 #

ki67_tests <- map_dfr(
    set_names(c('Cascon', density_vars)),
    ~ non_par_test(clin[[.]], clin[['ki67perc_t']]),
    .id = 'variable2')
ki67cat_tests <- map_dfr(
    set_names(c('Cascon', density_vars)),
    ~ non_par_test(clin[[.]], clin[['ki67_cat']]),
    .id = 'variable2')

bind_rows(list(ki67per_t = ki67_tests, ki67_cat = ki67cat_tests), .id = 'variable1') %>%
  write_xlsx('results/significance_ki67.xlsx')

# Immune cell aggregates #

ag_vars <- c('Cascon', 'COX2_status', 'Her2status', 'ER_status', 'fibrosis_yn', 'grade', 'ki67_cat')

ag_zone_tests <- map_dfr(
    set_names(ag_vars),
    ~ non_par_test(clin[[.]], clin[['ag_zone_t']]),
    .id = 'variable2')
TLS_GC_tests <- map_dfr(
    set_names(ag_vars),
    ~ non_par_test(clin[[.]], clin[['TLS_GC_t']]),
    .id = 'variable2')

bind_rows(list(ag_zone_t = ag_zone_tests, TLS_GC_t = TLS_GC_tests),
          .id = 'variable1') %>%
  write_xlsx('results/significance_ag_tls.xlsx')

pdf('plots/hist-ag.pdf')
hist(clin$ag_zone_t)
hist(clin$TLS_GC_t)
invisible(dev.off())

cat_ag <- function (ag) {
  med_ag <- median(ag, na.rm = T)
  labels <- c(glue('lte{med_ag}'), glue('gt{med_ag}'))
  factor(ag > med_ag, levels = c(F, T), labels = labels)
}

clin$ag_zone_cat <- cat_ag(clin$ag_zone_t)
clin$TLS_GC_cat <- cat_ag(clin$TLS_GC_t)

ag_zone_cat_tests <- map_dfr(
    set_names(ag_vars),
    ~ non_par_test(clin[[.]], clin[['ag_zone_cat']]),
    .id = 'variable2')
TLS_GC_cat_tests <- map_dfr(
    set_names(ag_vars),
    ~ non_par_test(clin[[.]], clin[['TLS_GC_cat']]),
    .id = 'variable2')

bind_rows(list(ag_zone_cat = ag_zone_cat_tests,
               TLS_GC_cat = TLS_GC_cat_tests),
          .id = 'variable1') %>%
  write_xlsx('results/significance_ag_tls_cat.xlsx')

pdf('plots/boxplot_tls.pdf')

dens_vars = c('TLS_GC_t', 'ag_zone_t')

clin_plot <- clin %>%
    select(all_of(dens_vars), all_of(ag_vars)) %>%
    pivot_longer(all_of(ag_vars),
                 names_to = "clinical_variable",
                 values_to = "clinical_value",
                 values_drop_na = TRUE) %>%
    pivot_longer(all_of(dens_vars),
                 names_to = "density_variable",
                 values_to = "density",
                 values_drop_na = TRUE) %>%
    mutate(clin_var_val = paste0(clinical_variable, '=', clinical_value) %>%
        factor(levels = c("COX2_status=High", "COX2_status=Low",
                          "ki67_cat=lt14", "ki67_cat=gte14",
                          "Her2status=Positive", "Her2status=Negative",
                          "ER_status=Positive", "ER_status=Negative",
                          "fibrosis_yn=present", "fibrosis_yn=absent",
                          "grade=3", "grade=2", "grade=1",
                          "Cascon=case", "Cascon=control")))


get_breaks <- function(lim) {
  print(lim)
  width = 10**floor(log10(lim[2]))
  print(width)
  last_break = width * floor(lim[2] / width)
  print(last_break)
  breaks_large = seq(width, last_break, width)
  breaks_small = seq(0, width - (width/10), width/2)
  c(breaks_small, breaks_large)
}

clin_plot %>%
    ggplot(aes(y=density, x=clin_var_val)) +
    geom_boxplot(outlier.alpha=0.0) +
    geom_jitter(aes(colour = clinical_variable), width = 0.2, height = 0, shape=1, size=0.7) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    scale_y_continuous(expand = expansion(c(0, 0.05), c(0, 0)), limits = c(0, NA),
                       trans = scales::sqrt_trans(),
                       breaks = get_breaks) +
    scale_x_discrete("") +
    scale_colour_brewer(palette='Dark2') +
    facet_wrap(vars(density_variable), scales = 'free_x', ncol=2) +
    guides(colour = 'none') +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing.x = unit(5, 'pt'),
          axis.ticks.y = element_blank())

dev.off()

##############################
# Compute summary statistics #
##############################

summary_stats <- clin_plot %>%
  group_by(density_variable, clinical_variable, clinical_value) %>%
  summarize(
    median = median(density, na.rm = T),
    p25 = quantile(density, probs = .25, na.rm = T),
    p75 = quantile(density, probs = .75, na.rm = T),
    n_obs = n()
  )

write_xlsx(summary_stats, 'results/summary_stats_ag_tls.xlsx')
