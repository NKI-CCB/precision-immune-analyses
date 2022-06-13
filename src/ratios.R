library(conflicted)

library(broom)
library(coin, quietly = T)
library(dplyr)
library(forcats)
library(ggplot2)
library(purrr)
library(RColorBrewer)
library(readr)
library(rprojroot)
library(tibble)
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

cell_type_labels_vectra <- c(
    lymphocytes = "Lymphocytes",
    all_t_cells = "All T-cells",
    `CD20+` = "CD20+ B-cells",
    `CD3+_CD8-` = "Helper T-cells",
    `CD3+_CD8+` = "CD3+CD8+ T-cells",
    `CD3+_FOXP3+` = "CD3+FOXP3+ T-cells",
    `CD68+` = "CD68+ cells")
cell_type_labels_ki67 <- c(
    `CD8+_Ki67+` = "CD8+Ki67+ T-cells")
cell_type_labels <- c(cell_type_labels_vectra, cell_type_labels_ki67)

parse_args <- function(args) {
    arg_names <- c('cell_density', 'area', 'clinical', 'subset', 'stats',
                   'tests')
    stopifnot(length(args) == length(arg_names))
    args <- as.list(args)
    names(args) <- arg_names
    args
}

args <- parse_args(commandArgs(T))

ki67_cd8_var <- paste0('density_', args$area, '_CD8._Ki67.')

#############################
# Read and prepare the data #
#############################


# Read and prepare densities from Vectra #

na_as_zero <- function (x) {
    x[is.na(x)] <- 0
    x
}

density_vectra <- args$cell_density %>%
    read_tsv(col_types = cols(
        .default = col_double(),
        batch = col_factor(),
        ID = col_character())) %>%
    select(-batch)  %>%
    mutate(
        all_t_cells = na_as_zero(`CD3+_CD8+`) + na_as_zero(`CD3+_CD8-`) + na_as_zero(`CD3+_FOXP3+`),
        lymphocytes = all_t_cells + na_as_zero(`CD20+`))  %>%
    pivot_longer(-ID, names_to = "cell_type", values_to = "density") %>%
    dplyr::filter(cell_type != 'panCK+', cell_type != 'Other', !is.na(density))

# Read and prepare clinical variables #

clin_boxplot_vars <- c('Cascon', 'ER_status', 'grade', 'Her2status', 'COX2_status', 'fibrosis_yn') 
clin_vars_cat <- c(clin_boxplot_vars,
                   'PR_status', 'clin_pres', 'Yearsto_iIBC_cat_cases', 'diameter_cat',
                   'margin', 'dom_growthpat', 'necrosis', 'calcs')                                                          
clin_col_spec <- c(structure(map(clin_vars_cat, ~ col_factor()), names=clin_vars_cat))
clin_col_spec[['ki67perc_t']] = col_double()
clin_col_spec[['ID']] <- col_character()
clin_col_spec[['Cascon']] <- col_factor(levels=c('0', '1'))
clin_col_spec[['fibrosis_yn']] <- col_factor(levels=c('0', '1'))
clin_col_spec[['grade']] <- col_factor(levels=c('1', '2', '3'))
clin_col_spec[[ki67_cd8_var]] <- col_double()
clin_col_spec[['Subtype_10']] <- col_integer()

clinical <- args$clinical %>%
    read_tsv(col_types = do.call(cols_only, clin_col_spec)) %>%
    mutate(
        Cascon=fct_recode(Cascon, control='0', case='1'),
        fibrosis_yn=fct_recode(fibrosis_yn, absent='0', present='1'),
        ki67_cat = factor(ki67perc_t >= 14,
                          levels = c(F, T),
                          labels = c('<14', '>=14'))) %>%
    select(-ki67perc_t)
clin_vars_cat <- c(clin_vars_cat, 'ki67_cat')
clin_boxplot_vars <- c(clin_boxplot_vars, 'ki67_cat')

if (args$subset == 'all') {
} else if (args$subset == 'erpos_her2neg') {
  clinical <- dplyr::filter(clinical, Subtype_10 == 1)
} else if (args$subset == 'fibrosis') {
  clinical <- dplyr::filter(clinical, fibrosis_yn == 'absent')
} else {
  stop('No such subset')
}

clin_cat <- clinical %>%
    select(ID, Cascon, all_of(clin_vars_cat)) %>%
    pivot_longer(all_of(clin_vars_cat),
                 names_to = "clinical_variable",
                 values_to = "clinical_value",
                 values_drop_na = TRUE)

## One cell type's data lives in the clinical data frame

density_ki67 <- clinical %>%
    select(ID, all_of(ki67_cd8_var)) %>%
    rename(density = !!ki67_cd8_var) %>%
    mutate(cell_type = 'CD8+_Ki67+') %>%
    dplyr::filter(!is.na(density))

density <- bind_rows(density_vectra, density_ki67) %>%
    mutate(cell_type = factor(cell_type, levels = names(cell_type_labels))) %>%
    dplyr::filter(ID %in% clinical$ID)

stopifnot(!is.na(density$cell_type))

cell_type_combinations <- levels(density$cell_type) %>%
  combn(2) %>%
  array_tree(c(2, 1)) %>%
  map_dfr(~ tibble(cell_type1=.[[1]], cell_type2=.[[2]]))
density_ratios <- cell_type_combinations %>%
  full_join(density, by=c(cell_type1='cell_type')) %>%
  full_join(density, by=c(cell_type2='cell_type', ID='ID'), suffix=c('1', '2')) %>%
  mutate(ratio=density1/density2) %>%
  dplyr::filter(!is.na(ratio))

clin_ratios <- left_join(clin_cat, density_ratios, by = 'ID')

#########################
# Test for significance #
#########################

# Perform a non-parametric test, automaticly choosing wilcox or kruskal based on the number of the
# number of levels in the variable.
non_par_test <- function(ratio, clinical_value) {
  clinical_value <- fct_drop(clinical_value)
  if (length(levels(clinical_value)) < 2) {
    tibble(
      p.value = NA,
      statistic = NA,
      method = 'only one level')
  } else if (length(levels(clinical_value)) == 2) {
    tidy(coin::wilcox_test(ratio ~ clinical_value))
  } else {
    tidy(coin::kruskal_test(ratio ~ clinical_value))
  }
}

tests_cat <- clin_ratios %>%
    group_by(cell_type1, cell_type2, clinical_variable) %>%
    summarise(test = non_par_test(ratio, clinical_value),
              n_obs = n()) %>% 
    unpack(test) %>%
    rename(nominal_p = p.value)

fdr_df <- tests_cat %>%
  ungroup() %>%
  dplyr::filter(clinical_variable == 'Cascon') %>%
  mutate(fdr =  p.adjust(nominal_p, 'BH')) %>%
  select(clinical_variable, cell_type1, cell_type2, fdr)
tests <- tests_cat %>%
  full_join(fdr_df, by=c('clinical_variable', 'cell_type1', 'cell_type2')) %>%
  arrange(cell_type1, cell_type2, clinical_variable)
write_xlsx(tests, args$tests)


##############################
# Compute summary statistics #
##############################

summary_stats_by_clin_var <- clin_ratios %>%
  group_by(cell_type1, cell_type2, clinical_variable, clinical_value) %>%
  summarize(
    median = median(ratio, na.rm = T),
    p25 = quantile(ratio, probs = .25, na.rm = T),
    p75 = quantile(ratio, probs = .75, na.rm = T),
    n_obs = n()
  )

summary_stats_overall <- density_ratios %>%
  group_by(cell_type1, cell_type2) %>%
  summarize(
    median = median(ratio, na.rm = T),
    p25 = quantile(ratio, probs = .25, na.rm = T),
    p75 = quantile(ratio, probs = .75, na.rm = T),
    n_obs = n()
  )

summary_stats <- bind_rows(
  summary_stats_overall,
  summary_stats_by_clin_var)

write_xlsx(summary_stats, args$stats)
