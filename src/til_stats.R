library(conflicted)

library(broom)
library(coin, quietly = T)
library(dplyr)
library(forcats)
library(purrr)
library(RColorBrewer)
library(readr)
library(readxl)
library(tibble)
library(tidyr)
library(writexl)

options(dplyr.summarise.inform = F)
options(warn = 2)

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
    arg_names <- c('cell_density', 'area', 'clinical', 'subset', 
                   'til_scores', 'tests')
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
                   'age_cat1', 'PR_status', 'clin_pres', 'Yearsto_iIBC_cat_cases', 'diameter_cat',
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

til_vars <- c('avgTIL', 'avgTIL_cat', 'avgTIL_cat_1', 'avgTIL_cat_5')
til_scores <- read_tsv(args$til_scores, col_types=cols_only(
    ID = col_character(),
    avgTIL = col_double(),
    avgTIL_cat = col_factor(levels = c('0', '1', '2')),
    avgTIL_cat_1 = col_factor(levels = c('0', '1')),
    avgTIL_cat_5 = col_factor(levels = c('0', '1'))))
clinical <- left_join(clinical, til_scores, by='ID')

if (args$subset == 'all') {
} else if (args$subset == 'erpos_her2neg') {
  clinical <- dplyr::filter(clinical, Subtype_10 == 1)
} else if (args$subset == 'fibrosis') {
  clinical <- dplyr::filter(clinical, fibrosis_yn == 'absent')
} else {
  stop('No such subset')
}

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

clin_dens <- left_join(clinical, density, by = 'ID')

tests_cor <- clin_dens %>%
  select(ID, cell_type, density, avgTIL) %>%
  pivot_longer(avgTIL, names_to = "variable", values_to = "value") %>%
  group_by(cell_type, variable) %>%
  summarise(
    test = tidy(coin::spearman_test(value ~ density)),
    estimate = cor(value, density, method = 'spearman'),
    n=n()) %>%
  unpack(test) %>%
  rename(nominal_p = p.value) %>%
  ungroup() %>%
  mutate(fdr =  p.adjust(nominal_p, 'BH')) %>%
  arrange(cell_type)


cascon_test <- function(cascon, y) {
  stopifnot(length(cascon) == length(y))
  excluded <- is.na(cascon) | is.na(y)
  cascon <- cascon[!excluded]
  y <- y[!excluded]
  stopifnot(is.factor(cascon))
  stopifnot(levels(cascon) == c('control', 'case'))
  if (is.double(y)) {
    cascon <- fct_drop(cascon)
    if (length(levels(cascon)) == 2) {
      test <- tidy(coin::wilcox_test(y ~ cascon))
    } else {
      test <- tidy(coin::kruskal_test(y ~ cascon))
    }
  } else if (is.factor(y)) {
    cascon <- fct_drop(cascon)
    y <- fct_drop(y)
    test <- tidy(coin::chisq_test(cascon ~ y))
  } else {
    str(cascon)
    str(y)
    stop("Variable types not supported.")
  }
  test$n_obs <- length(cascon)
  test$n_control <- sum(cascon == 'control')
  test$n_case <- sum(cascon == 'case')
  rename(test, nominal_p = p.value)
}

tests_surv <- clin_dens %>%
  select(ID, Cascon, all_of(til_vars)) %>%
  distinct() %>%
  summarise(across(all_of(til_vars), ~ cascon_test(Cascon, .))) %>%
  pivot_longer(all_of(til_vars), names_to='til_variable',  values_to='test') %>%
  unpack(test) %>%
  ungroup() %>%
  mutate(fdr =  p.adjust(nominal_p, 'BH'))

write_xlsx(list(density_correlation=tests_cor, case_control=tests_surv, density), args$tests)
