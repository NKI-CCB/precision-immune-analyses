library(conflicted)

library(broom)
library(coin, quietly = T)
library(dplyr)
library(forcats)
library(purrr)
library(RColorBrewer)
library(readr)
library(rprojroot)
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
    arg_names <- c('cell_density', 'cell_density_ki67', 'area', 'clinical', 'subset', 'tests')
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
        t_number = col_character())) %>%
    select(-batch)  %>%
    mutate(
        all_t_cells = na_as_zero(`CD3+_CD8+`) + na_as_zero(`CD3+_CD8-`) + na_as_zero(`CD3+_FOXP3+`),
        lymphocytes = all_t_cells + na_as_zero(`CD20+`))  %>%
    pivot_longer(-t_number, names_to = "cell_type", values_to = "density") %>%
    dplyr::filter(cell_type != 'panCK+', cell_type != 'Other', !is.na(density))

# Read and prepare clinical variables #

clin_boxplot_vars <- c('Cascon', 'ER_status', 'grade', 'Her2status', 'COX2_status', 'fibrosis_yn') 
clin_vars_cat <- c(clin_boxplot_vars,
                   'age_cat1', 'PR_status', 'clin_pres', 'Yearsto_iIBC_cat_cases', 'diameter_cat',
                   'margin', 'dom_growthpat', 'necrosis', 'calcs')                                                          
clin_col_spec <- c(structure(map(clin_vars_cat, ~ col_factor()), names=clin_vars_cat))
clin_col_spec[['ki67perc_t']] = col_double()
clin_col_spec[['t_number']] <- col_character()
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

## One cell type's data lives in the clinical data frame

density_ki67 <- clinical %>%
    select(t_number, all_of(ki67_cd8_var)) %>%
    rename(density = !!ki67_cd8_var) %>%
    mutate(cell_type = 'CD8+_Ki67+') %>%
    dplyr::filter(!is.na(density))

density <- bind_rows(density_vectra, density_ki67) %>%
    mutate(cell_type = factor(cell_type, levels = names(cell_type_labels))) %>%
    dplyr::filter(t_number %in% clinical$t_number)

stopifnot(!is.na(density$cell_type))

##

clin_dens <- left_join(clinical, density, by = 't_number') %>%
  group_by(cell_type) %>%
  mutate(
    grade_dens = factor(case_when(
      grade == 1 ~ 'low',
      grade == 3 ~ 'high',
      grade == 2 & density <= median(density) ~ 'low',
      grade == 2 & density > median(density) ~ 'high',
      TRUE ~ 'ERROR'), levels = c('low', 'high'))) %>%
  ungroup() %>%
  select(t_number, cell_type, Cascon, grade_dens)

tests <- clin_dens %>%
  group_by(cell_type) %>%
  summarise(
    test = tidy(coin::chisq_test(Cascon ~ grade_dens)),
    n_control_low = sum(Cascon == 'control' & grade_dens == 'low'),
    n_control_high = sum(Cascon == 'control' & grade_dens == 'high'),
    n_case_low = sum(Cascon == 'case' & grade_dens == 'low'),
    n_case_high = sum(Cascon == 'case' & grade_dens == 'high')) %>%
  unpack(test) %>%
  rename(nominal_p = p.value) %>%
  ungroup() %>%
  mutate(fdr =  p.adjust(nominal_p, 'BH')) %>%
  arrange(cell_type) %>%
  write_xlsx(args$tests)
