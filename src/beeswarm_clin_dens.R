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
    arg_names <- c('cell_density', 'cell_density_ki67', 'area', 'clinical', 'subset', 'plot',
                   'stats', 'tests')
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

clin_cat <- clinical %>%
    select(t_number, Cascon, all_of(clin_vars_cat)) %>%
    pivot_longer(all_of(clin_vars_cat),
                 names_to = "clinical_variable",
                 values_to = "clinical_value",
                 values_drop_na = TRUE)

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

clin_dens <- left_join(clin_cat, density, by = 't_number')


#########################
# Test for significance #
#########################

# Perform a non-parametric test, automaticly choosing wilcox or kruskal based on the number of the
# number of levels in the variable.
non_par_test <- function(density, clinical_value) {
  clinical_value <- fct_drop(clinical_value)
  if (length(levels(clinical_value)) == 1) {
    tibble(
      p.value = NA,
      statistic = NA,
      method = 'only one level')
  } else if (length(levels(clinical_value)) == 2) {
    tidy(coin::wilcox_test(density ~ clinical_value))
  } else {
    tidy(coin::kruskal_test(density ~ clinical_value))
  }
}

tests_cat <- clin_dens %>%
    group_by(cell_type, clinical_variable) %>%
    summarise(test = non_par_test(density, clinical_value),
              n_obs = n()) %>% 
    unpack(test) %>%
    rename(nominal_p = p.value)

fdr_df <- tests_cat %>%
  ungroup() %>%
  dplyr::filter(clinical_variable == 'Cascon') %>%
  mutate(fdr =  p.adjust(nominal_p, 'BH')) %>%
  select(clinical_variable, cell_type, fdr)
tests <- tests_cat %>%
  full_join(fdr_df, by=c('clinical_variable', 'cell_type')) %>%
  arrange(cell_type, clinical_variable)
write_xlsx(tests, args$tests)

##############################
# Compute summary statistics #
##############################

summary_stats_by_clin_var <- clin_dens %>%
  group_by(cell_type, clinical_variable, clinical_value) %>%
  summarize(
    median = median(density, na.rm = T),
    p25 = quantile(density, probs = .25, na.rm = T),
    p75 = quantile(density, probs = .75, na.rm = T),
    n_obs = n()
  )

summary_stats_overall <- density %>%
  group_by(cell_type) %>%
  summarize(
    median = median(density, na.rm = T),
    p25 = quantile(density, probs = .25, na.rm = T),
    p75 = quantile(density, probs = .75, na.rm = T),
    n_obs = n()
  )

summary_stats <- bind_rows(
  summary_stats_overall,
  summary_stats_by_clin_var)

write_xlsx(summary_stats, args$stats)

############
# Plotting #
############


# Prepare Data #

# clin_var_val is a variable for the x axis and needs to be in a specific order.
clin_density_plot <- clin_dens %>%
    dplyr::filter(clinical_variable %in% clin_boxplot_vars) %>%
    mutate(clin_var_val = paste0(clinical_variable, '=', clinical_value) %>%
        factor(levels = c("COX2_status=High", "COX2_status=Low",
                          "ki67_cat=>=14", "ki67_cat=<14",
                          "Her2status=Positive", "Her2status=Negative",
                          "ER_status=Positive", "ER_status=Negative",
                          "fibrosis_yn=present", "fibrosis_yn=absent",
                          "grade=3", "grade=2", "grade=1",
                          "Cascon=case", "Cascon=control")))

# We need to to know maximum density for significance star placement
dens_summary <- clin_dens %>%
    group_by(cell_type) %>%
    summarize(max_dens=max(density, na.rm=T))

significance_annot <- left_join(
    transmute(tests_cat, cell_type, clinical_variable, nominal_p = nominal_p) %>%
        dplyr::filter(clinical_variable %in% clin_boxplot_vars),
    clin_density_plot %>%
        select(clinical_variable, clin_var_val) %>%
        distinct() %>%
        group_by(clinical_variable) %>%
        summarize(clin_var_val = median(unclass(clin_var_val))),
    by = 'clinical_variable')
significance_annot <- significance_annot %>%
    left_join(dens_summary, by = 'cell_type') %>%
    mutate(stars = case_when(
        nominal_p < 0.001 ~ '***',
        nominal_p < 0.01 ~ '**',
        nominal_p < 0.05 ~ '*'))

# Make the plot #

pdf(args$plot, 7, 10)
clin_density_plot %>%
    dplyr::filter(!is.na(density)) %>%
    ggplot(aes(y=density, x=clin_var_val)) +
    geom_boxplot(outlier.alpha=0.0) +
    geom_jitter(aes(colour = clinical_variable), width = 0.2, height = 0, shape=1, size=0.7) +
    geom_text(aes(y=0.95*max_dens, label=stars),
              data=dplyr::filter(significance_annot, !is.na(stars))) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    facet_wrap(vars(cell_type), scales = "free_x", ncol=4,
               labeller = labeller(cell_type = cell_type_labels)) +
    scale_y_continuous(expand = expansion(c(0, 0.05), c(0, 0)), limits = c(0, NA)) +
    scale_x_discrete("") +
    scale_colour_brewer(palette='Dark2') +
    guides(colour = 'none') +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing.x = unit(5, 'pt'),
          axis.ticks.y = element_blank())
invisible(dev.off())
