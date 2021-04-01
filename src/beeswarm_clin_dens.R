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

# Careful, the order should match that of celltype_labels
density_vars <- c("lymphocytes", "all_t_cells", "CD20+", "CD3+_CD8-",
                  "CD3+_CD8+", "CD3+_FOXP3+_NA", "CD68+", "density_Tissue_CD8._Ki67.")
clin_boxplot_vars <- c('Cascon', 'ER_status', 'grade', 'Her2status', 'COX2_status', 'fibrosis_yn')
clin_vars_cat <- c(clin_boxplot_vars,
    'age_cat1', 'PR_status', 'clin_pres', 'Yearsto_iIBC_cat_cases', 'diameter_cat', 'margin',
    'dom_growthpat', 'necrosis', 'calcs')
clin_vars_cont <- character()

clin_vars <- c(clin_vars_cat, clin_vars_cont)

# Labels for cell types in the plot
celltype_labels <- c("lymphocytes", "All T-cells", "CD20+ B-cells", "Helper T-cells",
                     "CD3+CD8+ T-cells", "CD3+FOXP3+ T-cells", "CD68+ cells", "CD8+Ki67+ T-cells")
names(celltype_labels) <- density_vars

#############################
# Read and prepare the data #
#############################

col_spec <- c(
    structure(map(density_vars, ~ col_double()), names=density_vars),
    structure(map(clin_vars_cat, ~ col_character()), names=clin_vars_cat),
    structure(map(clin_vars_cont, ~ col_double()), names=clin_vars_cont))

col_spec[['t_number']] <- col_character()
col_spec[['Cascon']] <- col_factor(levels=c('0', '1'))
col_spec[['fibrosis_yn']] <- col_factor(levels=c('0', '1'))
col_spec[['grade']] <- col_factor(levels=c('1', '2', '3'))
clin_density <- read_tsv(
        'data/clin_DBL_v_str_dens.tsv',
        col_types = do.call(cols_only, col_spec)) %>%
    mutate(
        Cascon=fct_recode(Cascon, control='0', case='1'),
        fibrosis_yn=fct_recode(fibrosis_yn, absent='0', present='1'))

# Move the density and clinical variables into a single column.
clin_density <- clin_density %>%
    pivot_longer(all_of(density_vars),
                 names_to = "cell_type",
                 values_to = "density")

clin_density_cat <- clin_density %>%
    select(-all_of(clin_vars_cont)) %>%
    pivot_longer(all_of(clin_vars_cat),
                 names_to = "clinical_variable",
                 values_to = "clinical_value") %>%
    dplyr::filter(!is.na(clinical_value)) %>%
    mutate(cell_type = factor(cell_type, levels=density_vars))


clin_density_cont <- if (length(clin_vars_cont) > 1) {
    clin_density %>%
        select(-all_of(clin_vars_cat)) %>%
        pivot_longer(all_of(clin_vars_cont),
                     names_to = "clinical_variable",
                     values_to = "clinical_value") %>%
        dplyr::filter(!is.na(clinical_value)) %>%
        mutate(cell_type = factor(cell_type, levels=density_vars))
} else {
    NULL
}


#########################
# Test for significance #
#########################

# Perform a non-parametric test, automaticly choosing wilcox or kruskal based on the number of the
# number of levels in the variable.
non_par_test <- function(density, clinical_value) {
  clinical_value <- fct_drop(clinical_value)
  if (length(levels(clinical_value)) == 2) {
    tidy(coin::wilcox_test(density ~ clinical_value))
  } else {
    tidy(coin::kruskal_test(density ~ clinical_value))
  }
}

tests_cat <- clin_density_cat %>%
    group_by(cell_type, clinical_variable) %>%
    summarise(test = non_par_test(density, clinical_value)) %>%
    unpack(test) %>%
    rename(nominal_p = p.value)

tests_cont <- if (length(clin_vars_cont) > 1) {
    clin_density_cont %>%
        group_by(cell_type, clinical_variable) %>%
        summarise(test = tidy(cor.test(density, clinical_value, method = 'spearman',
                                       exact = F))) %>%
        unpack(test) %>%
        rename(nominal_p = p.value)
} else {
    NULL
}

tests <- bind_rows(tests_cat, tests_cont) %>%
    mutate(p = p.adjust(nominal_p, 'bonferroni'),
           fdr =  p.adjust(nominal_p, 'BH')) %>%
    arrange(cell_type, clinical_variable)
write_xlsx(tests, 'results/significance_stroma.xlsx')

############
# Plotting #
############


# Prepare Data #

# clin_var_val is a variable for the x axis and needs to be in a specific order.
clin_density_plot <- clin_density_cat %>%
    dplyr::filter(clinical_variable %in% clin_boxplot_vars) %>%
    mutate(clin_var_val = paste0(clinical_variable, '=', clinical_value) %>%
        # FIXME: generate this, but it needs to be in order, and input validation
        factor(levels = c("COX2_status=High", "COX2_status=Low",
                          "Her2status=Positive", "Her2status=Negative",
                          "ER_status=Positive", "ER_status=Negative",
                          "fibrosis_yn=present", "fibrosis_yn=absent",
                          "grade=3", "grade=2", "grade=1",
                          "Cascon=case", "Cascon=control")))

# We need to to know maximum density for significance star placement
dens_summary <- clin_density_cat %>%
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

pdf('plots/boxplot_density_clinical_variables.pdf', 7, 10)
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
               labeller = labeller(cell_type = celltype_labels)) +
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

##############################
# Compute summary statistics #
##############################

summary_stats <- clin_density_cat %>%
  group_by(cell_type, clinical_variable, clinical_value) %>%
  summarize(
    median = median(density, na.rm = T),
    p25 = quantile(density, probs = .25, na.rm = T),
    p75 = quantile(density, probs = .75, na.rm = T)
  )

write_xlsx(summary_stats, 'results/summary_stats.xlsx')
