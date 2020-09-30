# load packages
library(readr)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(beeswarm)
library(tidyr)
library(ggplot2)
library(broom)
library(tibble)
library(purrr)


# select data for plot
clin_str_density <- read_tsv('20200921_clin141_DBL_v_str_dens.tsv',
                             col_types=cols(
                                 .default = col_double(),
                                 Cascon=col_factor(levels=c(0, 1)),
                                 grade=col_factor(levels=c(1, 2, 3)),
                                 tnumber_string = col_character(),
                                 fibrosis_yn=col_factor(levels=c(0, 1)),
                                 t_number = col_character(),
                                 caco_string = col_character(),
                                 Invasief_string = col_character(),
                                 margin = col_character(),
                                 ER_status = col_character(),
                                 PR_status = col_character(),
                                 Her2status = col_character(),
                                 COX2_status = col_character(),
                                 invasive = col_character(),
                                 age_cat = col_character(),
                                 age_cat1 = col_character(),
                                 ER_1perc = col_character(),
                                 PR_1perc = col_character(),
                                 DCISdiag_cat = col_character(),
                                 diam_cat = col_character())) %>%
    mutate(
        Cascon=fct_recode(Cascon, control='0', case='1'),
        fibrosis_yn=fct_recode(fibrosis_yn, absent='0', present='1'))

# data voor plots maken
density_plot <- clin_str_density %>%
    select(t_number, Cascon, ER_status, grade, Her2status, COX2_status, fibrosis_yn,
            `CD3+_FOXP3+_NA`, `CD3+_CD8+`, `CD20+`, `CD68+`, all_t_cells, lymphocytes,
           density_Tissue_CD8._Ki67., `CD3+_CD8-`) %>%
    pivot_longer(c(`CD3+_FOXP3+_NA`, `CD3+_CD8+`, `CD20+`, `CD68+`, all_t_cells, lymphocytes,
                   density_Tissue_CD8._Ki67., `CD3+_CD8-`),
                 names_to = "cell_type",
                 values_to = "density") %>%
    pivot_longer(c(Cascon, ER_status, grade, Her2status, COX2_status, fibrosis_yn),
                 names_to = "clinical_variable",
                 values_to = "clinical_value") %>%
    filter(!is.na(clinical_value))

# New facet labels
celltype_labels <- c("lymphocytes", "All T-cells", "CD20+ B-cells", "Helper T-cells", "CD3+CD8+ T-cells", "CD3+FOXP3+ T-cells", "CD68+ cells", "CD8+Ki67+ T-cells")
names(celltype_labels) <- c("lymphocytes", "all_t_cells", "CD20+", "CD3+_CD8-", "CD3+_CD8+",  "CD3+_FOXP3+_NA", "CD68+", "density_Tissue_CD8._Ki67.")

# concatenate clinical variable with result in 1 factor variable ("combivar") and specify levels later to be used for x axis
density_plot$combivar <- paste0(density_plot$clinical_variable, '=', density_plot$clinical_value)
density_plot$combivar <- factor(density_plot$combivar, levels = c("COX2_status=High", "COX2_status=Low", "Her2status=Positive", "Her2status=Negative", "ER_status=Positive", "ER_status=Negative", "fibrosis_yn=present", "fibrosis_yn=absent", "grade=3", "grade=2", "grade=1", "Cascon=case", "Cascon=control")) 

# right order of facet labels
density_plot$cell_type <- factor(density_plot$cell_type, levels= names(celltype_labels))

# test function
non_par_test <- function(density, clinical_value) {
  clinical_value <- fct_drop(clinical_value)
  if (length(levels(clinical_value)) == 2) {
    tidy(wilcox.test(density ~ clinical_value))
  } else {
    tidy(kruskal.test(density ~ clinical_value))
  }
}

# tests
tests <- group_by(density_plot, cell_type, clinical_variable) %>% summarise(test = non_par_test(density, clinical_value))
density_plot <- left_join(density_plot, tests)
density_plot <- mutate(density_plot, significant = test$p.value < 0.05)

# make plot
pdf('boxplot_density_clinical_variables2.pdf', 7, 10)
ggplot(density_plot, aes(y=density, x=combivar)) +
    geom_boxplot(aes(colour = significant), outlier.alpha=0.0) +
    geom_jitter(width = 0.2, height = 0, shape=1, size=0.7) +
    coord_flip() +
    facet_wrap(vars(cell_type), scales = "free_x", ncol=4, labeller = labeller(cell_type = celltype_labels)) +
    scale_y_continuous(expand = expansion(c(0, 0.05), c(0, 0)), limits = c(0, NA)) +
    scale_x_discrete("", limits=clin_var_val_order) +
    scale_colour_manual(values = c("black", "red")) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.line.y =element_blank(),  
          axis.ticks.y = element_blank())
dev.off()




  








 

    
    




 
