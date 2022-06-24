if(!dir.exists(odir)) dir.create(odir)
setwd(paste0(wd, odir))

# tracking sheet -----
df_trk_2 = read_rds(glue("{wd}/df_trk_2_old.rds"))
df_trk_3 = df_trk_2 %>% dplyr::select(case_type, recurrence_type, m_if_staining_id, case_final)

# read mda file -----
df_mda = read_rds(glue("{wd}/01_read_files_mda/df_grouped_sum.rds"))
df_mda$source = "mda"
df_mda = df_mda[, -1]
df_mda$sample_id = rm_between(df_mda$sample_name, 'P5 ', '_[', extract=TRUE) %>% unlist()
df_mda$sample_id = df_mda$sample_id %>% str_replace("-", " ")
df_mda$sample_id = df_mda$sample_id %>% str_replace("-R", "")
df_mda = left_join(df_mda, df_trk_3, by = c("sample_id" = "m_if_staining_id"))
write_rds(df_mda, "df_mda.rds");write.csv(df_mda, "df_mda.csv")


# read nki file -----
df_nki = read_rds(glue("{wd}/01_read_files_nki/df_grouped_sum.rds"))
df_nki$source = "nki"
df_nki = df_nki[, -1]
unique(df_nki$sample_num)
df_nki$sample_id = rm_between(df_nki$sample_name, 'P5 ', '_[', extract=TRUE) %>% unlist()
df_nki = left_join(df_nki, df_trk_3, by = c("sample_id" = "m_if_staining_id"))
write_rds(df_nki, "df_nki.rds");write.csv(df_nki, "df_nki.csv")


# barplots ---
## only cases (Invasive recurrence) and controls (no_Inv rec or DCIS recurrence)
df_mrg = bind_rows(df_mda, df_nki)
df_mrg$recurrence_type = ifelse(df_mrg$recurrence_type == "invasive", "Invasive", df_mrg$recurrence_type)
df_mrg$recurrence = ifelse(df_mrg$recurrence_type == "Invasive", "case", "control")
write_rds(df_mrg, "df_mrg.rds");write.csv(df_mrg, "df_mrg.csv")

# barplot
df_bar = df_mrg %>% dplyr::filter(case_final == "primary")
df_bar = df_bar %>% ungroup() %>% dplyr::mutate(sample_id = factor(sample_id, levels = unique(sample_id[order(recurrence)]), ordered = TRUE))
p1 = ggplot(df_bar) +
  geom_bar(aes(x = sample_id, fill = phenotype0), position = "fill") + 
  xlab("") + ylab("") + scale_fill_npg() +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank());p1
ggsave(p1, file = "bar_primary.pdf", width = 6, height = 3)

# barplot wiht only myeloid immune cells
df_bar1 = df_mrg %>% dplyr::filter(case_final == "primary", recurrence_type %in% c("Invasive","none"), 
                                   phenotype0 %in% c("Macrophages", "Monocytes",
                                                                              "M-MDSC", "G-MDSC", "TAM2","Granulocytes"))

df_bar2 = df_bar1 %>% group_by(sample_id) %>% mutate(all_cells = n()) %>% ungroup() %>% group_by(sample_id, phenotype0) %>% mutate(tot_cells = n()) %>% ungroup()
df_meta_tum_cells = df_bar2 %>% dplyr::filter(phenotype0 == "G-MDSC") %>% group_by(sample_id) %>% mutate(tum_cells_prop = tot_cells/all_cells) %>% dplyr::select(sample_id, tum_cells_prop) %>% unique() 
df_bar2 = left_join(df_bar2, df_meta_tum_cells)
df_bar3 = df_bar2 %>% mutate(sample_id =  factor(sample_id, levels = unique(sample_id[order(recurrence, tum_cells_prop,decreasing = T)]), ordered = TRUE))

p1 = ggplot(df_bar3) +
  geom_bar(aes(x = sample_id, fill = phenotype0), position = "fill") + 
  xlab("") + ylab("") + scale_fill_npg() +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(.~recurrence);p1
ggsave(p1, file = "fig_b.pdf", width = 6, height = 3)

# barplot wiht only myeloid immune cells
df_bar1 = df_mrg %>% dplyr::filter(case_final == "primary", recurrence_type %in% c("Invasive","none"))
df_clin1 = df_bar1 %>% ungroup() %>% dplyr::select(sample_id, phenotype0, recurrence) %>% group_by(sample_id, phenotype0) %>% mutate(num_cells = n()) %>% unique() %>% ungroup() 
df_clin2 = df_clin1 %>% spread(phenotype0, num_cells) 
df_clin2[is.na(df_clin2)] <- 0
df_clin2 <- df_clin2 %>% clean_names()
df_clin2 = df_clin2 %>% mutate(total = (g_mdsc+granulocytes+m_mdsc+macrophages+m_cs+mis+monocytes+tam2+undef_mdsc+undefined_cells),
                               tumor_cells = m_cs, 
                               g_mdsc_prop = g_mdsc/total, 
                               granulocytes_prop = granulocytes/total, 
                               m_mdsc_prop = m_mdsc/total,
                               macrophages_prop = macrophages/total, 
                               monocytes_prop = monocytes/total,
                               tam2_prop = tam2/total,
                               
                               # within tumor proportion
                               g_mdsc_tum_prop = g_mdsc/tumor_cells, 
                               granulocytes_tum_prop = granulocytes/tumor_cells, 
                               m_mdsc_tum_prop = m_mdsc/tumor_cells,
                               macrophages_tum_prop = macrophages/tumor_cells, 
                               monocytes_tum_prop = monocytes/tumor_cells,
                               tam2_tum_prop = tam2/tumor_cells)
#df_clin3 = df_clin2[df_clin2 %>% complete.cases(),]
df_clin2$recurrence = factor(df_clin2$recurrence, levels = c("control", "case"))


my_comparisons <- list( c("control", "case"))
theme_sel = theme(legend.position = "none",
                  legend.background = element_blank(),
                  legend.title = element_text(size=12),
                  legend.text = element_text(size=12),
                  legend.key = element_blank(),
                  panel.grid = element_blank(),
                  panel.border = element_rect(colour = "grey", fill=NA, size=1),
                  axis.line = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_text(size = 12),
                  axis.ticks = element_blank(),
                  axis.title=element_text(size=12))
stat.test = compare_means(formula = g_mdsc_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!recurrence == 0), 
                          p.adjust.method = "fdr")
p0 = ggboxplot(df_clin2 %>% dplyr::filter(!g_mdsc_prop == 0), x = "recurrence", y = "g_mdsc_prop",title = "g_mdsc_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of g_mdsc_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p0


stat.test = compare_means(formula = granulocytes_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!granulocytes_prop == 0), 
                          p.adjust.method = "fdr")
p1 = ggboxplot(df_clin2 %>% dplyr::filter(!granulocytes_prop == 0), x = "recurrence", y = "granulocytes_prop",title = "granulocytes_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of granulocytes_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p1


stat.test = compare_means(formula = m_mdsc_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!m_mdsc_prop == 0), 
                          p.adjust.method = "fdr")
p2 = ggboxplot(df_clin2 %>% dplyr::filter(!m_mdsc_prop == 0), x = "recurrence", y = "m_mdsc_prop",title = "m_mdsc_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of m_mdsc_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p2

stat.test = compare_means(formula = macrophages_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!macrophages_prop == 0), 
                          p.adjust.method = "fdr")
p3 = ggboxplot(df_clin2 %>% dplyr::filter(!macrophages_prop == 0), x = "recurrence", y = "macrophages_prop",title = "macrophages_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of macrophages_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p3

stat.test = compare_means(formula = monocytes_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!monocytes_prop == 0), 
                          p.adjust.method = "fdr")
p4 = ggboxplot(df_clin2 %>% dplyr::filter(!monocytes_prop == 0), x = "recurrence", y = "monocytes_prop",title = "monocytes_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of monocytes_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p4

stat.test = compare_means(formula = tam2_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!tam2_prop == 0), 
                          p.adjust.method = "fdr")
p5a = ggboxplot(df_clin2 %>% dplyr::filter(!tam2_prop == 0), x = "recurrence", y = "tam2_prop",title = "tam2_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of tam2_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p5a

stat.test = compare_means(formula = tam2_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!tam2_tum_prop == 0), 
                          p.adjust.method = "fdr")
p5 = ggboxplot(df_clin2 %>% dplyr::filter(!tam2_tum_prop == 0), x = "recurrence", y = "tam2_tum_prop",title = "tam2_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of tam2_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p5

stat.test = compare_means(formula = g_mdsc_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!g_mdsc_tum_prop == 0), 
                          p.adjust.method = "fdr")
p6 = ggboxplot(df_clin2 %>% dplyr::filter(!g_mdsc_tum_prop == 0), x = "recurrence", y = "g_mdsc_tum_prop",title = "g_mdsc_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of g_mdsc_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p6

stat.test = compare_means(formula = granulocytes_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!granulocytes_tum_prop == 0), 
                          p.adjust.method = "fdr")
p7 = ggboxplot(df_clin2 %>% dplyr::filter(!granulocytes_tum_prop == 0), x = "recurrence", y = "granulocytes_tum_prop",title = "granulocytes_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of granulocytes_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p7

stat.test = compare_means(formula = m_mdsc_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!m_mdsc_tum_prop == 0), 
                          p.adjust.method = "fdr")
p8 = ggboxplot(df_clin2 %>% dplyr::filter(!m_mdsc_tum_prop == 0), x = "recurrence", y = "m_mdsc_tum_prop",title = "m_mdsc_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of m_mdsc_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p8

stat.test = compare_means(formula = macrophages_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!macrophages_tum_prop == 0), 
                          p.adjust.method = "fdr")
p9 = ggboxplot(df_clin2 %>% dplyr::filter(!macrophages_tum_prop == 0), x = "recurrence", y = "macrophages_tum_prop",title = "macrophages_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of macrophages_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p9
stat.test = compare_means(formula = monocytes_tum_prop ~ recurrence, data = df_clin2 %>% dplyr::filter(!monocytes_tum_prop == 0), 
                          p.adjust.method = "fdr")
p10 = ggboxplot(df_clin2 %>% dplyr::filter(!monocytes_tum_prop == 0), x = "recurrence", y = "monocytes_tum_prop",title = "monocytes_tum_prop",
               fill = "recurrence", palette = c("orange", "purple"),
               add = "jitter") + ylab("(Fraction of monocytes_tum_prop)") +
  theme_sel +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, y.position = 0.015, color = "black", linetype = 1);p10

library(cowplot)
p = plot_grid(p0,p1, p2, p3, p4, p5a, nrow = 1)
pdf(file = paste0("boxplots_allcell_prop.pdf"), width = 11, height = 3.5)
print(p)
dev.off()
p_tum = plot_grid(p6, p7, p8, p9, p10, p5, nrow = 1)
pdf(file = paste0("boxplots_tumcell_prop.pdf"), width = 11, height = 3.5)
print(p_tum)
dev.off()
