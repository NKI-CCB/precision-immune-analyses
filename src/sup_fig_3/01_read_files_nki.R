library(qdapRegex)
#install.packages("stringi")
library(tidyverse)
library(janitor)

if(!dir.exists(odir)) dir.create(odir)
setwd(paste0(wd, odir))

wd = getwd()

# list files ----------
# "/home/rstudio10/rsrch2_home/projects/vectra/data/dcis_panel_sept19//P5 DPRE Arg1/P5 DPRE 01 Arg1/P5 DPRE 01_[43155,17659]_cell_seg_data.txt"
ls.files = list.files(path = paste0(dat_path1), pattern = "*_cell_seg_data.txt", full.names = TRUE, 
                      recursive = TRUE)
write_sheet(x = ls.files, file = paste0(wd, "ls.files.xlsx"))
head(ls.files);length(ls.files)

# "P5 DPRE Arg1/P5 DPRE 01 Arg1/P5 DPRE 01_[43155,17659]_cell_seg_data.txt"
ls.files.names = list.files(path = paste0(dat_path1), pattern = "*_cell_seg_data.txt", full.names = FALSE, 
                            recursive = TRUE, include.dirs = FALSE)
head(ls.files.names);length(ls.files.names)

# Remove restain string "P5 DPRE 15 RESTAIN_[47559,21258]_cell_seg_data.txt" 
# "P5 DPRE 14_[41481,9046]_cell_seg_data.txt"
# ls.files = gsub(' RESTAIN','', ls.files)

# parse filename i.e # [1] "P5 DPRE 01"
samplenums_0 = rm_between(ls.files, '/P5', '_[', extract=TRUE)
samplenums = str_sub(samplenums_0,-10,-1)
unique(samplenums) # 22 unique values

# ihc_sample_id: id used by the vectra facility
# <panel> <disease sample id> <location>
#ihc_sampleid = rm_between(ls.files, 'P5 ', '_[', extract=TRUE) %>% unlist()
#ihc_sampleid = gsub("-R$", "", ihc_sampleid)
# move to lower case, and NO spaces etc
#ihc_sampleid = tolower(ihc_sampleid) %>% str_replace("-", "_") %>%
#       str_replace(" ", "_")

# get samples num i.e # [1] "01" and sample name,ihc_sampleid DPRE 01, ihc_sampleid_orig "P5 DPRE 18 CK CD68" and section cooridnates "47559,21258"
num = str_sub(samplenums,9,10)
ihc_sampleid = str_sub(samplenums,4,10)
ihc_sampleid_orig = rm_between(ls.files, '/', '/', extract=TRUE) %>% lapply(., `[[`, 10) %>% as.character()
sec = rm_between(ls.files, '_[', ']_', extract=TRUE) %>% unlist()
marker = rm_between(samplenums_0, 'DPRE ', '/P5', extract=TRUE) %>% lapply(., `[[`, 1) %>% as.character()
unique(marker) #"Arg1"    "CD11B"   "CD14"    "CD33"    "CD66B"   "CK CD68"

# put as a df
meta_df = cbind(ls.files, num, ihc_sampleid, ihc_sampleid_orig, sec, marker) %>% as_tibble()
# meta_df$rec_type = if_else(meta_df$rec_type == "R", "R","NR")
colnames(meta_df) = c("file_path","sample_num", "ihc_sampleid", "ihc_sampleid_orig", "sec", "marker")
meta_df = dplyr::group_by(meta_df, sample_num) %>% dplyr::mutate(section = length(unique(sec)), count = n(), 
                                                                 sc = rep(1:section, count%/%section, len = count)) %>%
  ungroup() %>%
  mutate(panel = 5)
#meta_df$ihc_sampleid = ihc_sampleid

# Save sample tracker
#write.csv(meta_df, paste0(wd, "sample_tracker.csv"), quote = FALSE)
write_sheet(meta_df, paste0(wd, "sample_tracker.xlsx"))
write_rds(meta_df, paste0(wd, "sample_tracker.rds"))

# read all files ------
lst = list()
# force all to be character!
# so that they are read quickly and correctly
# can make a few of them numeric later
dim(meta_df) # 690 10
concat.mat = lapply(1:nrow(meta_df), function(i){
  message("working on ",i)
  df = read_tsv(meta_df$file_path[i], col_types = cols(.default = col_character())) %>%
    clean_names()
  df$sample_num = meta_df$sample_num[i]
  #df$repeated = meta_df$repeated[i]
  df$panel = "5"
  df$section = meta_df$sc[i]
  df$ihc_sample_id = meta_df$ihc_sampleid[i]
  df$ihc_sample_id_orig = meta_df$ihc_sampleid_orig[i]
  df
  #lst[[i]] = df no need, not using lst
  #do.call(rbind, lst)
})

df = bind_rows(concat.mat)
# get rid fo restain from P5 DPRE 22 RESTAIN_[54204,15216].im3 to P5 DPRE 22_[54204,15216].im3
df$path = gsub(' RESTAIN','', df$path)
df$sample_name = gsub(' RESTAIN','', df$sample_name)

dim(df) #[1] 1146834      27
write_sheet(x = df, file = paste0(wd, "/intial_df.xlsx"))
write_rds(x = df, path = paste0(wd, "/intial_df.rds"))

df_sel = df %>% dplyr::select(-nucleus_area_square_microns, -cytoplasm_area_square_microns, -membrane_area_square_microns, -entire_cell_area_square_microns,
                              -cytoplasm_cd33_opal_520_mean_normalized_counts_total_weighting, -nucleus_cd33_opal_520_mean_normalized_counts_total_weighting,
                              -entire_cell_cd33_opal_520_mean_normalized_counts_total_weighting)
df_sel$source = "nki"

write_sheet(x = df_sel, file = paste0(wd, "/intial_df_sel.xlsx"))
write_rds(x = df_sel, path = paste0(wd, "/intial_df_sel.rds"))

# df_sel$phenotype = coalesce(df_sel$pheno2, "NA")
# df_sel$phenotype0 = rename_marks(df_grouped_sum$pheno2)


unique(df_sel$sample_num)
for (i in 1:length(unique(df$sample_num))) {
  message("working on-", unique(df_sel$sample_num)[i])
  df_sel_samp = df_sel %>% dplyr::filter(sample_num == unique(df_sel$sample_num)[i])
  write_sheet(x = df_sel_samp, file = paste0(wd, "/df_sel_samp_", unique(df_sel$sample_num)[i], ".xlsx"))
}

# STOP HERE -------------
# 
# library(caret)  
# library(Rtsne)
# mat_phens_wd = df_mat_phens_wd
# mat_phens_wd[is.na(mat_phens_wd)] <- 0
# mat_phens_wd = data.matrix(mat_phens_wd)
# tsne_model_1 = Rtsne(as.matrix(mat_phens_wd), check_duplicates=FALSE, pca=TRUE, perplexity=50, theta=0.5, dims=2)
# d_tsne_1 = as.data.frame(tsne_model_1$Y)
# ggplot(d_tsne_1, aes(x=V1, y=V2)) +  
#   geom_point(size=0.25) +
#   guides(colour=guide_legend(override.aes=list(size=6))) +
#   xlab("") + ylab("") +
#   ggtitle("t-SNE") +
#   theme_light(base_size=20) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_colour_brewer(palette = "Set2")
# fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))
# ## setting 3 clusters as output
# d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3))  

# group by cell id, position and check what all phenotypes are expressed - results as single phenotypes, i.e each cell-each phenotype
df_grouped = df_sel %>% filter(!phenotype %in% c("Others"), !is.na(phenotype)) %>% 
  select(path, tissue_category,sample_name, phenotype, cell_id, cell_x_position, cell_y_position) %>% 
  group_by(cell_id, cell_x_position, cell_y_position) %>% 
  mutate(pheno2 = paste0(phenotype)) %>% ungroup();dim(df_grouped) #89062     8

# Group cells and summarise ith unique phenotype - 53 combination of phenotypes
df_grouped_sum = df_sel %>% 
  select(path, tissue_category,sample_name, sample_num, phenotype, cell_id, cell_x_position, cell_y_position) %>% 
  group_by(sample_num, cell_id, cell_x_position, cell_y_position) %>% 
  mutate(pheno2 = paste0(unique(phenotype), collapse = ",")) %>% unique();length(unique(df_grouped_sum$pheno2));dim(df_grouped_sum) #1131327       9

# count the number of cells in each sample expressing set of 1/53 phenotypes
df_grouped_norm1 = df_grouped_sum %>%  
  group_by(pheno2, sample_name) %>% 
  summarise(num_cells = n());dim(df_grouped_norm1) #1472    3

# count the number of cells xpressing set of 1/69 phenotypes
df_grouped_norm2 = df_grouped_sum %>%  
  group_by(pheno2) %>% 
  summarise(num_cells = n());dim(df_grouped_norm2) #53  2

write.csv(df_grouped_norm1, paste0(wd, "cell_pheno_per_sample.csv"), quote = FALSE)
write_rds(df_grouped_norm1, paste0(wd, "cell_pheno_per_sample.rds"))
write.csv(df_grouped_norm2, paste0(wd, "cell_pheno.csv"), quote = FALSE)
write_rds(df_grouped_norm2, paste0(wd, "cell_pheno.rds"))
write.csv(unique(df_grouped_sum$pheno2), paste0(wd, "53marker_combinations.csv"), quote = FALSE)

rename_marks <- function(x){
  y = fct_collapse(x, 
                   mis = c("Others", "NA"),
                   "MCs" = c("Others,CK+"),
                   
                   "Macrophages" = c("Others,CD11B+,CD33+,CD68+", "Others,CD11B+,CD68+",  "Others,CD33+,CD68+", 
                                     "Others,CD11B+,CD14+,CD68+", "Others,CD11B+,CD14+,CD33+,CD68+"), # controversial macro, macro with CD14+
                   
                   "Monocytes" = c("Others,CD11B+,CD14+", "Others,CD14+,CD33+"),
                   
                   "TAM2" = c("Arg1+,CD11B+,Others,CD68+", "Arg1+,CD11B+,CD14+,CD33+,Others,CD68+", "Arg1+,CD11B+,CD14+,Others,CD68+", 
                              "Arg1+,Others,CD33+,CD68+"),
                   
                   "Granulocytes" = c("Others,CD11B+,CD33+,CD66B+", "Others,CD11B+,CD66B+", "Others,CD33+,CD66B+"),
                   
                   "M-MDSC" = c("Others,CD11B+,CD14+,CD33+"),		
                   "G-MDSC" = c("Arg1+,CD11B+,Others,CD33+,CD66B+", "Arg1+,CD11B+,Others,CD66B+", "Arg1+,Others,CD33+,CD66B+"),
                   "Undef-MDSC" = c("Arg1+,CD11B+,Others", "Arg1+,CD11B+,Others,CD33+"))
  y
}

# rename_marks_v2 <- function(x){
#   y = fct_collapse(x,
#                    mis = c("Others", "NA"),
#                    "MCs" = c("Others,MC"),
#                    "Macrophages" = c("Others,CD11b+,CD68"),
#                    "Total Macrophages" =c("Others,CD11b+,CD33,CD68"),
#                    "TAM2" = c("Arg1,CD11b+,Others,CD68"),
#                    "Total Granulocytes" = c("Others,CD11b+,CD66b"),
#                    "MDSC_monocytic" = c("Others,CD11b+,CD14,CD33", "Others,CD11b,CD14,CD33"),
#                    "MDSC_granulocytic" = c("Others,CD11b+,CD33", "Others,CD11b,CD33"),
#                    "Monocytes" = c("Others,CD11b+,CD14", "Others,CD11b,CD14", "Others,CD11b+,CD14,Other"))
#   y
# }
df_grouped_sum$phenotype = coalesce(df_grouped_sum$pheno2, "NA")
df_grouped_sum$phenotype0 = rename_marks(df_grouped_sum$pheno2)
#df_grouped_sum$phenotype_tmp = rename_marks_v2(df_grouped_sum$pheno2)

# initial df
# df_sel$phenotype = coalesce(df_sel$phenotype, "NA")
# df_sel$phenotype0 = rename_marks(df_sel$pheno2)


# additional step
pheno_grps = c("mis", "MCs", "Macrophages", "Monocytes", "TAM2", "Granulocytes", "M-MDSC", "G-MDSC", "Undef-MDSC")
df_grouped_sum =  df_grouped_sum %>% mutate(phenotype0 = if_else(condition = phenotype0 %in% pheno_grps, 
                                                                 true = as.character(phenotype0),
                                                                 false =  "undefined cells", missing = "missing"))

write_csv(df, glue("{wd}df_merged.csv.gz"))
write_rds(df, glue("{wd}df_merged.rds"))

write_csv(df_grouped_sum, glue("{wd}df_grouped_sum.csv.gz"))
write_csv(df_grouped_sum, glue("{wd}df_grouped_sum.csv"))
write_rds(df_grouped_sum, glue("{wd}df_grouped_sum.rds"))

# checks
unique(df_grouped_sum$phenotype0)
table(df_grouped_sum$phenotype0)
table(df_grouped_sum$phenotype) %>% View()
