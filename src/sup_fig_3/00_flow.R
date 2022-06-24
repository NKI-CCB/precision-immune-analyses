# setup
library(pacman)
p_load(tidyverse, janitor, glue)
p_load(cowplot)
#install.packages("spatstat")
p_load(spatstat)
p_load(qdapRegex)
p_load(stringr)
p_load(ggpubr)
library(ggsci)
options(scipen=999)
p_load(ggplot2)
p_load(params)
p_load(openxlsx)


if(FALSE){
        wd = "~/rsrch2_tseth/projects/vectra/analysis/dcis_v3/"
        dat_path = "~/rsrch2_tseth/projects/vectra/data/dcis_panel/"
}

# wd
wd = "~/sea_rsrch3_home/rsrch2_backup/projects/vectra/analysis/myeloid_paper_v1/"
#wd = "/home/rstudio10/.rsrch2/genomic_med/tseth/projects/vectra/analysis/dcis_v3/"
dat_path1 = "~/sea_rsrch3_home/rsrch2_backup/projects/vectra/data/dcis_panel_sept19/"
dat_path2 = "~/sea_rsrch3_home/rsrch2_backup/projects/vectra/data/dcis_panel/panel5/"

setwd(wd)
#setwd("/home/rstudio10/.rsrch2/genomic_med/tseth/projects/vectra/analysis/dcis_v3/")
# read data -----
odir = "data"

#' if we change something, will save in the same dir with new versions
odir = "01_read_files_mda"
source("01_read_files_mda.R")
odir = "01_read_files_nki/"
source("01_read_files_nki.R")

odir = "01a_plot_freq"
source("01a_plot_freq.R")

odir = "01_read_tissue_mda"
source("01_read_tissue_mda.R")
odir = "01_read_tissue_nki"
source("01_read_tissue_nki.R")# runs smoothly

# read meta (only read final df_trk_3) ------
trk_fl = glue("~/sea_rsrch3_home/rsrch2_backup/projects/vectra/data/dcis_panel/DCIS vectra progress 12-12-2017.xlsx")
df_trk = params::read_sheet(trk_fl, id_column = "MRN") %>% 
        clean_names()
# added to only include new samples
df_trk = df_trk %>% dplyr::filter(filt_accession == "99")

df_trk$ihc_sample_id = df_trk$m_if_staining_id %>% 
        #tolower() %>% 
        #str_replace("-", "_") %>% 
        str_replace_all(" ", "_")
df_trk$case_type = tolower(df_trk$case_type)
df_trk = df_trk %>% mutate(recurrence_type = case_when(recurrence_type == "Invasive" ~ "invasive", 
                                                       recurrence_type == "invasive" ~ "invasive",
                                                       recurrence_type == "DCIS" ~ "DCIS",
                                                       recurrence_type == "none" ~ "none"))
        #gsub(" ", "_", .)
df_trk = group_by(df_trk, mrn) %>% # i changed this to mrn
        mutate(rec_type = recurrence_type) # since for nki, there is only 1 sample/pt, take just 1 case_type row
df_trk$lbl = glue_data(df_trk, "pt{orid_pt_id}_{recurrence_type}");df_trk$lbl
rownames(df_trk) = df_trk$ihc_sample_id
write_rds(df_trk, paste0(wd, "df_trk.rds") )

df_trk_2 = df_trk %>% mutate(case_final = if_else(case_type == "dcis", "primary", "recurrence"))
write_rds(df_trk_2, paste0(wd, "df_trk_2.rds") )
df_trk_2 = read_rds("~/sea_rsrch3_home/rsrch2_backup/projects/vectra/analysis/dcis_v3_sept19_v2/df_trk_2.rds")
df_trk_3 = df_trk_2 %>% dplyr::select(patient_id, orid_pt_id, accession, case_type, recurrence_type, m_if_staining_id, institute,
                                      ihc_sample_id, rec_type, lbl, case_final)
write_rds(df_trk_3,paste0(wd, "df_trk_3.rds") )
write.csv(df_trk_3, paste0(wd, "final_tracker.csv") )

df_trk_3 = read_rds("~/sea_rsrch3_home/rsrch2_backup/projects/vectra/analysis/dcis_v3_sept19_v2/df_trk_3.rds")


# tracking sheet for old panel ---------
trk_fl = glue("~/sea_rsrch3_home/rsrch2_backup/projects/vectra/data/dcis_panel/DCIS vectra progress 12-12-2017.xlsx")
df_trk = params::read_sheet(trk_fl, id_column = "MRN") %>% 
        clean_names()
df_trk$ihc_sample_id = df_trk$m_if_staining_id %>% 
        tolower() %>% str_replace("-", "_") %>% 
        str_replace_all(" ", "_")
df_trk$case_type = tolower(df_trk$case_type) %>% 
        gsub(" ", "_", .)
df_trk = group_by(df_trk, patient_id) %>% 
        mutate(rec_type = sort(case_type)[2])
df_trk$lbl = glue_data(df_trk, "pt{patient_id}_{case_type}");df_trk$lbl
rownames(df_trk) = df_trk$ihc_sample_id
write_rds(df_trk, paste0(wd, "df_trk_old.rds") )

df_trk_2 = df_trk %>% mutate(case_final = if_else(case_type == "dcis", "primary", "recurrence"))
write_rds(df_trk_2, paste0(wd, "df_trk_2_old.rds") )

# simple plots ------
odir = "02_plots"
source('02_plots.R')

odir = "02_plots_mda"
source('02_plots_mda.R')

