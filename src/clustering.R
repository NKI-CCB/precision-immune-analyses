library(RColorBrewer)
library(gplots)
library(ComplexHeatmap)

clinical <- read.delim("data/clin_DBL_v_str_dens.tsv", stringsAsFactors = F, check.names = F)
clinical$ki67_cat = factor(clinical$ki67perc_t >= 14, levels = c(F, T), labels = c('<14', '>=14'))
hmcolors <- colorRampPalette(c("green","black","red"))(17)
hmcolors <- brewer.pal(9, "Purples")

labels <- c("Stroma", "DCIS")

for (label in labels) {
  dfSummary <- read.delim(paste("data/cell_density_", label, ".tsv", sep=""), stringsAsFactors = F, check.names = F)
  
  dfKi <- read.delim("data/ki67_density.tsv", stringsAsFactors = F, check.names = F)
  ki_mapping <- match(dfSummary$t_number, dfKi$t_number)
  if (label == "DCIS") {
    dfSummary$`CD3+_KI67+` <- dfKi$`density_DCIS_CD8+_Ki67+`[ki_mapping]
  } else {
    dfSummary$`CD3+_KI67+` <- dfKi$`density_Tissue_CD8+_Ki67+`[ki_mapping]
  }
  
  data <- log(1+data.matrix(dfSummary[,-c(1,2)]))
  data <- data[,-which(colnames(data)=='Other')]
  data <- data[,-which(colnames(data)=='panCK+')]
  
  colnames(data)[which(colnames(data) == "CD3+_CD8-")] <- "Helper T-cells"
  colnames(data)[which(colnames(data) == "CD3+_CD8+")] <- "CD8+ T-cells"
  colnames(data)[which(colnames(data) == "CD20+")] <- "CD20+ B-cells"
  colnames(data)[which(colnames(data) == "CD68+")] <- "CD68+ cells"
  colnames(data)[which(colnames(data) == "CD3+_FOXP3+")] <- "FOXP3+ T-cells"
  colnames(data)[which(colnames(data) == "CD3+_KI67+")] <- "CD8+Ki67+ T-cells"
 
  mapping <- match(dfSummary$t_number, clinical$t_number)
  
  pdf(paste("plots/heatmap_density_", label, ".pdf", sep=""), width=15, height=5)
  
  set_colors <- brewer.pal(9, "Set1")
  paired_colors <- brewer.pal(12,"Paired")
  gradient_colors <- brewer.pal(3, "Reds")
  batch_colors <- brewer.pal(7,"Set3")
  
  dfAnnotation <- data.frame(
    case_control=clinical$Cascon[mapping],
    ER_status=clinical$ER_status[mapping],
    HER2_status=clinical$Her2status[mapping],
    COX2_status=clinical$COX2_status[mapping],
    grade=clinical$grade[mapping],
    fibrosis=clinical$fibrosis_yn[mapping],
    Ki67_status=clinical$ki67_cat[mapping]
  )
  ha = HeatmapAnnotation(
    df = dfAnnotation,
    col = list(grade = c('1' = gradient_colors[1], '2' = gradient_colors[2], '3' = gradient_colors[3]),
               ER_status = c('Negative' = paired_colors[1],'Positive' = paired_colors[2]),
               HER2_status = c('Negative' = paired_colors[5],'Positive' = paired_colors[6]),
               COX2_status = c('Low' = paired_colors[3],'High' = paired_colors[4]),
               case_control = c('0' = paired_colors[7],'1' = paired_colors[8]),
               fibrosis = c('0' = paired_colors[11],'1' = paired_colors[12]),
               Ki67_status = c('<14' = paired_colors[9],'>=14' = paired_colors[10])),
    na_col = 'white'
  )
  ht1 = Heatmap(t(data), name = "log(1+cells/mm^2)",
                clustering_distance_rows = "euclidean",
                column_title = label, top_annotation = ha,
                col=hmcolors, column_names_gp=gpar(cex=0.5))
  draw(ht1)
  res <- dev.off()
}
