library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(reshape2)
library(gplots)
require(gtools)
library(psych)
library(RColorBrewer)
library(mclust)
library(factoextra)
library(ggrepel)
library(GGally)
library(Hmisc)
library(corrplot)
library(network)
library(sna)
library(intergraph)
library(destiny)
library(ggbiplot)
library(monocle3)
library(viridis)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(igraph)
library(intergraph)
library(ggVennDiagram)
library(stringr)
library(topGO)
library(WGCNA)
library("openxlsx")
library(lsmeans)
library(emmeans)
library(multcomp)
library(ggcorrplot)
library("ggpubr")
library(scran)
library(ggplot2)                  
library(ggalluvial)
library(igraph)
library(limma)
library(xgboost)
library(caTools)
library(gprofiler2)
library(DoubletFinder)
library(harmony)

#Fig.3a
EPI<-readRDS("Epithelial.rds")
data<-EPI
data$orig.ident <- factor(data$orig.ident,levels = c("EcO1","EcO2","EcO3","EcO4","EcO5",
                                                     "EcOA1","EcOA2","EcOA3","EcOA4","EcOA5",
                                                     "Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5",
                                                     "CA1","CA2","CA3"))
data$annotation_new_final <- data$cellclass
data$annotation_new_final <- as.character(data$annotation_new_final)
Epi_in_T <- data$diseasestate == "CA" & data$annotation_new_final == "Epithelial"
data$annotation_new_final[Epi_in_T] <- "Malignant"

Tumor_Epithelial <- subset(data, annotation_new_final %in% c("Malignant","Epithelial"))
data <- Tumor_Epithelial
normal_data <- subset(data, tissue %in% c("Ctrl"))
normal_data$orig.ident <- "Normal"

sample <- as.data.frame(table(data$orig.ident))
sample_filtered <- filter(sample, Freq >= 2)
data_filtered <- subset(data, orig.ident %in% sample_filtered$Var1)

diff <- NULL
for (i in unique(data_filtered$orig.ident)) {
  tmp <- subset(data, orig.ident == i)
  tmp_data <- merge(normal_data, tmp,add.cell.ids = c("Normal","i"))
  tmp_data <- JoinLayers(tmp_data)
  Idents(tmp_data) <- tmp_data$orig.ident
  differential_test <- FindMarkers(tmp_data, ident.1 = i, ident.2 = "Normal", verbose = FALSE, min.pct = 0.15,
                                   logfc.threshold = 0.25, min.cells.feature =0, min.cells.group = 0,max.cells.per.ident = 200,
                                   test.use = "MAST")
  colnames(differential_test) <- paste0(colnames(differential_test), "_", i)
  print(paste0("Producing diff genes for ", i))
  if (is.null(diff)) {
    diff <- differential_test
  } else {
    diff <- merge(diff, differential_test, by=0, all=TRUE)
    rownames(diff) <- diff$Row.names
    diff <- diff[, !colnames(diff) %in% c("Row.names")]
  }
}
write.csv(diff, "Tumor_normal_tissue_findmarkers_match.csv")

diff <- read.csv("Tumor_normal_tissue_findmarkers_match.csv", row.names = 1)
diff_avg_logFC <- diff[, grepl("avg_log2FC", colnames(diff))]
diff_p_val_adj <- diff[, grepl("p_val_adj", colnames(diff))]
num = 1
diff_avg_logFC_significant <- diff_avg_logFC[(rowSums(diff_avg_logFC > num & diff_p_val_adj < 0.05) > 1 | rowSums(diff_avg_logFC < (-num) & diff_p_val_adj < 0.05) > 1), ]
diff_p_val_adj_significant <- diff_p_val_adj[(rowSums(diff_avg_logFC > num & diff_p_val_adj < 0.05) > 1 | rowSums(diff_avg_logFC < (-num) & diff_p_val_adj < 0.05) > 1), ]
clean_data <- diff_avg_logFC_significant[complete.cases(diff_avg_logFC_significant), ]
partial_data <- diff_avg_logFC_significant[!apply(is.na(diff_avg_logFC_significant), 1, all), ]
pcs <- prcomp(clean_data)
pc_df <- data.frame(pcs$rotation)
pc_df$sample <- substr(rownames(pc_df), 12, nchar(rownames(pc_df)))
table(pc_df$sample)
pc_df$DiseaseState <- ifelse(
  grepl("CA", pc_df$sample), "CA",  
  ifelse(
    grepl("EcOA", pc_df$sample), "EcOA",
    ifelse(
      grepl("EcO", pc_df$sample), "EcO",
        "Normal"
      )
    )
  )
table(pc_df$DiseaseState)
pc_df <- pc_df[, c("PC1", "PC2", "DiseaseState")]
scalefactor2 <- -1
pc_df$PC2 <- pc_df$PC2 * scalefactor2
pc_df$DiseaseState <- factor(pc_df$DiseaseState, levels = c("Normal","EcOA", "EcO", "CA"))
ggplot(pc_df, aes(x=PC1, y=PC2, color=DiseaseState)) +
  geom_point(size=2) + scale_color_manual(values=c("#e8c559", "darkgreen","darkblue", "#89288F", "darkred"))

fit <- lm(PC2 ~ bs(PC1,df=3), data = pc_df)  # df=4提供合理的灵活性

age.grid <- seq(from = min(pc_df$PC1), to = max(pc_df$PC1), length.out = 10000)

splinefit <- data.frame(
  PC1 = age.grid,
  PC2 = predict(fit, newdata = data.frame(PC1 = age.grid))
)

pc1_var <- 100 * summary(pcs)$importance["Proportion of Variance", "PC1"]
pc2_var <- 100 * summary(pcs)$importance["Proportion of Variance", "PC2"]

minDisToSpline <- vector(mode = "numeric", length = nrow(pc_df))
minXs <- vector(mode = "numeric", length = nrow(pc_df))
minYs <- vector(mode = "numeric", length = nrow(pc_df))

for (i in 1:nrow(pc_df)) {
  x <- pc_df[i, "PC1"]
  y <- pc_df[i, "PC2"]

  distances <- sqrt((splinefit$PC1 - x)^2 + (splinefit$PC2 - y)^2)

  min_idx <- which.min(distances)
  minDisToSpline[i] <- distances[min_idx]
  minXs[i] <- splinefit$PC1[min_idx]
  minYs[i] <- splinefit$PC2[min_idx]
}

pc_df$NearestPointX <- minXs
pc_df$NearestPointY <- minYs

min_x <- min(splinefit$PC1)
max_x <- max(splinefit$PC1)
pc_df$Normalized_Rank <- ((pc_df$NearestPointX - min_x) / (max_x - min_x)) * 100

p<-ggplot(pc_df, aes(x = PC1, y = PC2)) +
  geom_line(data = splinefit, 
            aes(x = PC1, y = PC2), 
            color = "black", 
            linewidth = 1.2,
            linetype = "solid") +
  geom_point(aes(color = DiseaseState), size = 3, alpha = 0.8) +
  scale_color_manual(values=c("#8dd3c7","#80b1d3","#bebada","#fb8072"), name = "")+
  xlab(paste0("PC1 (", round(pc1_var, 1), "%)")) + 
  ylab(paste0("PC2 (", round(pc2_var, 1), "%)")) +
  ggtitle("Epithelial Cells") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    aspect.ratio = 1,
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  ) +
  coord_cartesian(clip = 'off') 
ggsave("curve.png",p,width = 6,height = 5)

#Fig.3b,d,Fig.S3a
library(msigdbr)
go <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)
markers <- read.csv("markers_cnv.csv")
markers <- subset(markers, p_val<0.05 & avg_log2FC>0.5)
genes <- markers$gene

gene_sets <- list(
  "Pre-neoplastic Stress Response Gene Signature" = c("SCNN1G","H2AFX","ATR","CHEK1","ATM","CCND1","CCNE1","CDKN1A","CDK2","CDKN1B",
                                                      "CYP1A1","CYP1B1","HDAC1","HDAC2","MMP2","MKI67","JUN","FOS"),
  CNV_high_Gene_Program=genes,
  POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION = go %>%
    filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION") %>%
    pull(gene_symbol) %>%
    unique(),
  
  RESPONSE_TO_MECHANICAL_STIMULUS = go %>%
    filter(gs_name == "GOBP_RESPONSE_TO_MECHANICAL_STIMULUS") %>%
    pull(gene_symbol) %>%
    unique(),
  
  ACTIN_FILAMENT_ORGANIZATION = go %>%
    filter(gs_name == "GOBP_ACTIN_FILAMENT_ORGANIZATION") %>%
    pull(gene_symbol) %>%
    unique(),
  
  CELL_SUBSTRATE_ADHESION = go %>%
    filter(gs_name == "GOBP_CELL_SUBSTRATE_ADHESION") %>%
    pull(gene_symbol) %>%
    unique()
  
)

Epithelial_cells<-EPI
for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  Epithelial_cells <- AddModuleScore(Epithelial_cells, features = list(genes), name = set_name)
}
score_data <- FetchData(Epithelial_cells, 
                        vars = c("diseasestate", paste0(names(gene_sets), "1"))) 

score_data$diseasestate <- factor(score_data$diseasestate, levels = c("Ctrl","EuE","EcO","EcOA","CA"))

comparisons <- list(
  c("Ctrl", "EcO"),
  c("Ctrl", "EcOA"),
  c("Ctrl", "CA")
)

score_data <- score_data[!is.na(score_data$Diseasestate), ]

for (set_name in names(gene_sets)) {
  score_col <- paste0(set_name, "1")
  
  p <- ggplot(score_data, aes(x = Diseasestate, y = .data[[score_col]], fill = Diseasestate)) +
    geom_violin(, width = 1,trim = F, color = NA) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                       label = "p.forms", tip.length = 0.01) +
    scale_y_continuous(limits = c(0, NA))+
    theme_classic(base_size = 14) +
    labs(title = paste0(set_name),#, " Score in Epi"),
         y = "Module Score", x = "") +
    scale_fill_manual(values = c("#8dd3c7","#E4C755","#80b1d3","#bebada","#fb8072")) +
    theme_classic(base_size = 12) +   # 换成classic去掉背景线
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      , legend.position = "none"
    )
  
  
  # 保存图像
  ggsave(filename = paste0("Violin_", set_name, "_Epithelial_by_Tissue.pdf"),
         plot = p, width = 5, height = 4)
}

#Fig.3c,f,h
Epi2<-subset(EPI,Diseasestate%in%c("Ctrl","EcO","EcOA","CA"))
Epi2$diseasestate<-factor(Epi2$Diseasestate,levels = c("Ctrl","EcO","EcOA","CA"))
p151 <- DotPlot(Epi2, features = c("FOS","JUN","WFDC2"),group.by = "Diseasestate") + RotatedAxis()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.title = element_blank(),axis.line = element_blank())
ggsave("Fig.3c.pdf",p151,width =4.5,height = 3.5)

Epi4<-subset(EPI,Diseasestate%in%c("Ctrl","EcOA"))
Epi4$diseasestate<-factor(Epi4$Diseasestate,levels = c("Ctrl","EcOA"))
p151 <- DotPlot(Epi2, features = c("YAP1","CCN1","CCN2"),group.by = "Diseasestate") + RotatedAxis()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.title = element_blank(),axis.line = element_blank())
ggsave("Fig.3f.pdf",p151,width =5,height = 3.5)

p151 <- DotPlot(Epi3, features = c("SPHK2","S1PR2"),group.by = "Diseasestate") + RotatedAxis()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.title = element_blank(),axis.line = element_blank())+
  labs(title = "S1P–S1P2")
ggsave("S1P2.pdf",p151,width =4,height = 3.5)
p151 <- DotPlot(Epi3, features = c("RHOA","ROCK1","ROCK2","MYL9","MYL12A","MYL12B"),group.by = "Diseasestate") + RotatedAxis()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5,size=12),
        axis.title = element_blank(),axis.line = element_blank())+
  labs(title = "RHOA–Myosin")
ggsave("RHOA.pdf",p151,width =6.5,height = 3.5)

#Fig.3e
library(org.Hs.eg.db)
library(clusterProfiler)
Idents(scRNA)<-"diseasestate"
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, logfc.threshold =0.2,min.pct=0.2)
write.csv(scRNA.markers,file="markers_diseasestate.csv")

markers <- read.csv("markers_diseasestate.csv")
markers <- subset(markers, p_val<0.05 & avg_log2FC>0.25)

for(cluster2 in unique(markers$cluster)){
  input_gene = as.vector(subset(markers,cluster==cluster2)[,1])
  entrezIDS = mget(input_gene,org.Hs.egSYMBOL2EG,ifnotfound = NA)
  entrezIDS = as.character(entrezIDS)
  gene=entrezIDS[entrezIDS!="NA"]
  gene=gsub("c\\(\"(\\d+)\".*","\\1",gene)
  pvalueFilter=0.05
  qvalueFilter=0.1
  if(pvalueFilter==0.05){colorSel="pvalue"}else{colorSel="qvalue"}
  
  #KEGG
  kk<-enrichKEGG(gene = gene,organism = "hsa")
  KEGG=as.data.frame(kk)
  KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(input_gene[match(strsplit(x, "/")[[1]],as.character(entrezIDS))],collapse = "/")))
  write.table(KEGG,file=paste0(cluster2,".KEGG.txt"),sep="\t",quote = F,row.names = F)
  showNum=20
  pdf(file=paste0(cluster2,".KEGGbarplot.pdf"),width = 10,height = 7)
  bar=barplot(kk,drop=TRUE,showCategory=showNum,label_format=130,color=colorSel)
  print(bar)
  dev.off()
}

#Fig.3g,Fig.S3b


calc_pseudotime <- function(pc_df, splinefit) {
  min_x <- min(splinefit$PC1)
  max_x <- max(splinefit$PC1)
  
  pseudotime <- numeric(nrow(pc_df))
  
  for (i in 1:nrow(pc_df)) {
    x <- pc_df[i, "PC1"]
    y <- pc_df[i, "PC2"]
    distances <- sqrt((splinefit$PC1 - x)^2 + (splinefit$PC2 - y)^2)
    min_idx <- which.min(distances)
    pseudotime[i] <- ((splinefit$PC1[min_idx] - min_x) / (max_x - min_x)) * 100
  }
  
  pc_df$pseudotime <- pseudotime
  return(pc_df)
}

aggregate_pseudotime <- function(pc_df) {
  pc_df %>%
    group_by(sample) %>%
    summarise(mean_pseudotime = mean(pseudotime, na.rm = TRUE))
}

correlate_pseudotime_celltypes <- function(sample_pseudo, celltype_prop, method="spearman") {
  merged_df <- merge(sample_pseudo, celltype_prop, by = "sample")
  celltypes <- colnames(merged_df)[-(1:2)]
  
  cor_res <- data.frame(CellType = character(),
                        Correlation = numeric(),
                        Pvalue = numeric(),
                        stringsAsFactors = FALSE)
  
  for (ct in celltypes) {
    test <- cor.test(merged_df$mean_pseudotime, merged_df[[ct]], method = method)
    cor_res <- rbind(cor_res,
                     data.frame(CellType = ct,
                                Correlation = test$estimate,
                                Pvalue = test$p.value))
  }
  return(list(cor_res=cor_res, merged_df=merged_df))
}

plot_correlation_heatmap <- function(cor_res) {
  mat <- as.matrix(cor_res[, "Correlation", drop=FALSE])
  rownames(mat) <- cor_res$CellType
  pheatmap(mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Correlation between pseudotime and cell type proportions")
}

pc_df <- calc_pseudotime(pc_df, splinefit)
pc_df$sample <- substr(rownames(pc_df), 12, nchar(rownames(pc_df)))

sample_pseudo <- aggregate_pseudotime(pc_df)

celltype_prop <- read_excel("msc_prop.table_idents.xlsx")
celltype_prop <- celltype_prop %>%
  dplyr::select(sample = Var2, celltype = Var1, Freq) %>%
  pivot_wider(names_from = celltype, values_from = Freq)

res <- correlate_pseudotime_celltypes(sample_pseudo, celltype_prop, method="spearman")

res$merged_df$DiseaseState <- ifelse(
  grepl("CA", res$merged_df$sample), "CA",  
  ifelse(
    grepl("EcOA", res$merged_df$sample), "EcOA",
    ifelse(
      grepl("EcO", res$merged_df$sample), "EcO",
        "Normal"
      )
    )
  )

merged_df <- res$merged_df
merged_df$DiseaseState <- ifelse(
  grepl("CA", merged_df$sample), "CA",  
  ifelse(
    grepl("EcOA", merged_df$sample), "EcOA",
    ifelse(
      grepl("EcO", merged_df$sample), "EcO",
        "Ctrl"
      )
    )
  )
merged_df$diseaseState<-factor(merged_df$diseaseState,levels = c("Ctrl","EcO","EcOA","CA"))
# 获取细胞类型列
celltypes <- colnames(merged_df)[!(colnames(merged_df) %in% c("sample", "mean_pseudotime", "diseaseState"))]

# 循环绘图
for (ct in celltypes) {
  
  # 计算 Spearman 相关性
  cor_test <- cor.test(merged_df$mean_pseudotime, merged_df[[ct]], method = "spearman")
  corr_val <- round(cor_test$estimate, 2)
  p_val <- signif(cor_test$p.value, 2)
  
  # 注释文本
  annot_text <- paste0("Spearman r: ", corr_val, "\nP-value: ", p_val)
  
  # 绘图
  p <- ggplot(merged_df, aes_string(x = "mean_pseudotime", y = ct)) +
    geom_point(aes(color = DiseaseState), size = 3, alpha = 0.8) +
    scale_color_manual(values=c("#8dd3c7","#80b1d3","#bebada","#fb8072"),
                       name = "") +
    labs(x = "Malignancy continuum", y = paste0(ct, " proportion")) +
    ggtitle(ct) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      aspect.ratio = 1,
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
    ) +
    coord_cartesian(clip = 'off') +
    annotate("text", x = max(merged_df$mean_pseudotime, na.rm = TRUE)*0.7,
             y = max(merged_df[[ct]], na.rm = TRUE)*0.9,
             label = annot_text,hjust = 0.3,vjust=1,size = 4)# hjust = 0.5,vjust=4
  
  ggsave(paste0(ct,"_pseudotime_correlation.pdf"),p,width = 4.5,height = 3.5,dpi=800)
}

#Fig.3i,Fig.S3b
epi_yap <- FetchData(
  EPI,
  vars = c("orig.ident", "diseasestate", "YAP_Score1")
)
epi_yap_sample <- epi_yap %>%
  group_by(orig.ident, diseasestate) %>%
  summarise(
    YAP_score = mean(YAP_Score1, na.rm = TRUE),
    YAP_sd = sd(YAP_Score1, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )
cor_data <- merge(
  pc_sample,
  epi_yap_sample,
  by.x = "sample",
  by.y = "orig.ident"
)
cor_data <- cor_data %>%
  filter(DiseaseState %in% c("Ctrl", "EcO", "EcOA", "CA"))
cor_test <- cor.test(
  cor_data$Normalized_Rank,
  cor_data$YAP_score,
  method = "spearman",
  exact = FALSE
)

corr_val <- round(as.numeric(cor_test$estimate), 2)
p_val <- signif(cor_test$p.value, 2)

annot_text <- paste0(
  "Spearman ρ = ", corr_val,
  "\nP = ", p_val
)
p <- ggplot(
  cor_data,
  aes(
    x = Normalized_Rank,
    y = YAP_score,
    color = DiseaseState,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  geom_text(
    size = 3,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Ctrl" = "#8dd3c7",
      "EcO"  = "#80b1d3",
      "EcOA" = "#bebada",
      "CA"   = "#fb8072"
    ),
    name = ""
  ) +
  annotate(
    "text",
    x = min(cor_data$Normalized_Rank),
    y = max(cor_data$YAP_score),
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold",
    label = annot_text
  ) +
  labs(
    x = "Epithelial progression score",
    y = "Epithelial YAP activity"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )
ggsave("Epithelial_progression_vs_YAP.pdf",p,width = 5,height = 3.8)

msc<-readRDS("MSCs.rds")
MSC<-subset(msc,cell_type=="MSC_NR4A3")
msc_data <- GetAssayData(MSC, assay = "RNA", layer = "data")["MIF", ]
msc_pseudobulk <- data.frame(
  sample = MSC$orig.ident,
  state  = MSC$Diseasestate,
  expr   = as.numeric(msc_data)
) %>%
  group_by(sample, state) %>%
  summarise(pseudobulk = mean(expr), .groups = "drop")

epi_genes <- c("YAP1","CCN1","CCN2")

epi_data <- GetAssayData(EPI, assay = "RNA", layer = "data")[epi_genes, ]

epi_pseudobulk <- data.frame(
  sample = EPI$orig.ident,
  state  = EPI$Diseasestate,
  expr   = colMeans(epi_data) 
) %>%
  group_by(sample, state) %>%
  summarise(pseudobulk = mean(expr), .groups = "drop")

cor_data <- merge(
  msc_pseudobulk,
  epi_pseudobulk,
  by = c("sample", "state"),
  suffixes = c("_MIF", "_YAP")
)

cor_data <- subset(cor_data, state %in% c("CA","Ctrl","EcOA","EcO"))

cor_test<-cor.test(cor_data$pseudobulk_MIF, cor_data$pseudobulk_YAP, method = "spearman")
corr_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 2)
annot_text <- paste0("r = ", corr_val, "\nP = ", p_val)
p<-ggplot(cor_data, aes(pseudobulk_MIF, pseudobulk_YAP, color = state, label = sample)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values=c("#8dd3c7","#80b1d3","#bebada","#fb8072"),
                     name = "") +
  labs(
    x = "Expression level of MIF",
    y = "Expression level of YAP target genes",
    title = "MIF in MSC_NR4A3 / YAP target genes in Epi"
  )+
  theme_classic()+
  annotate("text", x = 0.5,
           y = 1.5,hjust = 0,vjust=1,size = 4,face="bold",
           label = annot_text)+
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(hjust = 0.5,size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
  )
ggsave("msc-mif_epi-yap.pdf",p,width = 5,height=3.5)