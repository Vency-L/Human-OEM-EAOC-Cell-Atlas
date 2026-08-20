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

# Read in RNA data for full cohort
scRNA <- readRDS("scRNA.rds")

scRNA$diseasestate<-factor(scRNA$diseasestate,levels=c("Ctrl","EuE","EcO","EcOA","CA"))
scRNA$orig.ident=factor(scRNA$orig.ident,c("Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5",
                                           "EuE1","EuE2","EuE3","EuE4","EuE5",
                                           "EcO1","EcO2","EcO3","EcO4","EcO5",
                                           "EcOA1","EcOA2","EcOA3","EcOA4","EcOA5",
                                           "CA1","CA2","CA3"))

#Fig.1b
cols<-c("#67001f","#e2796e","#f4a582",
        "#e7d36d","#993404","#97cac1", "#3860a0","#8c96c6","darkgrey")
celltype_levels <- c("FibC7", "MSC", "eF","Endo-vascular", "Endo-lymphatic",
                     "Epithelial","Lymphoid", "Myeloid", "Cycling")
celltype_ids <- setNames(seq_along(celltype_levels), celltype_levels)
plotData <- as.data.frame(scRNA[["umap"]]@cell.embeddings)
plotData$celltype <- factor(scRNA@meta.data$cell_type, levels = celltype_levels)
scRNA@meta.data$celltype_id <- celltype_ids[as.character(scRNA@meta.data$cell_type)]
plotData$cluster <- scRNA@meta.data$celltype_id
label <- plotData %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  ) %>%
  mutate(label_text = paste0(celltype_ids[celltype]))

p1<-ggplot(plotData, aes(x = umap_1, y = umap_2, fill = celltype, color = celltype)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 6)))+
  geom_text(data = label,
            aes(x = umap_1, y = umap_2, label = label_text),
            fontface = "bold",
            color = 'black', size = 3.5)
arrow_style <- arrow(length = unit(0.2, "cm"), type = "closed")
umap_coord <- ggplot(
  tibble(
    group = c("UMAP_1", "UMAP_2"),
    x = c(0, 0), xend = c(1, 0),
    y = c(0, 0), yend = c(0, 1),
    lx = c(0.55, -0.15), ly = c(-0.15, 0.55),
    angle = c(0, 90)
  )
) +
  geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
               arrow = arrow_style, size = 0.8, lineend = "round") +
  geom_text(aes(lx, ly, label = group, angle = angle),
            size = 2.5) +#, fontface = "bold"
  theme_void() +
  coord_fixed(xlim = c(-0.3, 1.2), ylim = c(-0.3, 1.2))

p_final <- ggdraw() +
  draw_plot(p1, 0, 0.02, 1, 1) +             
  draw_plot(umap_coord,
            x = -0.02,  
            y = -0.02,  
            width = 0.25,  
            height = 0.25) 
ggsave("celltype.png", p_final, width = 4, height = 3.5)
#Fig.S1b
p <- DimPlot(scRNA,reduction="umap",group.by = "cell_type",split.by = "diseasestate",cols =cols,pt.size = 1)
ggsave("celltype_ds.png",p,width =15,height = 4)

#Fig.1c
scRNA$cellclass<-factor(scRNA$cellclass,levels = c("Mesenchymal","Endothelial","Epithelial","Immune","Cycling"))
markers = c("DCN","PDGFRA","ACTA2",
            "VWF","PECAM1",
            "EPCAM", "CLDN3", "ASRGL1",
            "PTPRC",
            "MKI67", "TOP2A")
Vln.cols = c("#e2796e","#e7d36d","#97cac1", "#3860a0")
Vln.idents.cols = c(rep("#e2796e", 3), rep("#e7d36d", 2), 
                    rep("#97cac1", 3), rep("#3860a0", 1),
                    rep("darkgrey", 2))
p <- VlnPlot(scRNA, features = markers,group.by = "cellclass", 
             pt.size = 0, cols =Vln.idents.cols,
             stack = T, flip = T) + 
  NoLegend() +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 60,size=14),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(50, 20, 40, 20)
  )
ggsave("markers_by_cellclass.png",p,width = 4, height = 5.5)

#Fig.1d
cell_prop <- as.data.frame(prop.table(table(scRNA@meta.data$cell_type,scRNA@meta.data$diseasestate)))
colnames(cell_prop) <- c("celltype","tissue","proportion")
p = ggplot(cell_prop,aes(tissue,proportion,fill=celltype))+
  geom_bar(stat="identity",position = "fill")+
  guides(fill = guide_legend(title = NULL))+
  ggtitle("")+
  theme_bw()+ 
  scale_fill_manual(values = cols)+
  theme(plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  labs(title = "Fraction of Clusters")
ggsave("celltype_prop.png",p,width = 4,height=3)

#Fig.1e,Fig.S1c
# Run t.test on all celltypes of the subset
t.test_celltype <- function(celltype.x, prop.df, save.label, set.y.max = NULL) {
  
  # Subset the specific celltype
  prop.df = prop.df[prop.df$Var1 == celltype.x,]
  
  # Shorten the sample ID
  prop.df$Var2 <- sub("^.*_(KS\\d+)_.*$", "\\1", prop.df$Var2)
  
  # Run Kruskal-Wallis test
  kruskal.x = kruskal.test(Freq ~ sample.group, data = prop.df)
  
  kruskal.x$p.value
  
  # Perform Pairwise Wilcoxon Rank Sum Tests between groups
  wilcox.x <- pairwise.wilcox.test(prop.df$Freq, 
                                   prop.df$sample.group, 
                                   p.adjust.method="none",
                                   exact = TRUE)$p.value
  wilcox.x.padj <- pairwise.wilcox.test(prop.df$Freq, 
                                        prop.df$sample.group, 
                                        p.adjust.method="BH",
                                        exact = TRUE)$p.value
  
  # Save the data in a dataframe to be returned
  wilcox.df = data.frame(expand.grid(dimnames(wilcox.x)), array(wilcox.x), array(wilcox.x.padj))
  wilcox.df$Celltype = celltype.x
  
  print(paste("p-value = ", wilcox.df[1,3], "in", celltype.x))
  

  comp.order <- list(c("EcO", "EuE"),c("EcO", "EcOA"), c("Ctrl", "EcO"), c("Ctrl", "EcOA"),c("Ctrl", "CA"))
  # comp.order <- ""
  # Change freq to Percentage
  prop.df$Percentage = prop.df$Freq*100
  
  # Get max proportion and round to single digit for Y-axis limit if not set
  if (is.null(set.y.max) == TRUE) {
    y.max = max(prop.df$Percentage) + 5
  } else if (!is.null(set.y.max) == TRUE) {
    y.max = set.y.max
  }
  
  plot.x = ggbarplot(prop.df, x = "sample.group", y = "Percentage",
                     fill = "sample.group", palette = Group.cols, 
                     order = c("Ctrl", "EuE", "EcO", "EcOA","CA"),legend = "none",
                     add = "mean_sd",xlab = "", 
                     ylab = "Cell type proportion (%)", #ylim = c(0,y.max),
                     size = 0.5, width = 0.9,main = celltype.x) +
    stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                       label = "p.format", label.y = y.max-2.5,
                       bracket.size = 0.5) + theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1),legend.title=element_text("none"),
          plot.title = element_text(hjust = 0.5), legend.position = "none")
  
  plot.y = ggboxplot(prop.df, x = "sample.group", y = "Percentage",
                     fill = "sample.group", palette = Group.cols, 
                     order = c("Ctrl", "EuE", "EcO", "EcOA","CA"),
                     add = "dotplot",
                     bxp.errorbar = FALSE,legend = "none",
                     #legend="right"
                     xlab = "",ylab = "", #ylim = c(0,y.max), 
                     size = 0.5, width = 0.9,main = celltype.x) +
    stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                       label = "p.format",
                       bracket.size = 0.5) + theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.title=element_text("none"), legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  ggsave2(plot = plot.y, filename = paste0(set.ident,"_", celltype.x, "_", save.label, "_boxplot_wilcox.png"), 
          dpi = 700, width = 3,height=4.5)
  
  #saveRDS(object = wilcox.x$p.value, file = paste0(set.ident,"_", celltype.x, "_", save.label, "_wilcox_test_padj.rds"))
  return(wilcox.df)
  
}
wilcox.list = lapply(order.celltypes, function(i) t.test_celltype(celltype.x = i, prop.df = prop.df_idents, save.label = set.ident, set.y.max = 100))

#Fig.1f
library(Startrac)
idents.order = c("Ctrl","EuE","EcO","EcOA","CA")
data <- scRNA@meta.data
data$diseasestate <- factor(data$diseasestate, levels = idents.order)
data$cell_type <- droplevels(factor(data$celltype))
R_oe <- calTissueDist(
  data,
  byPatient = FALSE,
  colname.cluster="celltype",
  colname.patient="orig.ident",
  colname.tissue="diseasestate",
  method="chisq",
  min.rowSum=0
)

pdf(paste0("roe.pdf"), width = 3, height = 4.5)
p <- Heatmap(
  as.matrix(R_oe),
  name = "Ro/e",                      
  column_title = "",                
  column_title_gp = gpar(fontsize=12, fontface="bold"),
  show_heatmap_legend = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = 'right',
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = c("white","#E95C59"),
  row_names_gp = gpar(fontsize=10),
  column_names_gp = gpar(fontsize=10),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Ro/e",
    at = seq(0, 2, by=0.5),
    labels = seq(0, 2, by=0.5),
    legend_gp = gpar(fill=c("white","#E95C59"))
  ),
  cell_fun = function(j,i,x,y,width,height,fill){
    grid.text(sprintf("%.2f", R_oe[i,j]), x, y, gp=gpar(fontsize=8))
  }
)
draw(p)
dev.off()

#Fig.1g
comparisons <- list(
  c("Ctrl", "EcOA"),
  c("Ctrl", "EcO"),
  c("Ctrl", "CA")
)
for (pair in comparisons) {
  group1 <- pair[1]
  group2 <- pair[2]
  sub_scRNA <- subset(scRNA, diseasestate %in% c(group1, group2))
  augur_result <- readRDS(paste0("augur_", group1, "_vs_", group2, ".rds"))
  
  plot_l <- plot_lollipop(augur_result) +
    geom_segment(aes(xend = cell_type, yend = 0.5), size = 1) +
    labs(title=paste0(group1, " vs ", group2))+
    theme(legend.position = "none",
          axis.text.y = element_text(size=8),
          plot.title = element_text(size=10),
          axis.text.x = element_text(size=8),
          geom = element_geom(fontsize = 30,color = c("FibC7"="#67001f",
                                                      "MSCs"="#e2796e",
                                                      "eF"="#f4a582",
                                                      "Endo-vascular"="#e7d36d",
                                                      "Endo-lymphatic"="#993404",
                                                      "Epithelial"="#97cac1", 
                                                      "Lymphoid"="#3860a0",
                                                      "Myeloid"="#8c96c6",
                                                      "Cycling"="darkgrey"),pointsize = 20))
  pdf_file1 <- paste0("Lollipop_augur_", group1, "_vs_", group2, ".pdf")
  ggsave(pdf_file1, width = 4, height = 3, plot = plot_l)
}


