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

#Fig.7a-b
scRNA <- readRDS("myeloid_mouse.rds")#"lympoid_mouse.rds"
scRNA$group<-factor(scRNA$group,levels=c("CTRL","MIF"))

cols<-c("#8E24AA","#7B3FA7","#AB47BC","#CE93D8","#E6C1EC",
        "#9e9ac8","#9779C3", "#6a51a3","#5E3A8C","#542788","#45206E")#"#bdd7e7","#6baed6","#2171b5","#08519c","#08306b","#3b6cb0","#4292c6","#9ecae1","#41b6c4","#1b9e77","#66c2a5","#99d8c9","#a6dba0","#d9f0d3"
celltype_levels <- c("Macro_C1qa","Macro_Trem2","Macro_Mrc1","Macro_Arg1","Macro_Cxcl9","Macro_Gpnmb",
                     "DC_Cd209a","DC_Xcr1","DC_Ccr7","pDC","Neutrophils")#"Cd4_Treg_Foxp3","Cd4_Tn_Ccr7","Cd8_T_Ccl5","T_Cd247","T_ISG","NKT_Zbtb16","NK_Klra4","NK_Gzmc","NK_Ncr1","NK_Prf1","T_NK_cycling","gdT17_Rorc","ILC2_Gata3","B"
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

#Fig.7c
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
                     ylab = "Cell type proportion (%)",
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
                     xlab = "",ylab = "",
                     size = 0.5, width = 0.9,main = celltype.x) +
    stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                       label = "p.format",
                       bracket.size = 0.5) + theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.title=element_text("none"), legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  ggsave2(plot = plot.y, filename = paste0(set.ident,"_", celltype.x, "_", save.label, "_boxplot_wilcox.png"), 
          dpi = 700, width = 3,height=4.5)

  return(wilcox.df)
  
}
wilcox.list = lapply(order.celltypes, function(i) t.test_celltype(celltype.x = i, prop.df = prop.df_idents, save.label = set.ident, set.y.max = 100))

#Fig.7d-e
Idents(scRNA)<-"celltype"
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, logfc.threshold =0.2,min.pct=0.2)
write.csv(scRNA.markers,file="cell_identify/single_cell_marker/markers_celltype.csv")

markers <- read.csv("cell_identify/single_cell_marker/markers_celltype.csv")
markers <- subset(markers, p_val<0.05 & avg_log2FC>0.25)
custom_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black")
)
for(cluster2 in unique(markers$cluster)){

  input_gene = as.vector(subset(markers, cluster == cluster2)[, 1])

  entrezIDS = mget(input_gene, org.Mm.egSYMBOL2EG, ifnotfound = NA)
  entrezIDS = as.character(entrezIDS)
  gene = entrezIDS[entrezIDS != "NA"]
  gene = gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  
  pvalueFilter = 0.05
  qvalueFilter = 0.1

  kk = enrichGO(gene = gene, 
                OrgDb = org.Mm.eg.db, 
                ont = "ALL", 
                readable = TRUE, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1)
  GO = as.data.frame(kk)
  GO = GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ]
  write.table(GO, file = paste0(cluster2, ".GO.txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  showNum = 10

  pdf(file = paste0(cluster2, ".GObarplot.pdf"), width = 7.5, height = 8)
  bar = barplot(kk, drop = TRUE, showCategory = showNum, 
                label_format = 130, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~., scales = "free") +
    custom_theme
  print(bar)
  dev.off()

  pdf(file = paste0(cluster2, ".GObubble.pdf"), width = 7.5, height = 8)
  bub = dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
                label_format = 130, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~., scales = "free") +
    custom_theme
  print(bub)
  dev.off()

  kk = enrichKEGG(gene = gene, 
                  organism = "mmu",   
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1)
  KEGG = as.data.frame(kk)

  if(nrow(KEGG) > 0){
    KEGG$geneID = as.character(sapply(KEGG$geneID, function(x) {
      paste(input_gene[match(strsplit(x, "/")[[1]], as.character(entrezIDS))], 
            collapse = "/")
    }))
  }
  write.table(KEGG, file = paste0(cluster2, ".KEGG.txt"), sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  showNum = 20

  pdf(file = paste0(cluster2, ".KEGGbarplot.pdf"), width = 9, height = 7.5)
  bar = barplot(kk, drop = TRUE, showCategory = showNum, 
                label_format = 130) +
    custom_theme
  print(bar)
  dev.off()

  pdf(file = paste0(cluster2, ".KEGGbubble.pdf"), width = 9, height = 7.5)
  bub = dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
                label_format = 130) +
    custom_theme
  print(bub)
  dev.off()
}

#Fig.7f
macro = subset(scRNA,celltype %in% c("Macro_C1qa","Macro_Trem2","Macro_Mrc1",
                                     "Macro_Arg1","Macro_Cxcl9","Macro_Gpnmb"))
gene_sets <- list(
  M1_like = c(
    "Il1b", "Il6",  "Cxcl9",  "Nos2", "Cd40", "Cd86", "Fcgr4", 
    "Irf5", "Fpr1", "Fcgr1", "Cd80", "Irf1", "Ccl5","Cd38","Fpr2","Gpr18"
  ),
  M2_like = c(
    "Arg1","Arg2","Il10","Cd163","Cd200r1","Pdcd1lg2","Tgfb1","Tgfb2","Tgfb3",
    "Fn1","Ccl4","Myc","Egr2","Vegfa","Il1rn","Il1r2","Il4r","Mrc1","Chi3l3","Retnla","Irf4","Trem2"
  )
)
for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  macro <- AddModuleScore(macro, features = list(genes), name = set_name)
}

score_data <- FetchData(macro, 
                        vars = c("group", paste0(names(gene_sets), "1")))
score_data$group <- factor(score_data$group, levels = c("CTRL","MIF"))

comparisons <- list(
  c("CTRL", "MIF")
)

score_data <- score_data[!is.na(score_data$group), ]

for (set_name in names(gene_sets)) {
  score_col <- paste0(set_name, "1")
  
  p <- ggplot(score_data, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_violin(, width = 1,trim = F, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",color="black") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                       label = "p.form", tip.length = 0.01) +
    scale_y_continuous(limits = c(0, NA))+
    theme_classic(base_size = 15) +
    labs(title = set_name,
         y = "Module Score", x = "") +
    scale_fill_manual(values = c("#1F77B4","#FF7F0E")) +
    theme_classic(base_size = 12) +  
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.position = ""
    )

  ggsave(filename = paste0("Violin_", set_name, "_macro_by_Tissue.png"),
         plot = p, width = 2.5, height = 3.3,dpi = 600)
}

t = subset(scRNA,celltype %in% c("Cd4_Treg_Foxp3","Cd4_Tn_Ccr7","Cd8_T_Ccl5","T_Cd247","T_ISG",
                                 "NKT_Zbtb16","gdT17_Rorc"))
gene_sets <- list(
  Cytotoxic = c("Prf1","Gzma","Gzmk","Nkg7","Gzmb",
                "Ifng","Tnfsf10","Eomes","Cx3cr1","Ccl3", "Ccl4", "Ccl5", "Xcl1", "Cxcr3"),
  Exhausted = c("Pdcd1", "Ctla4", "Havcr2", "Lag3", "Maf")
)
for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  t <- AddModuleScore(t, features = list(genes), name = set_name)
}

score_data <- FetchData(t, 
                        vars = c("group", paste0(names(gene_sets), "1")))
score_data$group <- factor(score_data$group, levels = c("CTRL","MIF"))

comparisons <- list(
  c("CTRL", "MIF")
)

score_data <- score_data[!is.na(score_data$group), ]

for (set_name in names(gene_sets)) {
  score_col <- paste0(set_name, "1")
  
  p <- ggplot(score_data, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_violin(, width = 1,trim = F, color = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",color="black") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                       label = "p.form", tip.length = 0.01) +
    scale_y_continuous(limits = c(0, NA))+
    theme_classic(base_size = 15) +
    labs(title = set_name,
         y = "Module Score", x = "") +
    scale_fill_manual(values = c("#1F77B4","#FF7F0E")) +
    theme_classic(base_size = 12) +   # 换成classic去掉背景线
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.position = ""
    )

  ggsave(filename = paste0("Violin_", set_name, "_T_by_Tissue.png"),
         plot = p, width = 2.5, height = 3,dpi = 600)
}

#Fig.7g-h,Fig.S7e
m <- readRDS("D:/scRNA/mouse/total/myeloid/harmony_annotation.rds")
l <- readRDS("D:/scRNA/mouse/total/lymphoid/harmony_annotation.rds")
scRNA3<-merge(m,l,add.cell.ids = c("m","l"))
scRNA3[["RNA"]] <- JoinLayers(scRNA3[["RNA"]])
colnames(scRNA3@assays[["RNA"]]@layers[["data"]]) <- rownames(scRNA3@meta.data)

Idents(scRNA3) = "group"
scRNA.CTRL = scRNA3[,Idents(scRNA3) %in% c("CTRL")]
scRNA.MIF = scRNA3[,Idents(scRNA3) %in% c("MIF")]

seurat.list <- list(
  "CTRL"  = scRNA.CTRL,
  "MIF"   = scRNA.MIF
)

seurat.list <- lapply(seurat.list, function(x) { 
  Idents(x) <- "celltype"
  return(x)
})
celltypes.list <- lapply(seurat.list, function(x) levels(Idents(x)))
common.celltypes <- Reduce(intersect, celltypes.list)

CellChat.subset <- function(seurat.subset) {
  Idents(seurat.subset) <- "celltype"
  seurat.subset <- subset(seurat.subset, celltype %in% common.celltypes)
  CC.input <- GetAssayData(seurat.subset, assay = "RNA", slot = "data")
  CC.labels <- Idents(seurat.subset)
  CC.meta <- data.frame(labels = CC.labels, row.names = colnames(seurat.subset))
  CC.x <- createCellChat(object = CC.input, meta = CC.meta, group.by = "labels")
  CC.x@DB <- CellChatDB.mouse                     
  CC.x <- subsetData(CC.x)              
  CC.x <- identifyOverExpressedGenes(CC.x)
  CC.x <- identifyOverExpressedInteractions(CC.x)
  CC.x <- computeCommunProb(CC.x, 
                            population.size = TRUE, 
                            type = "triMean")      
  CC.x <- filterCommunication(CC.x, min.cells = 10)
  CC.x <- computeCommunProbPathway(CC.x)
  CC.x <- aggregateNet(CC.x)
  return(CC.x)
}

CC.list <- list()
for (group.name in names(seurat.list)) {
  seurat.subset <- seurat.list[[group.name]]
  DefaultAssay(seurat.subset) <- "RNA"
  
  print(paste("Running CellChat on", group.name))
  CC.save <- CellChat.subset(seurat.subset)
  
  df_comm <- CC.save@net$weight
  write.csv(df_comm, file = paste0("CellChat_LRpairs_", group.name, ".csv"), row.names = TRUE)
  
  saveRDS(CC.save, file = paste0("CellChat_object_", group.name, ".rds"))
  CellChat_object_group.name <- CC.save
  CC.list[[group.name]] <- CC.save
}

CC.list<-list(  "CTRL"  = CellChat_object_CTRL,
                "MIF"   = CellChat_object_MIF
)
CC.merged <- mergeCellChat(CC.list[c("CTRL","MIF")], 
                           add.names = c("CTRL","MIF"))

comparisons <- list(
  "CTRL_vs_MIF"  = c(1,2)
)

for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  
  pdf(paste0("compare_ranknet_", name, ".pdf"), width = 4.5, height = 4)
  g1 <- rankNet(CC.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE, 
                comparison = comp,color.use = c("#1F77B4","#FF7F0E"),
                sources.use = c("Macro_Trem2","Macro_Mrc1",
                                "Macro_Arg1","Macro_Cxcl9","Macro_Gpnmb"),
                targets.use = "CD4_Treg_Foxp3",
                font.size = 12)
  g2 <- rankNet(CC.merged, mode = "comparison", stacked = FALSE, do.stat = TRUE, 
                comparison = comp,color.use = c("#1F77B4","#FF7F0E"),
                sources.use = c("Macro_Trem2","Macro_Mrc1",
                                "Macro_Arg1","Macro_Cxcl9","Macro_Gpnmb"),
                targets.use = "CD4_Treg_Foxp3",
                font.size = 12)
  print(g1)
  print(g2)
  dev.off()
  
  p <- netVisual_bubble(
    CC.merged, 
    sources.use = c("Macro_Trem2","Macro_Mrc1",
                    "Macro_Arg1","Macro_Cxcl9","Macro_Gpnmb"), 
    targets.use = "CD4_Treg_Foxp3",
    color.text = c("#1F77B4","#FF7F0E"),
    comparison = comp, angle.x = 45, remove.isolate = F
  )
  ggsave(paste0("compare_bubble_", name, "_pathways.png"), p, width = 6.5, height = 4.5,dpi=800)

  # heatmap
  pdf(paste0("compare_heatmap_", name, ".pdf"), width = 6, height = 5.5)
  par(xpd=TRUE)
  print(netVisual_heatmap(CC.merged, comparison = comp,
  ))
  print(netVisual_heatmap(CC.merged, measure = "weight", comparison = comp
  ))
  dev.off()
}

