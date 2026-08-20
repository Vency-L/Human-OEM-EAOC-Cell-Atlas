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
scRNA <- readRDS("myeloid.rds")

scRNA$diseasestate<-factor(scRNA$diseasestate,levels=c("Ctrl","EuE","EcO","EcOA","CA"))
scRNA$orig.ident=factor(scRNA$orig.ident,c("Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5",
                                           "EuE1","EuE2","EuE3","EuE4","EuE5",
                                           "EcO1","EcO2","EcO3","EcO4","EcO5",
                                           "EcOA1","EcOA2","EcOA3","EcOA4","EcOA5",
                                           "CA1","CA2","CA3"))

#Fig.6a
cols<-c("#8E24AA","#AB47BC","#CE93D8","#E6C1EC",
        "#9e9ac8","#9779C3", "#6a51a3","#542788")
celltype_levels <- c("Macro_LYVE1","DC","pDC","Macro_APOE","Monocytes",
                     "Neutrophils","Mast","Myeloid_cycling")
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

#Fig.6b
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

#Fig.6c,Fig.S5a
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
  
  # Testing only EcOA vs. PCOS BL 
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
  
  return(wilcox.df)
  
}
wilcox.list = lapply(order.celltypes, function(i) t.test_celltype(celltype.x = i, prop.df = prop.df_idents, save.label = set.ident, set.y.max = 100))

#Fig.6d
p <- DotPlot(scRNA, "CD74",group.by = "diseasestate") + RotatedAxis()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),axis.line = element_blank())+
  labs(title = "CD74")+scale_x_discrete(labels = "CD74")
ggsave("CD74_UMAP.png",p,width =4,height = 4)

#Fig.6e
macro<-subset(scRNA,celltype %in% c("Macro_LYVE1","Macro_APOE"))

df <- FetchData(macro, vars = c("CD74", "celltype"))
p15 <- ggplot(df, aes(x = celltype, y = CD74, fill = celltype)) +
  geom_violin(color = NA, width = 1,trim = F) +   # violin更宽
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("Macro_LYVE1","Macro_APOE")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "CD74 Expression in Macro", x = "", y = "CD74 Expression") +
  scale_fill_manual(values = c("#8E24AA","#E6C1EC")) +
  theme_classic(base_size = 12) + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("CD74_macro.pdf",p15,width = 3,height=3)

#Fig.6f,Fig.S5b
scRNA2<-subset(scRNA,celltype %in% c("Macro_LYVE1"))#"Macro_APOE","DC"
df <- FetchData(scRNA2, vars = c("CD74", "Diseasestate"))
p15 <- ggplot(df, aes(x = Diseasestate, y = CD74, fill = Diseasestate)) +
  geom_violin(color = NA, width = 1,trim = F) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("Ctrl", "EcO"),
      c("Ctrl", "EcOA"),
      c("Ctrl", "CA")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "CD74 Expression in Macro_APOE", x = "", y = "CD74 Expression") +
  scale_fill_manual(values = c("#8dd3c7","#E4C755","#80b1d3","#bebada","#fb8072")) +
  theme_classic(base_size = 13) + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("CD74_Macro_LYVE1.pdf",p15,width = 4,height=4,dpi=800)

#Fig.6g
scRNA <- readRDS("MSC.rds")
scRNA2<-subset(scRNA,celltype%in%c("MSC_NR4A3"))
scRNA <- readRDS("myeloid.rds")
scRNA3<-subset(scRNA,celltype %in% c("Macro_LYVE1","Macro_APOE"))
scRNA3<-merge(scRNA2,scRNA3,add.cell.ids = c("MSC","Macro"))
scRNA3<-merge(msc,scRNA3,add.cell.ids = c("MSC","Macro"))
scRNA3[["RNA"]] <- JoinLayers(scRNA3[["RNA"]])
colnames(scRNA3@assays[["RNA"]]@layers[["data"]]) <- rownames(scRNA3@meta.data)
library(CellChat)
library(patchwork)
library(cowplot)
Idents(scRNA3) = "diseasestate"
scRNA.EcO = scRNA3[,Idents(scRNA3) %in% c("EcO")]
scRNA.EcOA = scRNA3[,Idents(scRNA3) %in% c("EcOA")]
scRNA.Ctrl = scRNA3[,Idents(scRNA3) %in% c("Ctrl")]

seurat.list <- list(
  "Ctrl"  = scRNA.Ctrl,
  "EcO"   = scRNA.EcO,
  "EcOA"   = scRNA.EcOA
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
  CC.x@DB <- CellChatDB.human
  CC.x <- subsetData(CC.x)
  CC.x <- identifyOverExpressedGenes(CC.x)
  CC.x <- identifyOverExpressedInteractions(CC.x)
  PPI_matrix <- as(PPI.human, "sparseMatrix")
  CC.x <- smoothData(CC.x, adj = PPI_matrix)
  CC.x <- computeCommunProb(CC.x, population.size = TRUE)
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

CC.list<-list(  "Ctrl"  = CellChat_object_Ctrl,
                "EcO"   = CellChat_object_EcO,
                "EcOA"   = CellChat_object_EcOA
)
CC.merged <- mergeCellChat(CC.list[c("Ctrl","EcO","EcOA")], 
                           add.names = c("Ctrl","EcO","EcOA"))
comparisons <- list(
  "Ctrl_vs_EcO"  = c(1,2),
  "Ctrl_vs_EcOA"  = c(1,3)
)

for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  message("Running comparison: ", name, " -> ", paste(comp, collapse = "_"))
  # circle plot
  pdf(paste0("compare_circle_", name, ".pdf"), width = 4.5, height = 3.5)
  par(xpd=TRUE)
  netVisual_diffInteraction(CC.merged, weight.scale = TRUE, comparison = comp,
                            color.use=c("#990000","#8E24AA","#E6C1EC"))
  netVisual_diffInteraction(CC.merged, weight.scale = FALSE, measure = "weight", comparison = comp,
                            vertex.weight = 20,
                            edge.width.max = 8,
                            color.use=c("#990000","#8E24AA","#E6C1EC"))
  dev.off()
}

#Fig.6h
scRNA3@meta.data$samples <- scRNA3@meta.data$orig.ident
cellchat<-createCellChat(scRNA3,meta=scRNA3@meta.data,
                         group.by="cell_type")
cellchatDB<-CellChatDB.human
showDatabaseCategory(cellchatDB)
cellchat@DB
cellchat@DB<-cellchatDB
cellchat<-subsetData(cellchat)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human)
cellchat<-computeCommunProb(cellchat,raw.use=F)
cellchat<-filterCommunication(cellchat,min.cells=10)
cellchat<-computeCommunProbPathway(cellchat)
df.net<-subsetCommunication(cellchat)
df.netP<-subsetCommunication(cellchat,slot.name="netP")
cellchat<-aggregateNet(cellchat)

pdf("MIF.pdf",width = 4,height = 3.5)
netVisual_heatmap(cellchat, signaling = "MIF", color.heatmap = "Reds",
                  color.use=c("#E6C1EC","#8E24AA","#990000"))
dev.off()

#Fig.6i
comparisons <- list(
  "Ctrl_vs_EcO"  = c(1,2)
)
for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  p <- netVisual_bubble(
    CC.merged, sources.use = "MSC_NR4A3",targets.use = c("Macro_LYVE1","Macro_APOE"), color.text = c("#8dd3c7","#80b1d3"),
    angle.x = 45,comparison = comp, remove.isolate = TRUE
  )
  ggsave(paste0("compare_bubble_", name, ".png"), p, width = 4, height = 8,dpi=600)
  bubble_data <- netVisual_bubble(CC.merged, sources.use = "MSC_NR4A3",targets.use = c("Macro_LYVE1","Macro_APOE"), 
                                  comparison = comp, angle.x = 45, return.data = TRUE)
  df_bubble <- as.data.frame(bubble_data$communication)
  write.csv(df_bubble, file = paste0("compare_bubble_", name, ".csv"), row.names = FALSE)
}

#Fig.6j
scRNA<-readRDS("myeloid.rds")
scRNA2<-subset(scRNA,celltype=="Macro_LYVE1")

cd74_expr <- FetchData(scRNA2, "CD74")[,1]
scRNA2$CD74_status <- cut(
  cd74_expr,
  breaks = quantile(cd74_expr, probs = c(0, 0.25, 0.75, 1)),
  labels = c("CD74_low","CD74_mid","CD74_high"),
  include.lowest = TRUE
)

p<-DimPlot(scRNA2, group.by = "CD74_status",size=0.01,cols = c("#619CFF", "#FFEBCD", "#F8766D"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(hjust = 0.5),aspect.ratio = 1,
        axis.line = element_blank())
ggsave("cell_identify/CD74_status_umap.png",p,width = 4,height=3,dpi = 600)

cell_prop <- as.data.frame(prop.table(table(scRNA2@meta.data$CD74_status,scRNA2@meta.data$diseasestate)))
colnames(cell_prop) <- c("CD74_status","tissue","proportion")
p = ggplot(cell_prop,aes(tissue,proportion,fill=CD74_status))+
  geom_bar(stat="identity",position = "fill")+
  guides(fill = guide_legend(title = NULL))+
  ggtitle("")+
  theme_bw()+ 
  scale_fill_manual(values = c("#619CFF", "#FFEBCD", "#F8766D"))+
  theme(plot.title = element_text(size = 14,face = "bold",hjust = 0.5),
        axis.ticks.length = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.title.x = element_blank()#,legend.position = "none"
  )+
  labs(title = "Fraction of Clusters")
ggsave("CD74_status_prop.pdf",p,width = 4,height=3,dpi = 600)

#Fig.6k
gene_sets <- list(
  Proinflammatory=c("IL1B","TNF","CCL2","CCL3","CCL5","CCL7","CCL8","CCL13","CCL17","CCL22"),
  Immunoregulatory=c("ARG1","ARG2","IL10","CD32","CD163","CD23","CD200R1","PDCD1LG2","CD274","MARCO","CSF1R","CD206","IL1RN","IL1R2","IL4R",
                     "CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CSTC","CTSD",
                     "TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4")
)

for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  scRNA2 <- AddModuleScore(scRNA2, features = list(genes), name = set_name)
}
score_data <- FetchData(scRNA2, 
                        vars = c("CD74_status", paste0(names(gene_sets), "1")))  # AddModuleScore后缀加1

score_data$CD74_status <- factor(score_data$CD74_status, levels = c("CD74_high","CD74_low"))

comparisons <- list(
  c("CD74_high","CD74_low")
)

score_data <- score_data[!is.na(score_data$CD74_status), ]

for (set_name in names(gene_sets)) {
  score_col <- paste0(set_name, "1")
  
  p <- ggplot(score_data, aes(x = CD74_status, y = .data[[score_col]], fill = CD74_status)) +
    geom_violin(, width = 1,trim = F, color = NA) +
    geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                       label = "p.forms", tip.length = 0.01) +
    scale_y_continuous(limits = c(0, NA))+
    theme_classic(base_size = 14) +
    labs(title = paste0(set_name),
         y = "Module Score", x = "") +
    scale_fill_manual(values = c("#F8766D", "#619CFF")) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(filename = paste0("Violin_", set_name, "_scRNA2.png"),
         plot = p, width = 3, height = 3)
}

#Fig.6l
Idents(scRNA2)<-"CD74_status"
comparisons <- list(
  c("CD74_low","CD74_high")
)

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste0(group1, "_vs_", group2)
 
  CELLDEG <- FindMarkers(scRNA2,
                         ident.1 = paste0(group2),
                         ident.2 = paste0(group1),
                         verbose = FALSE,
                         test.use = "wilcox",
                         min.pct = 0.2)#min.pct = 0.1)
  
  CELLDEG <- CELLDEG %>%
    mutate(Type = if_else(p_val_adj > 0.05, "ns",
                          if_else(abs(avg_log2FC) < 0.5, "ns",
                                  if_else(avg_log2FC >= 0.5, "up", "down")))) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    rownames_to_column("Gene_Symbol")
  
  CELLDEG$Gene_Symbol <- sub("\\..*", "", CELLDEG$Gene_Symbol)
  
  CELLDEG_sig <- CELLDEG %>% filter(Type %in% c("up", "down"))
 
  write.csv(CELLDEG, paste0(comp_name, "_all_deg.csv"), row.names = FALSE)
  write.csv(CELLDEG_sig, paste0(comp_name, "_sig_deg.csv"), row.names = FALSE)

  genes <- c("SIGLEC10","SERPING1","HLA-DQA1", "HLA-DQB1", "HLA-DRB1",  "HLA-DMA", "HLA-DRA","CD74","C1QA","C1QB","C1QC","CX3CR1","SERPINA1",
             "CD36","ATP5E","ATP5L","ATP5G2","FABP5","H2AFY","H3F3B","H3F3A","ATP5G3","ATP5I")
  genes_to_label <- CELLDEG %>% filter(Gene_Symbol %in% genes)
 
  pdf(paste0("Macro_DEG_volcano_", comp_name, ".pdf"), width = 8.5, height = 5.5)
  p=CELLDEG %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Type), alpha = 0.9, size = 2.5) +
    scale_color_manual(
      name = "Expression",
      values = c(ns = "grey70", up = "#F8766D", down = "#619CFF")#up = "#8E24AA", down = "#E6C1EC")
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
    geom_label_repel(
      data = genes_to_label,
      aes(label = Gene_Symbol, color = Type),
      size = 4,
      box.padding = 0.8,
      segment.size = 0.3,
      show.legend = FALSE,
      max.overlaps = Inf,
      seed = 18
    ) +
    labs(
      title = paste0("DEGs (", comp_name, ")"),
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("Adjusted p-value"))
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.1)),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  print(p)
  dev.off()
}
