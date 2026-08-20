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
scRNA <- readRDS("mesenchymal.rds")

scRNA$diseasestate<-factor(scRNA$diseasestate,levels=c("Ctrl","EuE","EcO","EcOA","CA"))
scRNA$orig.ident=factor(scRNA$orig.ident,c("Ctrl1","Ctrl2","Ctrl3","Ctrl4","Ctrl5",
                                           "EuE1","EuE2","EuE3","EuE4","EuE5",
                                           "EcO1","EcO2","EcO3","EcO4","EcO5",
                                           "EcOA1","EcOA2","EcOA3","EcOA4","EcOA5",
                                           "CA1","CA2","CA3"))

#Fig.1a
cols<-c("#990000","#E95C59", "pink")
celltype_levels <- c("MSC_NR4A3","MSC_STING","MSC_MMP11")
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

#Fig.2b,Fig.S2a
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

#Fig.2c
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

#Fig.2d-f
scRNA$celltype.tissue <- paste(scRNA$celltype, scRNA$diseasestate, sep = "_")
Idents(scRNA) <- "celltype.tissue"

comparisons <- list(
  c("EuE", "EcO"))

for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste0(group1, "_vs_", group2)
  
  CELLDEG <- FindMarkers(scRNA,
                         ident.1 = paste0("MSC_NR4A3_", group2),
                         ident.2 = paste0("MSC_NR4A3_", group1),
                         verbose = FALSE,
                         test.use = "wilcox",
                         min.pct = 0.2)
  
  CELLDEG <- CELLDEG %>%
    mutate(Type = if_else(p_val_adj > 0.05, "ns",
                          if_else(abs(avg_log2FC) < 1, "ns",
                                  if_else(avg_log2FC >= 1, "up", "down")))) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    rownames_to_column("Gene_Symbol")
  
  CELLDEG$Gene_Symbol <- sub("\\..*", "", CELLDEG$Gene_Symbol)
  
  CELLDEG_sig <- CELLDEG %>% filter(Type %in% c("up", "down"))
  
  write.csv(CELLDEG, paste0("MSC_NR4A3_", comp_name, "_all_deg.csv"), row.names = FALSE)
  write.csv(CELLDEG_sig, paste0("MSC_NR4A3_", comp_name, "_sig_deg.csv"), row.names = FALSE)
  
  top_genes <- bind_rows(
    CELLDEG %>%
      filter(Type == "up") %>%
      arrange(p_val_adj) %>%
      slice_head(n = 10),
    CELLDEG %>%
      filter(Type == "down") %>%
      arrange(p_val_adj) %>%
      slice_head(n = 10)
  )
  
  pdf(paste0("MSC_NR4A3_DEG_volcano_", comp_name, ".pdf"), width = 7, height = 5)
  p=CELLDEG %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Type), alpha = 0.6, size = 2.5) +
    scale_color_manual(
      name = "Expression",
      values = c(ns = "grey70", up = "#E41A1C", down = "#377EB8")
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
    geom_label_repel(
      data = top_genes,
      aes(label = Gene_Symbol, color = Type),
      size = 3.8,
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
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.1)),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  print(p)
  dev.off()
  
  sig_gene_symbols <- CELLDEG_sig %>% filter(Type == "up") %>% pull(Gene_Symbol)
  
  if (length(sig_gene_symbols) > 0) {
    sig_entrez <- bitr(sig_gene_symbols,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
    
    if (!is.null(sig_entrez)) {
      go_enrich <- enrichGO(gene = sig_entrez$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            ont = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1)
      kegg_enrich <- enrichKEGG(gene = sig_entrez$ENTREZID,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.1)
      write.csv(as.data.frame(go_enrich), 
                file = paste0("MSCs_GO_enrichment_", comp_name, ".csv"), 
                row.names = FALSE)
      write.csv(as.data.frame(kegg_enrich), 
                file = paste0("MSCs_KEGG_enrichment_", comp_name, ".csv"), 
                row.names = FALSE)
      
    }
  }
}

kegg_plot <- function(up_kegg, down_kegg){

  dat <- rbind(up_kegg, down_kegg)

  dat$GeneRatio <- sapply(strsplit(dat$GeneRatio, "/"),
                          function(x) as.numeric(x[1]) / as.numeric(x[2]))

  dat$GeneRatio_signed <- ifelse(
    dat$diseasestate == "EcO",
    dat$GeneRatio,
    -dat$GeneRatio
  )

  dat$logp <- -log10(dat$pvalue)

  pathway_order <- dat %>%
    group_by(Description) %>%
    summarise(min_p = min(pvalue, na.rm = TRUE)) %>%
    arrange(min_p) %>%           # p 越小越靠前
    pull(Description)
  
  dat$Description <- factor(dat$Description, levels = pathway_order)
  
  ## 6. 分组
  eco_dat <- dat %>% filter(diseasestate == "EcO")
  eue_dat <- dat %>% filter(diseasestate == "EuE")
  
  p <- ggplot() +
  geom_col(
    data = eco_dat,
    aes(x = GeneRatio_signed,
        y = Description,
        fill = logp),
    width = 0.6
  ) +
    scale_fill_gradient(
      name = "EcO (-log10 p)",
      low = "#80b1d3",
      high = "#2c5d73",
      guide = guide_colorbar(order = 2)
    ) +
    
    new_scale_fill() +

  geom_col(
    data = eue_dat,
    aes(x = GeneRatio_signed,
        y = Description,
        fill = logp),
    width = 0.6
  ) +
    scale_fill_gradient(
      name = "EuE (-log10 p)",
      low = "#E4C755",
      high = "#FF8F00",
      guide = guide_colorbar(order = 1)
    ) +
    
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(
      labels = abs,
      name = "Gene ratio"
    ) +
    
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 13),
      legend.position = "right"
    )
  
  print(p)
}

library(ggplot2)
library(dplyr)
library(ggnewscale)  


go_plot <- function(up_go, down_go){
  
  dat <- rbind(up_go, down_go)
  
  dat$GeneRatio <- sapply(strsplit(dat$GeneRatio, "/"),
                          function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  dat$logp <- -log10(dat$pvalue)
  
  dat$Description <- factor(dat$Description, levels = unique(dat$Description))
  
  eco_dat <- dat %>% dplyr::filter(diseasestate == "EcO")
  eue_dat <- dat %>% dplyr::filter(diseasestate == "EuE")
  
  p <- ggplot() +
    
    geom_point(data = eco_dat,
               aes(x = diseasestate,
                   y = Description,
                   size = GeneRatio,
                   color = logp),
               alpha = 0.9) +
    scale_color_gradient(
      name = "EcO  (-log10 p)",
      low = "#80b1d3",
      high = "#2c5d73",
      aesthetics = "color",
      guide = guide_colorbar(order = 2)
    ) +
    
    new_scale_color() +  
    geom_point(data = eue_dat,
               aes(x = diseasestate,
                   y = Description,
                   size = GeneRatio,
                   color = logp),
               alpha = 0.9) +
    scale_color_gradient(
      name = "EuE  (-log10 p)",
      low = "#E4C755",
      high = "#FF8F00",
      aesthetics = "color",
      guide = guide_colorbar(order = 1)
    ) +
    
    scale_size_continuous(range = c(3, 8), name = "GeneRatio") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 11),
          legend.position = "right")
  
  print(p)
}

up_kegg   <- read_csv("MSC_NR4A3_KEGG_EuE_vs_EcO.csv")   %>% 
  mutate(group =  1)  %>% 
  arrange(pvalue) %>% 
  slice_head(n = 10)
down_kegg <- read_csv("MSC_NR4A3_KEGG_EcO_vs_EuE.csv")   %>% 
  mutate(group = -1)%>% 
  arrange(pvalue) %>% 
  slice_head(n = 10)

up_go   <- read_csv("MSC_NR4A3_GO_EuE_vs_EcO.csv")   %>% 
  mutate(group =  1)  %>% 
  mutate(diseasestate = "EcO")%>%
  arrange(pvalue) %>% 
  slice_head(n = 10)
down_go <- read_csv("MSC_NR4A3_GO_EcO_vs_EuE.csv")   %>% 
  mutate(group = -1)%>% 
  mutate(diseasestate = "EuE")%>%
  arrange(pvalue) %>% 
  slice_head(n = 10)

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = 'kegg_up_down.png')

g_go=go_plot(up_go,down_go)
print(g_go)
ggsave(g_go,filename = 'go_up_down.png')

#Fig.2g
gene_sets <- list(
  IRON_TRANSMEMBRANE_TRANSPORT=c("ABCB6","ABCB7","ABCC5","ATP7A","HAMP","IFNG","ISCU","MCOLN1","MCOLN2","MIR210","MMGT1","SCARA5","SLC11A1","SLC11A2","SLC25A28","SLC25A37","SLC39A14","SLC39A8","SLC40A1","SLC48A1","STEAP2","TTYH1"),
  RESPONSE_TO_IRON=c("ACO1","ATG5","BECN1","MAP1LC3A","PDX1","SNCA"),
  TREATING_IRON_OVERLOAD=c("CP","HAMP","HPX","LCN2","MT2A","SLC40A1","TFRC"),
  FERROPTOSIS=c("ACSL1","ACSL3","ACSL4","ACSL5","ACSL6","AIFM2","AKR1C1","AKR1C2","AKR1C3","ALOX15","ATG5","ATG7","BACH1","CBS","CHMP5","CHMP6","CISD1","COQ2","CP","CTH","CYBB","DPP4","FDFT1","FTH1","FTL","FTMT","GCH1","GCLC","GCLM","GPX4","GSS","HMGCR","HMOX1","HSPB1","IREB2","LPCAT3","MAP1LC3A","MAP1LC3B","MAP1LC3C","NCOA4","NOX1","NOX4","PCBP1","PCBP2","PHKG2","POR","PRNP","SAT1","SAT2","SLC11A2","SLC1A5","SLC38A1","SLC39A14","SLC39A8","SLC3A2","SLC40A1","SLC7A11","STEAP3","TF","TFRC","TP53","TXNRD1","VDAC2","VDAC3")
)
scRNA2<-subset(scRNA,cell_type=="MSC_NR4A3")
for (set_name in names(gene_sets)) {
  genes <- gene_sets[[set_name]]
  scRNA2 <- AddModuleScore(scRNA2, features = list(genes), name = set_name)
}
score_data <- FetchData(scRNA2, 
                        vars = c("diseasestate", paste0(names(gene_sets), "1")))

comparisons <- list(c("EcO", "EuE"),c("EcO", "Ctrl"), c("EcOA", "Ctrl"), c("Ctrl", "CA"))

score_data <- score_data[!is.na(score_data$diseasestate), ]

for (set_name in names(gene_sets)) {
  score_col <- paste0(set_name, "1")
  
  p <- ggplot(score_data, aes(x = diseasestate, y = .data[[score_col]], fill = diseasestate)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                       label = "p.signif", tip.length = 0.01, size = 2) +
    theme_classic(base_size = 14) +
    labs(title = set_name),
         y = "Module Score", x = "") +
    scale_fill_brewer(palette = "Set2") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold",size = 8),
      legend.position = "none"
    )

  ggsave(filename = paste0("Violin_", set_name, "_MSCs_by_diseasestate.png"),
         plot = p, width = 4, height = 3)
}


#Fig.2h-i,Fig.S2d
library(CellChat)
immune <- readRDS("immune.rds")
scRNA3<-merge(scRNA2,EPI,add.cell.ids = c("MSC","epi"))
scRNA3<-merge(scRNA2,y=c(EPI,immune),add.cell.ids = c("MSC","epi","immune"))
scRNA3[["RNA"]] <- JoinLayers(scRNA3[["RNA"]])
colnames(scRNA3@assays[["RNA"]]@layers[["data"]]) <- rownames(scRNA3@meta.data)

Idents(scRNA3) = "Diseasestate"
scRNA.EuE = scRNA3[,Idents(scRNA3) %in% c("EuE")]
scRNA.EcO = scRNA3[,Idents(scRNA3) %in% c("EcO")]
scRNA.EcOA = scRNA3[,Idents(scRNA3) %in% c("EcOA")]
scRNA.Ctrl = scRNA3[,Idents(scRNA3) %in% c("Ctrl")]

seurat.list <- list(
  "Ctrl"  = scRNA.Ctrl,
  "EcO"   = scRNA.EcO,
  "EcOA"   = scRNA.EcOA,
  "EuE"   = scRNA.EuE
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
  CC.input <- GetAssayData(seurat.subset, assay = "RNA", layer = "data")
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
                "EuE"   = CellChat_object_EuE,
                "EcO"   = CellChat_object_EcO,
                "EcOA"   = CellChat_object_EcOA
)
CC.merged <- mergeCellChat(CC.list[c("Ctrl","EuE","EcO","EcOA")],
                           add.names = c("Ctrl","EuE","EcO","EcOA"))
saveRDS(CC.merged, "CellChat_Merged_object.rds")

All_Celltype_2D = Celltype_2D_plot(CC.plot.list = CC.list, x.name = "All")

comparisons <- list(
  "Ctrl_vs_EcO"  = c(1,3),
  "Ctrl_vs_EcOA"  = c(1,4),
  "EuE_vs_EcO"   = c(2,3)
)

# 循环绘制 circle 和 heatmap
for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  message("Running comparison: ", name, " -> ", paste(comp, collapse = "_"))
  # circle plot
  pdf(paste0("compare_circle_", name, ".pdf"), width = 4, height = 4.5)
  par(xpd=TRUE)
  netVisual_diffInteraction(CC.merged, weight.scale = TRUE, comparison = comp,
                            color.use=c("#990000","#97cac1","#8c96c6", "#3860a0"))
  netVisual_diffInteraction(CC.merged, weight.scale = TRUE, measure = "weight",
                            color.use=c("#990000","#97cac1","#8c96c6", "#3860a0"), comparison = comp)
  dev.off()
  # rankNet plot
  pdf(paste0("compare_ranknet_", name, ".pdf"), width = 4.5, height = 5)
  g1 <- rankNet(CC.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE,
                comparison = comp,color.use = c("#E4C755","#80b1d3"),
                sources.use = "MSC_NR4A3",targets.use = "Epithelial",
                font.size = 12)
  g2 <- rankNet(CC.merged, mode = "comparison", stacked = FALSE, do.stat = TRUE,
                comparison = comp,color.use = c("#E4C755","#80b1d3"),
                sources.use = "MSC_NR4A3",targets.use = "Epithelial",
                font.size = 12)
  print(g1)
  print(g2)
  dev.off()
}

#Fig.2j
df <- FetchData(scRNA2, vars = c("MIF", "diseasestate"))

p15 <- ggplot(df, aes(x = diseasestate, y = MIF, fill = diseasestate)) +
  geom_violin(color = NA, width = 1,trim = F) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("EuE", "EcO"),
      c("Ctrl", "EcOA"),
      c("Ctrl", "CA")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "MIF Expression in MSC_NR4A3", x = "", y = "MIF Expression") +
  scale_fill_manual(values = c("#8dd3c7","#E4C755","#80b1d3","#bebada","#fb8072")) +
  theme_classic(base_size = 12) +  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("MIF.png",p15,width = 4,height=4)

EPI<-readRDS("Epithelial.rds")
df <- FetchData(EPI, vars = c("CD74", "diseasestate"))

p15 <- ggplot(df, aes(x = diseasestate, y = CD74, fill = diseasestate)) +
  geom_violin(color = NA, width = 1,trim = F) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("EuE", "EcO"),
      c("Ctrl", "EcOA"),
      c("Ctrl", "CA")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "CD74 Expression in Epi", x = "", y = "CD74 Expression") +
  scale_fill_manual(values = c("#8dd3c7","#E4C755","#80b1d3","#bebada","#fb8072")) +
  theme_classic(base_size = 12) +  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("CD74.png",p15,width = 4,height=4)

#Fig.2k
MSC<-readRDS("GEO/MSC.rds")
df <- FetchData(MSC, vars = c("MIF", "diseasestate"))

p15 <- ggplot(df, aes(x = diseasestate, y = MIF, fill = diseasestate)) +
  geom_violin(color = NA, width = 1,trim = F) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("EuE", "EcO")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "MIF Expression in MSC_NR4A3", x = "", y = "MIF Expression") +
  scale_fill_manual(values = c("#E4C755","#80b1d3")) +
  theme_classic(base_size = 12) +  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("GEO/MIF.png",p15,width = 4,height=4)

EPI<-readRDS("GEO/Epithelial.rds")
df <- FetchData(EPI, vars = c("CD74", "diseasestate"))

p15 <- ggplot(df, aes(x = diseasestate, y = CD74, fill = diseasestate)) +
  geom_violin(color = NA, width = 1,trim = F) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "white", color = "black") +  # 箱线图更窄
  stat_compare_means(
    comparisons = list(
      c("EuE", "EcO")
    ),
    method = "wilcox.test",
    label ="p.format" ,
    tip.length = 0.005,
    step.increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, NA))+
  labs(title = "CD74 Expression in Epi", x = "", y = "CD74 Expression") +
  scale_fill_manual(values = c("#E4C755","#80b1d3")) +
  theme_classic(base_size = 12) +  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = ""
  )
ggsave("GEO/CD74.png",p15,width = 4,height=4)
