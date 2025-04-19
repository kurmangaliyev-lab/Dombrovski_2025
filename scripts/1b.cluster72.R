library(tidyverse)
library(Seurat)
library(patchwork)

rm(list = ls())
setwd("~/projects1/projMD1/analysis1/")

# 1. Load and Prep
data0 <- CreateSeuratObject(Read10X("../data1/Mark4/outs/filtered_feature_bc_matrix/"), project = "APF_72h")
data0$Perc_MT <- PercentageFeatureSet(data0, "mt:")

data0 <- AddMetaData(
  data0,
  read.table("../data1/multiplex1/Mark4/demuxlet0.best", header = T, sep = "\t") %>% 
    select(cellID = BARCODE, is_sng = DROPLET.TYPE, genotype = SNG.BEST.GUESS) %>% 
    mutate(genotype = ifelse(is_sng == "SNG", genotype, is_sng)) %>%
    mutate(group = "X") %>%
    mutate(group = ifelse(genotype %in% c("line_320","line_399","line_535"), "early", group)) %>%
    mutate(group = ifelse(genotype %in% c("line_730","line_799","line_805"), "late", group)) %>%
    column_to_rownames("cellID")
)

saveRDS(data0, "data/data1_72h.rds")

# 2. Cluster
table(data0$genotype, data0$group)
VlnPlot(data0, group.by = "group", features = "nCount_RNA") + 
  scale_y_continuous(trans = "log10", breaks = c(500, 2000, 5000, 10000, 50000)) +
  geom_hline(yintercept = c(5000, 50000)) +
  NoLegend()
data0 <- subset(data0, nCount_RNA > 5000 & nCount_RNA < 50000 & group != "X")
dim(data0)

data0 <- data0 %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 1000) %>%
  ScaleData() %>%
  RunPCA()
ElbowPlot(data0, ndims = 50)

dims0 <- 1:5
data0 <- data0 %>%
  RunTSNE(dims = dims0) %>%
  FindNeighbors(dims = dims0) %>%
  FindClusters(resolution = 0.05)

data0 <- RenameIdents(data0, "0"="LPLC2", "1"="LPLC1", "2"="LC4", "3"="X")
data0$type1 <- factor(Idents(data0), levels = c("LPLC2", "LPLC1", "LC4", "X"))
Idents(data0) <- data0$type1

table(data0$type1, data0$genotype)

pgenes1 <- c("acj6", "fkh", "Sox102F", "salm", "salr", "Lim1", "br", "tup", "Lim3", "kn", "lov")
cowplot::plot_grid(
  DimPlot(data0, label = T, label.size = 4) + 
    NoLegend() + theme(axis.text = element_blank()),
  DotPlot(data0, pgenes1, scale = F) + coord_flip() +
    theme(axis.title = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
  rel_widths = c(1.6,1)
)
ggsave(filename = "results/cluster_72h.pdf", units = "in", width = 7, height = 3.5)

saveRDS(data0, file = "data/data2_72h.rds")
