library(tidyverse)
library(Seurat)
library(patchwork)

rm(list = ls())
setwd("~/projects1/projMD1/analysis1/")

# 1. Load and Prep
dlist0 <- list(
  APF_48h = readRDS("data/data2_48h.rds"),
  APF_72h = readRDS("data/data2_72h.rds"),
  APF_96h = readRDS("data/data2_96h.rds")
)

# 2. Run PCA2
runPCA2 <- function(data00, type00){
  data00 <- data00 %>%
    subset(type1 == type00) %>%
    DietSeurat(layers = c("counts", "data")) %>%
    FindVariableFeatures(nfeatures = 500) %>%
    ScaleData() %>%
    RunPCA(seed.use = 123)
}

##
dlist1 <- list(
  APF_48h = runPCA2(dlist0[["APF_48h"]], "LPLC1"),
  APF_72h = runPCA2(dlist0[["APF_72h"]], "LPLC1"),
  APF_96h = runPCA2(dlist0[["APF_96h"]], "LPLC1")
)
saveRDS(dlist1, file = "data/dlist3_LPLC1.rds")

dlist2 <- list(
  APF_48h = runPCA2(dlist0[["APF_48h"]], "LPLC2"),
  APF_72h = runPCA2(dlist0[["APF_72h"]], "LPLC2"),
  APF_96h = runPCA2(dlist0[["APF_96h"]], "LPLC2")
)
saveRDS(dlist2, file = "data/dlist3_LPLC2.rds")

dlist4 <- list(
  APF_48h = runPCA2(dlist0[["APF_48h"]], "LC4"),
  APF_72h = runPCA2(dlist0[["APF_72h"]], "LC4"),
  APF_96h = runPCA2(dlist0[["APF_96h"]], "LC4")
)
saveRDS(dlist4, file = "data/dlist3_LC4.rds")

# 3. Plot PCA2: covariates on PC1-6
plotPCA1 <- function(data00){
  p1 <- DimPlot(data00, group.by = "group", dims = 1:2)
  p2 <- DimPlot(data00, group.by = "group", dims = 3:4)
  p3 <- DimPlot(data00, group.by = "group", dims = 5:6)
  p01 <- ((p1 + p2 + p3) & scale_color_manual(values = c("blue3", "red3"))) + plot_layout(guides = "collect", ncol = 3)
  
  p1 <- DimPlot(data00, group.by = "genotype", dims = 1:2)
  p2 <- DimPlot(data00, group.by = "genotype", dims = 3:4)
  p3 <- DimPlot(data00, group.by = "genotype", dims = 5:6)
  p02 <- ((p1 + p2 + p3) & scale_color_brewer(palette = "Dark2")) + plot_layout(guides = "collect", ncol = 3)
  
  p1 <- FeaturePlot(data00, dims = 1:2, features = "lncRNA:roX1", max.cutoff = "q95", min.cutoff = "q5")
  p2 <- FeaturePlot(data00, dims = 3:4, features = "lncRNA:roX1", max.cutoff = "q95", min.cutoff = "q5")
  p3 <- FeaturePlot(data00, dims = 5:6, features = "lncRNA:roX1", max.cutoff = "q95", min.cutoff = "q5")
  p03 <- (p1 + p2 + p3) + plot_layout(guides = "collect", ncol = 3)
  
  p1 <- FeaturePlot(data00, dims = 1:2, features = "nCount_RNA", max.cutoff = "q95", min.cutoff = "q5")
  p2 <- FeaturePlot(data00, dims = 3:4, features = "nCount_RNA", max.cutoff = "q95", min.cutoff = "q5")
  p3 <- FeaturePlot(data00, dims = 5:6, features = "nCount_RNA", max.cutoff = "q95", min.cutoff = "q5")
  p04 <- (p1 + p2 + p3) + plot_layout(guides = "collect", ncol = 3)
  
  theme1 <- theme(
    plot.title = element_blank(), legend.text = element_text(size = 9),
    axis.text = element_blank(), axis.title = element_text(size = 9)
  )
  p01/p02/p03/p04 & theme1
}

##
for(t00 in c("48h", "72h", "96h")){
  plotPCA1(dlist1[[paste0("APF_", t00)]])
  ggsave(filename = paste0("results/pca1_LPLC1_", t00, ".pdf"), units = "in", width = 6.5, height = 8)
  ##
  plotPCA1(dlist2[[paste0("APF_", t00)]])
  ggsave(filename = paste0("results/pca1_LPLC2_", t00, ".pdf"), units = "in", width = 6.5, height = 8)
  ##
  plotPCA1(dlist4[[paste0("APF_", t00)]])
  ggsave(filename = paste0("results/pca1_LC4_", t00, ".pdf"), units = "in", width = 6.5, height = 8)
}

# 4. Plot PCA2: heatmaps for PC1
plotPCA2 <- function(data00){
  pcells1 <- TopCells(data00, ncells = ncol(data00), dim = 1)
  pgenes1 <- TopFeatures(data00, dim = 1, nfeatures = 30, balanced = T)
  pgenes1 <- c(pgenes1[[2]], rev(pgenes1[[1]]))
  
  pdata1 <- FetchData(data00, vars = c("PC_1", "PC_2", pgenes1), layer = "scale.data")
  pdata2 <- FetchData(data00, vars = c("PC_1", "PC_2", pgenes1), layer = "data")
  
  p1 <- SingleRasterMap(pdata1[,pgenes1], cell.order = pcells1, feature.order = pgenes1, raster = F)
  p2 <- SingleRasterMap(pdata2[,pgenes1], cell.order = pcells1, feature.order = pgenes1, disp.min = 0, disp.max = log1p(20), color = BlueAndRed(), raster = F)
  
  p1 + p2 & theme(legend.position = "bottom")
}

##
plotPCA2(dlist1[["APF_48h"]])
ggsave(filename = "results/pca2_LPLC1_48h.pdf", units = "in", width = 8.5, height = 6)

plotPCA2(dlist2[["APF_48h"]])
ggsave(filename = "results/pca2_LPLC2_48h.pdf", units = "in", width = 8.5, height = 6)

# 5. Plot PCA2: genes on PC1/PC2
plotPCA3 <- function(dlist00, pgenes00){
  p1 <- FeaturePlot(dlist00[["APF_48h"]], pgenes00, max.cutoff = "q95", min.cutoff = "q5", combine = T, ncol = length(pgenes00))
  p2 <- FeaturePlot(dlist00[["APF_72h"]], pgenes00, max.cutoff = "q95", min.cutoff = "q5", combine = T, ncol = length(pgenes00))
  p3 <- FeaturePlot(dlist00[["APF_96h"]], pgenes00, max.cutoff = "q95", min.cutoff = "q5", combine = T, ncol = length(pgenes00))
  
  theme1 <- theme(
    plot.title = element_text(size = 11), legend.text = element_text(size = 9), legend.key.size = unit("10", "pt"),
    axis.text = element_blank(), axis.title = element_text(size = 9)
  )
  p1 / p2 / p3 & theme1
}

##
plotPCA3(dlist2, c("dpr13", "CG17716", "beat-VI", "SiaT"))
ggsave(filename = "results/pca3_LPLC2.pdf", units = "in", width = 10, height = 6)

plotPCA3(dlist1, c("DIP-kappa", "CG33543"))
ggsave(filename = "results/pca3_LPLC1.pdf", units = "in", width = 5, height = 6)
