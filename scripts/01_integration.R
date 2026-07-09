library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

seurat_list <- readRDS(file.path(CACHE_DIR, "00_qc_seurat_list.rds"))

liver_combined <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                         add.cell.ids = names(seurat_list))
liver_combined <- JoinLayers(liver_combined)

liver_combined <- NormalizeData(liver_combined, normalization.method = "LogNormalize", scale.factor = 10000)
liver_combined <- FindVariableFeatures(liver_combined, selection.method = "vst", nfeatures = 2000)
liver_combined <- ScaleData(liver_combined)
liver_combined <- RunPCA(liver_combined, features = VariableFeatures(liver_combined), npcs = 30)

p_elbow <- ElbowPlot(liver_combined, ndims = 30)
ggsave(file.path(FIG_DIR, "01_elbow_preharmony.png"), p_elbow, width = 6, height = 5, dpi = 150)

# Harmony integration, batch = individual 10x sample (each folder is a separate run/chip)
liver_combined <- RunHarmony(liver_combined, group.by.vars = "sample", dims.use = 1:30)

saveRDS(liver_combined, file.path(CACHE_DIR, "01_liver_combined_harmony.rds"))
cat("Saved harmony-integrated object. Cells:", ncol(liver_combined), "Genes:", nrow(liver_combined), "\n")
cat("Inspect", file.path(FIG_DIR, "01_elbow_preharmony.png"), "to pick dims for 02_clustering.R\n")
