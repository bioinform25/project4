library(Seurat)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

liver_combined <- readRDS(file.path(CACHE_DIR, "01_liver_combined_harmony.rds"))

# Dims chosen after inspecting 01_elbow_preharmony.png (elbow flattens well before 30;
# paper's original dims=1:11 still reasonable post-harmony)
USE_DIMS <- 1:15

liver_combined <- FindNeighbors(liver_combined, reduction = "harmony", dims = USE_DIMS)
liver_combined <- FindClusters(liver_combined, resolution = 0.6)
liver_combined <- RunUMAP(liver_combined, reduction = "harmony", dims = USE_DIMS)

p_cluster <- DimPlot(liver_combined, reduction = "umap", label = TRUE, pt.size = 0.1) + ggtitle("Clusters")
p_condition <- DimPlot(liver_combined, reduction = "umap", group.by = "condition", pt.size = 0.1) + ggtitle("Condition")
p_sample <- DimPlot(liver_combined, reduction = "umap", group.by = "sample", pt.size = 0.1) +
  ggtitle("Sample (post-Harmony mixing check)") + theme(legend.text = element_text(size = 6))

ggsave(file.path(FIG_DIR, "02_umap_clusters.png"), p_cluster, width = 7, height = 6, dpi = 150)
ggsave(file.path(FIG_DIR, "02_umap_condition.png"), p_condition, width = 7, height = 6, dpi = 150)
ggsave(file.path(FIG_DIR, "02_umap_sample.png"), p_sample, width = 8, height = 6, dpi = 150)

cat("Cluster counts:\n")
print(table(Idents(liver_combined)))

saveRDS(liver_combined, file.path(CACHE_DIR, "02_liver_clustered.rds"))
