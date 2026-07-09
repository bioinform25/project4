library(Seurat)
library(readxl)
library(dplyr)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

liver_combined <- readRDS(file.path(CACHE_DIR, "02_liver_clustered.rds"))

# Parse the paper's own supplementary lineage signature table (Ramachandran et al. 2019,
# Nature, MOESM3) instead of hand-retyping marker genes.
marker_path <- file.path(RAW_DATA_DIR, "cell_lineage_signature_genes.xlsx")
raw <- read_excel(marker_path, col_names = FALSE)
raw <- as.data.frame(raw)

lineage_markers <- list()
for (i in seq_len(nrow(raw))) {
  lineage <- raw[i, 1]
  genes <- unlist(raw[i, 3:ncol(raw)])
  genes <- genes[!is.na(genes)]
  genes <- gsub("HLA\\+AC0-DRA", "HLA-DRA", genes)  # xlsx encoding artifact for HLA-DRA
  lineage_markers[[lineage]] <- unique(genes)
}
cat("Parsed lineages:\n")
print(lineage_markers)

# Keep only genes present in the object
present <- rownames(liver_combined)
lineage_markers <- lapply(lineage_markers, function(g) intersect(g, present))
missing <- lapply(lineage_markers, function(g) g)
for (ln in names(lineage_markers)) {
  if (length(lineage_markers[[ln]]) == 0) cat("WARNING: no markers found in data for", ln, "\n")
}

# Module score per lineage
for (ln in names(lineage_markers)) {
  if (length(lineage_markers[[ln]]) >= 1) {
    liver_combined <- AddModuleScore(liver_combined, features = list(lineage_markers[[ln]]),
                                      name = paste0("score_", make.names(ln)))
  }
}

score_cols <- grep("^score_", colnames(liver_combined@meta.data), value = TRUE)
avg_scores <- liver_combined@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean)) %>%
  as.data.frame()

rownames(avg_scores) <- avg_scores$seurat_clusters
score_matrix <- as.matrix(avg_scores[, score_cols])
colnames(score_matrix) <- gsub("^score_", "", colnames(score_matrix))
colnames(score_matrix) <- sub("1$", "", colnames(score_matrix))
colnames(score_matrix) <- gsub("\\.", " ", colnames(score_matrix))

cluster_call <- colnames(score_matrix)[apply(score_matrix, 1, which.max)]
names(cluster_call) <- rownames(score_matrix)

cat("\nPer-cluster average lineage scores:\n")
print(round(score_matrix, 3))
cat("\nCluster -> assigned cell type:\n")
print(cluster_call)

annotation_table <- data.frame(cluster = names(cluster_call), cell_type = cluster_call)
write.csv(annotation_table, file.path(TABLE_DIR, "03_cluster_annotation.csv"), row.names = FALSE)
write.csv(round(score_matrix, 4), file.path(TABLE_DIR, "03_cluster_lineage_scores.csv"))

liver_combined$cell_type <- unname(cluster_call[as.character(liver_combined$seurat_clusters)])

# Marker dotplot for visual cross-check against the assigned labels
all_markers_flat <- unique(unlist(lineage_markers))
p_dot_cluster <- DotPlot(liver_combined, features = all_markers_flat, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  scale_color_gradient(low = "lightgrey", high = "red")
ggsave(file.path(FIG_DIR, "03_marker_dotplot_by_cluster.png"), p_dot_cluster, width = 16, height = 7, dpi = 150)

p_umap_annot <- DimPlot(liver_combined, reduction = "umap", group.by = "cell_type", label = TRUE, pt.size = 0.1) +
  ggtitle("Annotated cell types")
ggsave(file.path(FIG_DIR, "03_umap_annotated.png"), p_umap_annot, width = 8, height = 6, dpi = 150)

saveRDS(liver_combined, file.path(CACHE_DIR, "03_liver_annotated.rds"))
cat("\nSaved annotated object. Cell type counts:\n")
print(table(liver_combined$cell_type))
