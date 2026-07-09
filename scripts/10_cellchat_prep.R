# CellChat ligand-receptor cell-cell communication analysis, healthy vs
# cirrhotic liver (GSE136103), deferred from the TF/miRNA regulator analysis
# session and now run as its own extension.
#
# Rationale: 06_gutliver_lps_pathway.R manually checked one specific pathway
# (LPS/TLR4 signaling into Mesenchyme). This script runs an unbiased
# ligand-receptor screen across ALL annotated cell types x all of
# CellChatDB.human (Secreted Signaling + ECM-Receptor + Cell-Cell Contact --
# ECM-Receptor especially relevant given LUM/THY1/THBS2 are ECM-associated),
# separately for healthy and cirrhotic liver, so 11_cellchat_compare.R can ask
# which signaling INTO the Mesenchyme (activated HSC) population changes with
# disease, complementing the single-pathway LPS check and the upstream TF/
# miRNA regulator analysis (07-09) with a third, cell-cell-communication axis.

source("scripts/config.R")
suppressMessages({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(dplyr)
})

obj <- readRDS(file.path(CACHE_DIR, "04_liver_final.rds"))
cat("cell_type x condition counts:\n")
print(table(obj$cell_type, obj$condition))

run_cellchat <- function(so, label) {
  cat("\n=== building CellChat object:", label, "===\n")
  data.input <- GetAssayData(so, assay = "RNA", layer = "data")
  meta <- so@meta.data
  meta$cell_type <- droplevels(as.factor(meta$cell_type))

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
  cellchat@DB <- CellChatDB.human  # all 3 categories: Secreted Signaling, ECM-Receptor, Cell-Cell Contact

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # population.size = TRUE: this cohort mixes CD45+/CD45- sorted fractions, so
  # cell-type proportions in the data do NOT reflect true tissue abundance --
  # correcting for group size is more appropriate than assuming equal-sized
  # populations (the default population.size = FALSE)
  cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  cat(label, "- significant pathways:", length(cellchat@netP$pathways), "\n")
  cellchat
}

seurat_healthy <- subset(obj, condition == "healthy")
seurat_cirrhotic <- subset(obj, condition == "cirrhotic")
cat("\nhealthy cells:", ncol(seurat_healthy), " cirrhotic cells:", ncol(seurat_cirrhotic), "\n")

cc_healthy <- run_cellchat(seurat_healthy, "healthy")
saveRDS(cc_healthy, file.path(CACHE_DIR, "cellchat_healthy.rds"))

cc_cirrhotic <- run_cellchat(seurat_cirrhotic, "cirrhotic")
saveRDS(cc_cirrhotic, file.path(CACHE_DIR, "cellchat_cirrhotic.rds"))

cat("\nDone. cellchat_healthy.rds / cellchat_cirrhotic.rds saved to data_cache/\n")
