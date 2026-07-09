library(Seurat)
library(dplyr)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

all_folders <- list.dirs(RAW_DATA_DIR, full.names = FALSE, recursive = FALSE)
all_folders <- all_folders[all_folders != ""]
liver_folders <- all_folders[!grepl("mouse|blood|Rproj", all_folders)]
cat("Human liver sample folders:", length(liver_folders), "\n")

parse_meta <- function(folder) {
  condition <- ifelse(grepl("cirrhotic", folder), "cirrhotic", "healthy")
  num <- regmatches(folder, regexpr("(healthy|cirrhotic)([0-9]+)", folder))
  num <- regmatches(num, regexpr("[0-9]+", num))
  if (length(num) == 0) num <- "1"  # GSM4041161_cirrhotic_cd45+ has no number -> patient 1
  patient <- paste0(condition, num)
  fraction <- ifelse(grepl("cd45\\+", folder), "CD45+", "CD45-")
  list(condition = condition, patient = patient, fraction = fraction)
}

seurat_list <- list()
qc_summary <- data.frame()

for (folder in liver_folders) {
  folder_path <- file.path(RAW_DATA_DIR, folder)
  meta <- parse_meta(folder)

  counts <- Read10X(data.dir = folder_path)
  seobj <- CreateSeuratObject(counts = counts, project = folder, min.cells = 3, min.features = 300)
  seobj[["percent.mt"]] <- PercentageFeatureSet(seobj, pattern = "^MT-")
  seobj$condition <- meta$condition
  seobj$patient   <- meta$patient
  seobj$fraction  <- meta$fraction
  seobj$sample    <- folder

  n_before <- ncol(seobj)
  seobj_filtered <- subset(seobj, subset = nFeature_RNA > 300 & percent.mt < 30)
  n_after <- ncol(seobj_filtered)

  qc_summary <- rbind(qc_summary, data.frame(
    sample = folder, condition = meta$condition, patient = meta$patient,
    fraction = meta$fraction, n_cells_before_qc = n_before, n_cells_after_qc = n_after
  ))

  seurat_list[[folder]] <- seobj_filtered
  cat(sprintf("%s: %d -> %d cells\n", folder, n_before, n_after))
}

write.csv(qc_summary, file.path(TABLE_DIR, "00_qc_summary.csv"), row.names = FALSE)
cat("\nTotal cells after QC:", sum(qc_summary$n_cells_after_qc), "\n")
cat("Any samples dropped to zero cells:", any(qc_summary$n_cells_after_qc == 0), "\n")

# QC violin plots (before filtering, per-sample) for sanity check
merged_unfiltered <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
                           add.cell.ids = names(seurat_list))
p_qc <- VlnPlot(merged_unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "sample", pt.size = 0, ncol = 1) &
  theme(axis.text.x = element_text(size = 6, angle = 90))
ggsave(file.path(FIG_DIR, "00_qc_violin_postfilter.png"), p_qc, width = 14, height = 12, dpi = 150)

saveRDS(seurat_list, file.path(CACHE_DIR, "00_qc_seurat_list.rds"))
cat("\nSaved:", file.path(CACHE_DIR, "00_qc_seurat_list.rds"), "\n")
