RAW_DATA_DIR <- "C:/Users/SAMSUNG/Desktop/5-1 1~2/GSE136103/GSE136103/"
CACHE_DIR    <- "C:/Users/SAMSUNG/Desktop/project4/data_cache/"
FIG_DIR      <- "C:/Users/SAMSUNG/Desktop/project4/results/figures/"
TABLE_DIR    <- "C:/Users/SAMSUNG/Desktop/project4/results/tables/"

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

CANDIDATE_GENES <- c("CCL21", "CXCL8", "CCL20", "EPCAM", "LUM", "THY1", "THBS2")
