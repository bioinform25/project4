# miRNA regulators of the Mesenchyme fibrogenic signature (LUM/THY1/THBS2).
#
# Rationale: continues the upstream-regulator question from 07_tf_regulon_
# analysis.R at the post-transcriptional level. Queries multiMiR for
# experimentally validated (miRTarBase/miRecords/TarBase) and predicted
# (TargetScan) miRNA-target interactions for the 3 genes that cleanly
# localize to activated Mesenchyme, then flags miRNAs that (a) target more
# than one of the 3 genes, and (b) belong to families with an established
# fibrosis literature (miR-29, -21, -192, -200, -214, -34a, -122, -155).

source("scripts/config.R")
suppressMessages({
  library(multiMiR)
  library(dplyr)
  library(tidyr)
  library(readr)
})

fibro_genes <- c("LUM", "THY1", "THBS2")

# ---------------------------------------------------------------------------
# validated interactions (miRTarBase / miRecords / TarBase)
# ---------------------------------------------------------------------------
validated <- get_multimir(org = "hsa", target = fibro_genes, table = "validated",
                           summary = TRUE)@data
cat("validated interactions:", nrow(validated), "\n")

validated_clean <- validated %>%
  distinct(mature_mirna_id, target_symbol, database, support_type) %>%
  group_by(mature_mirna_id, target_symbol) %>%
  summarise(databases = paste(sort(unique(database)), collapse = ";"),
            n_evidence = n(), .groups = "drop")

# ---------------------------------------------------------------------------
# predicted interactions (TargetScan only -- the most widely used conserved-
# site predictor; other predicted tables in multiMiR would add volume without
# much additional confidence for a screen like this)
# ---------------------------------------------------------------------------
predicted <- get_multimir(org = "hsa", target = fibro_genes, table = "targetscan",
                           summary = TRUE)@data
cat("predicted (TargetScan) interactions:", nrow(predicted), "\n")

predicted_clean <- predicted %>%
  distinct(mature_mirna_id, target_symbol) %>%
  mutate(databases = "targetscan", n_evidence = NA_integer_)

# ---------------------------------------------------------------------------
# miRNAs hitting >=2 of the 3 fibrogenic genes (validated evidence only --
# stronger, multi-gene "master regulator" candidates)
# ---------------------------------------------------------------------------
multi_target_validated <- validated_clean %>%
  group_by(mature_mirna_id) %>%
  summarise(n_genes_targeted = n_distinct(target_symbol),
            targets = paste(sort(unique(target_symbol)), collapse = ";"),
            databases = paste(sort(unique(unlist(strsplit(databases, ";")))), collapse = ";"),
            total_evidence = sum(n_evidence)) %>%
  filter(n_genes_targeted >= 2) %>%
  arrange(desc(n_genes_targeted), desc(total_evidence))

cat("\nmiRNAs with validated evidence for >=2 of LUM/THY1/THBS2:\n")
print(multi_target_validated)

# ---------------------------------------------------------------------------
# cross-check against established fibrosis-associated miRNA families
# (literature-curated, not database-derived)
# ---------------------------------------------------------------------------
fibrosis_mirna_families <- c("miR-29", "miR-21", "miR-192", "miR-200", "miR-141",
                              "miR-429", "miR-214", "miR-34a", "miR-122", "miR-155")
family_pattern <- paste0("hsa-", fibrosis_mirna_families, collapse = "|")
family_pattern <- gsub("miR-", "miR-", family_pattern)  # keep as-is; grepl below is family-prefix match

flag_family <- function(mirna_id) {
  hits <- fibrosis_mirna_families[sapply(fibrosis_mirna_families, function(fam) {
    grepl(paste0("^hsa-", fam, "[a-z]?(-[0-9]p)?$"), mirna_id) ||
      grepl(paste0("^hsa-", fam, "-"), mirna_id)
  })]
  if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
}

validated_clean <- validated_clean %>%
  mutate(known_fibrosis_family = sapply(mature_mirna_id, flag_family))
predicted_clean <- predicted_clean %>%
  mutate(known_fibrosis_family = sapply(mature_mirna_id, flag_family))

cat("\nvalidated interactions where the miRNA belongs to a known fibrosis-miRNA family:\n")
print(validated_clean %>% filter(!is.na(known_fibrosis_family)) %>% arrange(target_symbol))

# ---------------------------------------------------------------------------
# save
# ---------------------------------------------------------------------------
write_csv(validated_clean, file.path(TABLE_DIR, "08_mirna_validated_targets.csv"))
write_csv(predicted_clean, file.path(TABLE_DIR, "08_mirna_predicted_targetscan.csv"))
write_csv(multi_target_validated, file.path(TABLE_DIR, "08_mirna_multi_target_candidates.csv"))

cat("\nDone. Results in results/tables/08_*.csv\n")
