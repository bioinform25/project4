# Compare healthy vs cirrhotic CellChat networks, focused on signaling INTO
# the Mesenchyme (activated HSC) population that carries the LUM/THY1/THBS2
# fibrogenic signature (04_candidate_gene_localization.R), and cross-check
# whether LUM/THY1/THBS2 themselves act as ligands/receptors in CellChatDB.

source("scripts/config.R")
suppressMessages({
  library(CellChat)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(httr)
  library(jsonlite)
  library(patchwork)
})

cc_healthy <- readRDS(file.path(CACHE_DIR, "cellchat_healthy.rds"))
cc_cirrhotic <- readRDS(file.path(CACHE_DIR, "cellchat_cirrhotic.rds"))

cat("healthy: ", nrow(cc_healthy@net$count), "cell types,",
    sum(cc_healthy@net$count), "total inferred interactions (edges)\n")
cat("cirrhotic: ", nrow(cc_cirrhotic@net$count), "cell types,",
    sum(cc_cirrhotic@net$count), "total inferred interactions (edges)\n")

# ---------------------------------------------------------------------------
# do LUM / THY1 / THBS2 appear in CellChatDB as ligands or receptors at all?
# ---------------------------------------------------------------------------
fibro_genes <- c("LUM", "THY1", "THBS2")
fibro_pattern <- paste0("(^|_)(", paste(fibro_genes, collapse = "|"), ")($|_)")
db_interactions <- CellChatDB.human$interaction
hits <- db_interactions %>%
  filter(grepl(fibro_pattern, ligand) | grepl(fibro_pattern, receptor) |
           grepl(fibro_pattern, interaction_name_2))
cat("\nCellChatDB interactions involving LUM/THY1/THBS2 as ligand or receptor:", nrow(hits), "\n")
print(hits %>% select(interaction_name, pathway_name, ligand, receptor, annotation))
write_csv(hits %>% select(interaction_name, pathway_name, ligand, receptor, annotation),
          file.path(TABLE_DIR, "10_candidate_genes_in_cellchatdb.csv"))

# ---------------------------------------------------------------------------
# merged comparison object
# ---------------------------------------------------------------------------
cc_list <- list(healthy = cc_healthy, cirrhotic = cc_cirrhotic)
cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))

png(file.path(FIG_DIR, "10_interaction_count_strength.png"), width = 2200, height = 1400, res = 300)
p1 <- compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2), measure = "count")
p2 <- compareInteractions(cc_merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
print(p1 + p2)
dev.off()

# differential pathway "information flow" ranking
png(file.path(FIG_DIR, "11_pathway_information_flow.png"), width = 2000, height = 3200, res = 300)
print(rankNet(cc_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE))
dev.off()

# ---------------------------------------------------------------------------
# signaling INTO Mesenchyme, healthy vs cirrhotic -- the core question
# ---------------------------------------------------------------------------
mes_healthy <- subsetCommunication(cc_healthy, targets.use = "Mesenchyme")
mes_cirrhotic <- subsetCommunication(cc_cirrhotic, targets.use = "Mesenchyme")

mes_healthy$condition <- "healthy"
mes_cirrhotic$condition <- "cirrhotic"
mes_all <- bind_rows(mes_healthy, mes_cirrhotic)
write_csv(mes_all, file.path(TABLE_DIR, "11_mesenchyme_incoming_signaling_all.csv"))

cat("\nsignificant L-R pairs targeting Mesenchyme: healthy =", nrow(mes_healthy),
    " cirrhotic =", nrow(mes_cirrhotic), "\n")

# pairs (source, ligand-receptor) present in cirrhotic but absent in healthy
gained <- anti_join(
  mes_cirrhotic %>% distinct(source, interaction_name, pathway_name, ligand, receptor),
  mes_healthy %>% distinct(source, interaction_name, pathway_name, ligand, receptor),
  by = c("source", "interaction_name")
) %>% arrange(source, pathway_name)
cat("\nL-R pairs into Mesenchyme gained in cirrhotic (absent in healthy):", nrow(gained), "\n")
write_csv(gained, file.path(TABLE_DIR, "11_mesenchyme_incoming_gained_in_cirrhotic.csv"))

lost <- anti_join(
  mes_healthy %>% distinct(source, interaction_name, pathway_name, ligand, receptor),
  mes_cirrhotic %>% distinct(source, interaction_name, pathway_name, ligand, receptor),
  by = c("source", "interaction_name")
) %>% arrange(source, pathway_name)
cat("L-R pairs into Mesenchyme lost in cirrhotic (present in healthy only):", nrow(lost), "\n")
write_csv(lost, file.path(TABLE_DIR, "11_mesenchyme_incoming_lost_in_cirrhotic.csv"))

# bubble plot: incoming signaling to Mesenchyme, both conditions side by side
png(file.path(FIG_DIR, "11_mesenchyme_incoming_bubble.png"), width = 2600, height = 3200, res = 300)
print(netVisual_bubble(cc_merged, targets.use = "Mesenchyme", comparison = c(1, 2),
                        angle.x = 45, font.size = 8))
dev.off()

# ---------------------------------------------------------------------------
# druggability of the top receptor genes driving the *gained* signaling
# ---------------------------------------------------------------------------
top_gained_pathways <- gained %>% count(pathway_name, sort = TRUE) %>% head(15)
cat("\ntop pathways among gained Mesenchyme-incoming signaling:\n")
print(top_gained_pathways)

receptor_genes <- unique(unlist(strsplit(gained$receptor, "_")))
receptor_genes <- receptor_genes[!is.na(receptor_genes) & receptor_genes != ""]
cat("\nunique receptor genes in the 'gained in cirrhotic' set:", length(receptor_genes), "\n")

gql_query <- function(gene_names) {
  q <- list(
    query = "query genes($names: [String!]) { genes(names: $names) { nodes { name interactions { drug { name } } } } }",
    variables = list(names = as.list(gene_names))
  )
  r <- POST("https://dgidb.org/api/graphql",
            body = toJSON(q, auto_unbox = TRUE),
            content_type("application/json"))
  stopifnot(status_code(r) == 200)
  fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
}

if (length(receptor_genes) > 0) {
  dgidb_res <- gql_query(receptor_genes)
  gene_nodes <- dgidb_res$data$genes$nodes
  receptor_drugs <- purrr::map_dfr(gene_nodes, function(n) {
    drugs <- unique(purrr::map_chr(n$interactions, ~ .x$drug$name))
    tibble(receptor_gene = n$name, n_drugs = length(drugs),
           example_drugs = paste(head(drugs, 5), collapse = "; "))
  })
  # how many distinct gained pathways/sources does each receptor gene appear in?
  receptor_freq <- gained %>%
    separate_rows(receptor, sep = "_") %>%
    count(receptor, name = "n_gained_edges") %>%
    rename(receptor_gene = receptor)

  receptor_summary <- receptor_freq %>%
    left_join(receptor_drugs, by = "receptor_gene") %>%
    mutate(n_drugs = ifelse(is.na(n_drugs), 0, n_drugs),
           example_drugs = ifelse(is.na(example_drugs), "", example_drugs)) %>%
    arrange(desc(n_gained_edges))

  write_csv(receptor_summary, file.path(TABLE_DIR, "11_gained_receptor_druggability.csv"))
  cat("\nreceptor genes driving newly-gained cirrhotic Mesenchyme signaling, by druggability:\n")
  print(receptor_summary %>% head(20))
}

cat("\nDone. results/tables/10-11_*.csv, results/figures/10-11_*.png\n")
