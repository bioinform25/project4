library(Seurat)
library(dplyr)
library(ggplot2)

source("C:/Users/SAMSUNG/Desktop/project4/scripts/config.R")

liver_combined <- readRDS(file.path(CACHE_DIR, "04_liver_final.rds"))

comp_table <- liver_combined@meta.data %>%
  count(condition, cell_type) %>%
  group_by(condition) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

write.csv(comp_table, file.path(TABLE_DIR, "05_composition_by_condition.csv"), row.names = FALSE)

# Chi-square test of independence: cell_type x condition
cont_table <- table(liver_combined$cell_type, liver_combined$condition)
chisq_res <- chisq.test(cont_table)
cat("Chi-square test, cell type composition ~ condition:\n")
print(chisq_res)

# Per cell-type proportion test (healthy vs cirrhotic), one-vs-rest
cell_types <- unique(liver_combined$cell_type)
per_type <- lapply(cell_types, function(ct) {
  in_type <- liver_combined$cell_type == ct
  tab <- table(in_type, liver_combined$condition)
  pt <- prop.test(tab)
  data.frame(cell_type = ct,
             pct_healthy = 100 * tab["TRUE", "healthy"] / sum(tab[, "healthy"]),
             pct_cirrhotic = 100 * tab["TRUE", "cirrhotic"] / sum(tab[, "cirrhotic"]),
             p_value = pt$p.value)
})
per_type_df <- do.call(rbind, per_type)
per_type_df$p_adj <- p.adjust(per_type_df$p_value, method = "BH")
per_type_df <- per_type_df[order(per_type_df$p_adj), ]
write.csv(per_type_df, file.path(TABLE_DIR, "05_composition_proportion_test.csv"), row.names = FALSE)
cat("\nPer cell-type proportion test (BH-adjusted):\n")
print(per_type_df)

p_comp <- ggplot(comp_table, aes(x = cell_type, y = pct, fill = condition)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "% of cells", x = "Cell type", title = "Cell-type composition: healthy vs cirrhotic")
ggsave(file.path(FIG_DIR, "05_composition_barplot.png"), p_comp, width = 9, height = 6, dpi = 150)
