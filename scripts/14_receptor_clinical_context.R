# Literature-verified clinical/druggability context for the 4 candidate
# receptors, going beyond DGIdb's raw interaction counts (09/11). Compiled
# from live web search 2026-07-10 (sources noted per row); this is manual
# curation, not a database query, matching the approach already used for the
# miRNA compounds in 09_regulator_druggability.R.

source("scripts/config.R")
suppressMessages({
  library(dplyr)
  library(readr)
})

receptor_clinical_context <- tribble(
  ~gene,         ~compound,                         ~modality,                    ~indication_status,
  "ITGA1/ITGB1", "SAN-300",                          "anti-alpha1(VLA-1) antibody", "Phase 1 completed in rheumatoid arthritis (NCT02047604) -- NOT tested in fibrosis; no fibrosis trial identified for this specific alpha1beta1-targeting compound",
  "ITGA1/ITGB1", "(literature target, not our hit)", "ITGA11 (alpha11beta1) -- different alpha subunit", "Strongest stellate-cell-specific fibrosis literature is for ITGA11, not ITGA1 (hedgehog-ITGA11 axis, miR-12135/ITGA11 axis, 2023-2024) -- in THIS dataset ITGA11 is barely detected in Mesenchyme (5.2% of cirrhotic cells) vs ITGA1 (41.4%), a real discrepancy worth flagging rather than assuming the literature target applies here",
  "CD44",        "RG7356",                           "anti-CD44 humanized antibody", "Phase 1 in CD44+ solid tumors -- terminated early, no dose-response relationship observed; not a fibrosis compound and had its own developmental problems",
  "CD44",        "(related target)",                 "CD248 antibody-drug conjugate", "A different myofibroblast-specific receptor (CD248, not CD44) has an ADC in preclinical development for liver fibrosis -- an analogous approach, not a CD44 drug",
  "EDNRB",       "Bosentan",                          "dual ETA/ETB antagonist (approved, PAH)", "Lowers portal pressure in cirrhotic animal models and small clinical studies, but NOT approved for cirrhosis/portal hypertension; carries hepatotoxicity risk (LiverTox: enzyme elevation, rare acute liver injury) -- an important safety caveat for a liver-disease application",
  "EDNRB",       "Ambrisentan",                       "ETA-selective antagonist (approved, PAH)", "Better hepatic safety profile than bosentan in PAH populations, but is ETA-selective -- our specific hit was EDNRB (ETB), so ambrisentan would not directly engage the receptor that changed in this analysis",
  "EDNRA",       "Zibotentan",                        "ETA-selective antagonist (investigational)", "Actively being tested in an ongoing Phase 2b cirrhosis trial (ZEAL-UNLOCK, NCT06269484, AstraZeneca) combined with dapagliflozin -- but targets ETA (EDNRA), which was ALREADY significant in healthy Mesenchyme in this data, not the EDNRB signal that is new in cirrhosis"
)

write_csv(receptor_clinical_context, file.path(TABLE_DIR, "14_receptor_clinical_context.csv"))

cat("Receptor clinical/druggability context (literature-verified 2026-07-10):\n")
print(as.data.frame(receptor_clinical_context))

cat("\nKey honest takeaway: none of the 4 receptors has a compound that is both\n",
    "(a) currently in active clinical development AND (b) precisely matched to\n",
    "the specific gene/receptor subtype that changed in this analysis.\n",
    "ITGA1 != the literature's ITGA11; EDNRB != the ETA that zibotentan/ambrisentan target.\n",
    "This is reported as a real gap rather than papered over with a loose 'a drug exists' claim.\n", sep = "")
