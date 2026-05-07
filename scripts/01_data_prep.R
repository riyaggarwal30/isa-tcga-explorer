# scripts/01_data_prep.R
library(IsoformSwitchAnalyzeR)

# 1. Load raw example data
data("exampleSwitchListAnalyzed")

# 2. Save it to the data folder using a RELATIVE path
# This ensures it works on ANY computer
saveRDS(exampleSwitchListAnalyzed, "data/switchData.rds")

print("Data preparation complete! switchData.rds is now in the data/ folder.")