# Master script to run all preprocessing files in order
# Run from the project root

# Output folder for the processed IS2 data
fs::dir_create(fs::path("data", "is2"))

# Preprocessing clinical data
source(fs::path("analysis", "preprocessing", "preprocessing_is2_clinical.R"))

# Preprocessing gene expression data
source(fs::path("analysis", "preprocessing", "preprocessing_is2_expression.R"))

# Preprocessing immune response data
source(fs::path("analysis", "preprocessing", "preprocessing_is2_immuneresponse.R"))

# Merging dataframes
source(fs::path("analysis", "preprocessing", "preprocessing_is2_merging.R"))

# Preprocessing geneset data
source(fs::path("analysis", "preprocessing", "preprocessing_BG3M.R"))
