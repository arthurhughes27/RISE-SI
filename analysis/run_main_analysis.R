# Master file: runs the main analysis (gene expression at day 7) from the raw data, in order.
# Each script clears the environment when it ends. Run from the project root.
# The screen (01) is slow on the first run; its surrogate indices are saved and reused afterwards.

# ---- Processed data (data/) from the raw data (data-raw/) --------------------------------------
source(fs::path("analysis", "preprocessing", "make_analysis_data.R"))

# ---- Data description -----------------------------------------------------------------------------
source(fs::path("analysis", "descriptive", "descriptive.R"))

# ---- Main analysis -----------------------------------------------------------------------------
source(fs::path("analysis", "application", "main", "01_screen.R"))              # screen of all gene sets
source(fs::path("analysis", "application", "main", "02_m156_aggregation.R"))    # top gene set M156.0
source(fs::path("analysis", "application", "main", "03_m156_single_genes.R"))
source(fs::path("analysis", "application", "main", "04_m156_learners.R"))
source(fs::path("analysis", "application", "main", "05_m156_covariates.R"))
source(fs::path("analysis", "application", "main", "06_m156_interpretation.R"))
source(fs::path("analysis", "application", "main", "07_m156_prediction.R"))
