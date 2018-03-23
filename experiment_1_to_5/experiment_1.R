# =========================================================================== #
#                                EXPERIMENT 1
# =========================================================================== #

# 1 - Setup
# ============================================
# Load all functions related to the experiement
source("functions.R")

# load simulation parameters
source("simulation_setup.R")

# load model parameters
source("model_setup.R")

# 2 - Observations
# ============================================
# Generate observations 
source("observation_simulation.R")

# Compute observables
source("observables.R")

# 3 - Filter
# ============================================
# Run filter
source("filtering_procedure.R")
