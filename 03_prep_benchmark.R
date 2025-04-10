# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PREPARE DATA AND GENERATE SEQUENCE ANALYSIS OUTPUT FOR BENCHMARKING
# 1953 BIRTH COHORT - AGE RANGE 0-67 

# huenteler@demogr.mpg.de
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)


cohort <- 1953
max_age <- 66
ages <- as.character(c(0:max_age))
base_seed <- 123
folder <- "N:\\durable\\Data20\\Project_JW_FamWealth\\simgpt/"

# Generate folders to store results
# Need to be manually exported from Norwegian server containing the register data

# 1. Upper level folder based on simulation base_seed
folder.baseseed <- paste0(folder,"/sim_results_",base_seed, "/")
if (!dir.exists(folder.baseseed)) {
  # If not, create the new folder
  dir.create(folder.baseseed)
  cat("Folder created:", folder.baseseed, "\n")
} else {
  cat("Folder already exists:", folder.baseseed, "\n")
}

# 2. Lower level folder for cohort-max age-specific output 
graph.folder <- paste0(folder.baseseed, "benchmark/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}





#### MERGE SIMULATION OUTPUTS FOR 1953 COHORT WITH EMPIRICAL DATA ####

# load the gp-data for from the simulation
load("N:\\durable\\Data20\\Project_JW_FamWealth\\simgpt\\data\\gp195366.RData")
gpm <- gp


#### EMPIRICAL REGISTER DATA #####

# Load
library(foreign)
gp <- read.dta("N:\\durable\\Data20\\Project_JW_FamWealth\\simgpt\\data\\gp_wide.dta", 
               convert.factors = TRUE)

# store as gp-dataframe with suffix e (empirical)
gpe <- gp %>% 
  select(id, gp1953:gp2019, dage, bage_c, bage_gc, numkids, numgkids, dead_p)

# add new variable "group" to the df for later comparison with simulated data
gpe <- cbind(gpe, group = 2)



##### MERGE EMPIRICAL and MICROSIMULATED DATA ####
gp <- gpe %>% 
  rbind(gpm)


### last line ###