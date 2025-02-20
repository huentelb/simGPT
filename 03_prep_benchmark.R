# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PREPARE DATA AND GENERATE SEQUENCE ANALYSIS OUTPUT FOR BENCHMARKING
# 1960 BIRTH COHORT - AGE RANGE 0-67 

# huenteler@demogr.mpg.de
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cohort <- 1960
max_age <- 59
ages <- as.character(c(0:max_age))

# Generate folders to store results
# Need to be manually exported from Norwegian server containing the register data

# 1. Upper level folder based on simulation base_seed
folder.baseseed <- paste0(folder,"/sim_results_", supfile, "_",base_seed,"_/")
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



#### MERGE MULTIPLE SIMULATION OUTPUTS FOR 1960 COHORT ####

# load the gp-data for from each simulation
for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {

load(paste0(folder, "/sim_results_", supfile, "_",base_seed,i,"_/gp196059.RData"))
assign(paste0("gp_", i), gp)  

}

# merge the simulations into one gp-dataframe
gp <- rbind(gp_1, gp_2, gp_3, gp_4, gp_5, gp_6, gp_7, gp_8, gp_9, gp_10)

# store combined gp dataframe in base_seed folder
save(gp, file = paste0(folder.baseseed, "gp196059.RData"))



### last line ###