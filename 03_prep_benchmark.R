# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PREPARE DATA AND GENERATE SEQUENCE ANALYSIS OUTPUT FOR BENCHMARKING
# 1953 BIRTH COHORT - AGE RANGE 0-67 

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)


cohort <- 1960
max_age <- 59
ages <- as.character(c(0:max_age))

# load the gp-data for from each simulation from each birth cohorts
for (c in c(1960)) {
  
  for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {
    
    load(paste0(folder, "/sim_results_",base_seed,i,"_/gp", c, max_age, ".RData"))
    assign(paste0("gp_", i), gp)  
  }
  
  # merge the simulations into one gp-dataframe for each birth cohort
  gp <- rbind(gp_1, gp_2, gp_3, gp_4, gp_5, gp_6, gp_7, gp_8, gp_9, gp_10)
  
  
  # store combined gp dataframe in base_seed folder
  save(gp, file = paste0(folder.baseseed, "gp", c, max_age, ".RData"))
  
}



### last line ###