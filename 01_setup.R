# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SET-UP FOR MICROSIMULATION

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Set-up for microsimulations
# Always make necessary adjustments BEFORE starting a new round of simulations here

## 1. Working directory
setwd("/Users/bhuenteler/Library/CloudStorage/OneDrive-DIWBerlin/projects/GenerationalPlacements/SimGPT/analysis/simGPT")


# 2. Input files (usually no need to change)
input_hfd <- "Input/NORasfrRR_med.txt" # fertility
input_hmd_f <- "Input/fltper_1x1_med.txt" # female mortality
input_hmd_m <- "Input/mltper_1x1_med.txt" # male mortality

# 2.1 Scenario
scenario <-  "Fmed_Mmed" # Fertility_Mortality

## 3. Adjust to specifications you set in sup-file so that labels contain right info
het <- "no" # set to no or leave empty, depending on setting in sup
bint <- "birth int. 9"
alpha <- 0 # default 0
beta <- 1 # default 1

# 4. Base seed (within the code, we will run i different simulations that will then automatically be starting with the base_seed + i. 
# This will create one output folder for each i simulation)
base_seed <- "251113" # usually set to date you conduct the simulation

# 5. Set the size of the starting population
size_opop <-  10000


# 6. Specify the folder where the supervisory and the rate files are. 
# If the R session is running through the project, you can use the following command. 
folder <- getwd()



# Type the name of the supervisory file  stored in the above folder:
supfile <- "socsim_NOR.sup"


final_sim_year <- 2100

### last line ###