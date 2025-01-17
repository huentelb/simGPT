
## LOAD / INSTALL NECESSARY PACKAGES

#install.packages("devtools")
library(devtools)
# install.packages("tidyverse")
library(tidyverse) #For data wrangling
# install.packages("ggh4x")
library(ggh4x) #for extended facet plotting
# install.packages("data.table") #for large datasets
library(data.table)

# remove.packages("rsocsim")
# devtools::install_github("MPIDR/rsocsim")
library(rsocsim)

# Load functions to write SOCSIM rate files
source("functions.R")

# convert HFD data to SOCSIM format
write_socsim_rates_HFD(Country = "NOR")

#convert HMD data to SOCSIM format
write_socsim_rates_HMD(Country = "NOR")


#### SIMULATION ####


## CREATE INITIAL POPULATION AND MARRIAGE FILES
# Set size of initial population
size_opop <-  10000

# Create data.frame with 14 columns and nrows = size_opop
presim.opop <- setNames(data.frame(matrix(data = 0, ncol = 14, nrow = size_opop)), 
                        c("pid","fem","group","nev","dob","mom","pop",
                          "nesibm","nesibp","lborn","marid","mstat","dod","fmult"))

# Add pid 1:sizeopop
presim.opop$pid <- 1:size_opop

# Add sex randomly
presim.opop$fem <- sample(0:1, nrow(presim.opop), replace = T)

# Add group 1 for all individuals
presim.opop$group <- 1

# Add random dates of birth (max age around 50)
presim.opop$dob <- sample(600:1200, nrow(presim.opop), replace = T)

# Write initial population for pre-simulation (without fertility multiplier)
write.table(presim.opop, "presim.opop", row.names = F, col.names = F)


## Create an empty data frame for presim.omar
presim.omar <- data.frame()

# Write empty omar for pre-simulation
write.table(presim.omar, "presim.omar", row.names = F, col.names = F)



## RUN SOCSIM SIMULATION

# start loop for simulation rounds
for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {


# Automatically set a new seed for each simulation in i based on base_seed
  seed <- paste0(base_seed,i)


# Record starting time 
start <- Sys.time()
# Run a single SOCSIM-simulation with a given folder and the provided supervisory file
# using the "future" process method
rsocsim::socsim(folder, supfile, seed, process_method = "future")
# Record ending time
end <- Sys.time()

# Store duration
duration <- difftime(end, start, unit = "hours")

## IMPORT OUTPUTS TO R
# Read the opop file using the read_opop function
opop <- rsocsim::read_opop(folder = getwd(), supfile = supfile, 
                           seed = seed, suffix = "",  fn = NULL)


# Read the omar file using the read_opop function
omar <- rsocsim::read_omar(folder = getwd(), supfile = supfile, 
                           seed = seed, suffix = "",  fn = NULL)



## ESTIMATE AGE-SPECIFIC RATES FROM SIMULATED DATA (5-year intervals)
# Estimate age-specific fertility rates
asfr_sim <- rsocsim::estimate_fertility_rates(opop = opop,
                                              final_sim_year = 2100, #[Jan-Dec]
                                              year_min = 1850, # Closed [
                                              year_max = 2100, # Open )
                                              year_group = 5, 
                                              age_min_fert = 15, # Closed [
                                              age_max_fert = 50, # Open )
                                              age_group = 5) #[,)

# Estimate age-specific mortality rates
asmr_sim <- rsocsim::estimate_mortality_rates(opop = opop, 
                                              final_sim_year = 2100, # [Jan-Dec]
                                              year_min = 1850, # Closed
                                              year_max = 2100, # Open )
                                              year_group = 5,
                                              age_max_mort = 100, # Open )
                                              age_group = 5) #[,)



## COMPARE SIMULATED RATES TO INPUT RATES

# FERTILITY

# Extract year and age breaks used in the estimate_fertility_rates() to apply the same values to HFD data
# Year breaks. Extract all the unique numbers from the intervals. 
year_breaks_fert <- unique(as.numeric(str_extract_all(asfr_sim$year, "\\d+", simplify = T)))

# Year range to filter HFD data
year_range_fert <- min(year_breaks_fert):max(year_breaks_fert-1)

# Age breaks of fertility rates. Extract all the unique numbers from the intervals 
age_breaks_fert <- unique(as.numeric(str_extract_all(asfr_sim$age, "\\d+", simplify = T)))

# Read the same HFD file used to write the input fertility files
HFD <- read.table(file = input_hfd, 
                  as.is=T, header=T, skip=0, stringsAsFactors=F)

# Wrangle data and compute monthly fertility rates
HFD <- HFD %>% 
  mutate(Age = case_when(Age == "14-" ~ "14",
                         Age == "50+" ~ "50", 
                         TRUE ~ Age),
         Age = as.numeric(Age)) %>% 
  filter(Year %in% year_range_fert) %>% 
  mutate(year = cut(Year, breaks = year_breaks_fert, 
                    include.lowest = F, right = F, ordered_results = T),
         age = cut(Age, breaks = age_breaks_fert, 
                   include.lowest = F, right = F, ordered_results = T)) %>% 
  filter(!is.na(age)) %>% 
  group_by(year, age) %>%
  summarise(ASFR = mean(ASFR)) %>%
  ungroup() %>%
  mutate(Source = "HFC/HMD/UNWPP", 
         Rate = "ASFR")


# Wrangle SOCSIM data
SocsimF <- asfr_sim %>% 
  rename(ASFR = socsim) %>% 
  mutate(Source = "SOCSIM",
         Rate = "ASFR")

yrs_plot <- c("[1850,1855)", "[1880,1885)", "[1910,1915)", "[1940,1945)", "[1970,1975)", "[2000,2005)", "[2020,2025)") 


# MORTALITY

# Extract year and age breaks used in the estimate_mortality_rates() to apply the same values to HMD data
# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort <- unique(as.numeric(str_extract_all(asmr_sim$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort <- min(year_breaks_mort):max(year_breaks_mort-1)

# Age breaks of mortality rates. Extract all the unique numbers from the intervals 
age_breaks_mort <- unique(as.numeric(str_extract_all(asmr_sim$age, "\\d+", simplify = T)))


# Read the same HMD female life table 1x1 file used to write the input mortality files
ltf <- read.table(file=input_hmd_f,
                  as.is=T, header=T, skip=0, stringsAsFactors=F)

# Read the same HMD male life table 1x1 file used to write the input mortality files
ltm <- read.table(file=input_hmd_m,
                  as.is=T, header=T, skip=0, stringsAsFactors=F)

# Wrangle HMD life tables
HMD <- ltf %>%
  select(Year, Age, mx) %>%
  mutate(Sex = "Female") %>% 
  bind_rows(ltm %>% 
              select(Year, Age, mx) %>%  
              mutate(Sex = "Male")) %>% 
  mutate(Age = as.numeric(Age)) %>% 
  filter(Year %in% year_range_mort) %>% 
  mutate(year = cut(Year, breaks = year_breaks_mort, 
                    include.lowest = F, right = F, ordered_results = T),
         age = cut(Age, breaks = age_breaks_mort, 
                   include.lowest = F, right = F, ordered_results = T)) %>%
  filter(!is.na(age)) %>% 
  group_by(year, Sex, age) %>% 
  summarise(mx = mean(mx)) %>%
  ungroup() %>% 
  mutate(Source = "HFC/HMD/UNWPP",
         Rate = "ASMR")

# Wrangle SOCSIM data
SocsimM <- asmr_sim %>% 
  rename(mx = socsim) %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Source = "SOCSIM",
         Rate = "ASMR") %>% 
  select(year, Sex, age,  mx, Source, Rate)

# yrs_plot <- c("[1850,1851)", "[1880,1881)", "[1910,1911)", "[1940,1941)", "[1970,1971)", "[2000,2001)", "[2020,2021)") 
yrs_plot <- c("[1850,1855)", "[1880,1885)", "[1910,1915)", "[1940,1945)", "[1970,1975)", "[2050,2055)", "[2020,2025)", "[2095,2100)") 




## PLOT SIMULATED AND INPUT DATA TOGETHER

# create folder to store graphs
graph.folder <- paste0(folder,"/sim_results_", supfile, "_",seed,"_/graphs/") # Check if the folder already exists
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}

## Plotting ASFR and ASMR (for females) from HFC/HMD/UNWPP vs SOCSIM 
yrs_plot <- c("[1850,1855)", "[1880,1885)", "[1910,1915)", "[1940,1945)", "[1970,1975)", "[2050,2055)", "[2020,2025)", "[2095,2100)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_sim$age)

png(file = paste0(graph.folder,"rates_f.png"),
    width = 964, height = 556)
print(
  bind_rows(HFD %>% rename(Estimate = ASFR), 
            SocsimF %>% rename(Estimate = ASFR)) %>% 
    mutate(Sex = "Female") %>%   
    bind_rows(HMD %>% rename(Estimate = mx),
              SocsimM %>% rename(Estimate = mx)) %>% 
    # Filter rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
    filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
    filter(Sex == "Female") %>% 
    mutate(age = factor(as.character(age), levels = age_levels)) %>%
    filter(year %in% yrs_plot) %>% 
    ggplot(aes(x = age, y = Estimate, group = interaction(year, Source)))+
    facet_wrap(. ~ Rate, scales = "free") + 
    geom_line(aes(colour = year, linetype = Source), linewidth = 1)+
    scale_linetype_manual(values = c("HFC/HMD/UNWPP" = "solid","SOCSIM" = "11")) +
    facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                                 ASMR =  scale_y_continuous(trans = "log10")))+
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(title = paste0("Age-Specific Fertility and Mortality rates of Women in Norway"),
         subtitle = paste0("Retrieved from HFC, HMD and a SOCSIM simulation \n ", het," heterogeneous fertility, ", bint,  ", opop size = ", size_opop, "\n alpha = ", alpha, ", beta = ", beta, " (", seed, "), estimation time: ", round(duration, 2), " hours"), 
         x = "Age") + 
    theme_bw()
  )
dev.off()

png(file = paste0(graph.folder,"rates_m.png"),
    width = 564, height = 556)
print(
  bind_rows(HMD %>% rename(Estimate = mx),
            SocsimM %>% rename(Estimate = mx)) %>% 
    # Filter rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
    filter(Estimate != 0 & !is.infinite(Estimate) & !is.nan(Estimate)) %>% 
    filter(Sex == "Male") %>% 
    mutate(age = factor(as.character(age), levels = age_levels)) %>%
    filter(year %in% yrs_plot) %>% 
    ggplot(aes(x = age, y = Estimate, group = interaction(year, Source)))+
    facet_wrap(. ~ Rate, scales = "free") + 
    geom_line(aes(colour = year, linetype = Source), linewidth = 1)+
    scale_linetype_manual(values = c("HFC/HMD/UNWPP" = "solid","SOCSIM" = "11")) +
    facetted_pos_scales(y = list(ASMR =  scale_y_continuous(trans = "log10")))+
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(title = paste0("Age-Specific Mortality rates of Men in Norway"),
         subtitle = paste0("Retrieved from HMD and a SOCSIM simulation \n ", het," heterogeneous fertility, ", bint,  ", opop size = ", size_opop, "\n alpha = ", alpha, ", beta = ", beta, " (", seed, "), estimation time: ", round(duration, 2), " hours"), 
         x = "Age") + 
    theme_bw()
  )
dev.off()



## GETTING ESTIMATES -- IDENTIFY KIN! ####

# Use the final simulation year (January-December)
# and the last simulated month to convert our monthly dates into yearly ones
asYr <- function(month, last_month=last_month, final_sim_year=final_sim_year) {
  return(final_sim_year - trunc((last_month - month)/12))
}

## Read the opop file using the read_opop function
# opop <- rsocsim::read_opop(folder = getwd(), supfile = "socsim_NOR.sup", 
#                            seed = seed, suffix = "",  fn = NULL)
# 
# 
# ## Read the omar file using the read_opop function
# omar <- rsocsim::read_omar(folder = getwd(), supfile = "socsim_NOR.sup", 
#                            seed = seed, suffix = "",  fn = NULL)

# Parameters specific to this simulation: will need to be changed
last_month <- max(opop$dob) # Last simulated month

# Cleaning our population file
opop <- opop %>% 
  #Fixing dates of death for individuals still living at the end
  mutate(dod = if_else(dod == 0, last_month, dod)) %>%
  #Dates of birth and death in years for both individual and parents
  mutate(dob_year = asYr(dob, last_month=last_month, final_sim_year=final_sim_year),
         dod_year = asYr(dod, last_month=last_month, final_sim_year=final_sim_year))

# Saving opop data frame for usage in R on my Mac
save(opop, file = paste0(folder,"/sim_results_", supfile, "_",seed,"_/opop.RData"))




## PREPARATION SEQUENCE ANALYSIS ##

# load simulated register data
load(paste0(folder, "/sim_results_", supfile, "_", seed, "_/opop.RData"))

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

# Start loop for selecting different cohorts
for (c in c(1953, 2000)) {

# Set selected cohort and maximum age to be displayed in SA
# size of 2000 birth cohort in Norway: 59,234
cohort <- c
max_age <- 100 # (for censoring of life courses)

# Select sample (based on birth cohort)
sample <- opop %>% 
  # keep if born in 1953    
  # & survived at least to first birthday (may produce inconsistent SA results if not)
  filter(dob_year == cohort & (dod_year - dob_year > 0))

#### PARENTS ####
## Build register containing death dates for parents for sample individuals

# 1. Register containing mother ids and their vital dates
mothers <- sample %>% 
  # join to mother identifier (mom) in sample (left) data from opop (right) 
  left_join(opop, by = c("mom" = "pid"), suffix = c("", "_p"), keep = TRUE) %>% 
  # keep original variables plus parental death years 
  select(pid, mom, pop, pid_p, dod_year_p) %>% 
  mutate(mother = 1)

# 2. Register containing father ids and their vital dates
fathers <- sample %>% 
  # join to father identifier (pop) in sample (left) data from opop (right) 
  left_join(opop, by = c("pop" = "pid"), suffix = c("", "_p"), keep = TRUE) %>% 
  # keep original variables plus parental death years 
  select(pid, mom, pop, pid_p, dod_year_p) %>% 
  mutate(mother = 0)

## Build one parent register for ALL sample individuals
parents <- bind_rows(mothers, fathers)
# ascending order by pid
parents <- arrange(parents, pid)


## Build register containing the parent who died last
lateparent <- parents %>% 
  group_by(pid) %>% 
  slice_max(dod_year_p, with_ties = FALSE) %>% 
  ungroup %>% 
  # drop id of parent who died first
  dplyr::select(pid, pid_p, dod_year_p, mother) %>% 
  
  # generate indicator if parent has died during observation period (dead_p)
  # set dead_p to 0 if died after end of obs period (=cohort+max_age)
  # useful if life course is censored (e.g., for benchmarking against empirical register)
  mutate(dead_p = ifelse(dod_year_p > cohort+max_age, 0, 1))




#### CHILDREN ####
## Build register containing children of FEMALE sample individuals
children_f <- sample %>% 
  # select females
  filter(fem == 1) %>% 
  # merge children's vars to the sample individuals who became their mothers
  # i.e. sample individual's id 'pid' equals the samples children's mother id 'mom'
  left_join(opop, by = c("pid" = "mom"), suffix = c("", "_c"), keep = TRUE) %>% 
  # keep original variables plus child birth- and death-years 
  select(pid,  
         pid_c, fem_c, mom_c, pop_c, dob_year_c, dod_year_c)
      
## Build register containing children of MALE sample individuals
children_p <- sample %>% 
  # select males
  filter(fem == 0) %>% 
  # merge children's vars to the sample individuals who became their fathers
  # i.e. sample individual's id 'pid' equals the samples children's father id 'pop'
  left_join(opop, by = c("pid" = "pop"), suffix = c("", "_c"), keep = TRUE) %>% 
  # keep original variables plus child birth- and death-years 
  select(pid,  
         pid_c, fem_c, mom_c, pop_c, dob_year_c, dod_year_c) 

## Build one child register for ALL sample individuals
children <- bind_rows(children_f, children_p)
# ascending order by pid
children <- arrange(children, pid)

# sample individuals without children
nochildren <- children %>% 
  filter(is.na(pid_c)) %>% 
  select(pid)

## Build register containing oldest child of sample individuals
oldestc <- children %>%
  group_by(pid) %>% 
  # add number of all children per sample individual before deleting younger siblings
  add_count(pid, name = "numkids") %>% 
  # select child that was born first (do not allow multiple choices per sample individual)
  slice_min(dob_year_c, with_ties = FALSE) %>%
  ungroup() %>% 
  # set numkids to NA if person has no children (i.e. child's yob is NA)
  mutate(numkids = ifelse(is.na(dob_year_c), NA, numkids)) 



#### GRANDCHILDREN ####
## Build register containing children of FEMALE children --> sample's grandchildren
gchildren_f <- children %>% 
  # select female children
  filter(fem_c == 1) %>% 
  # merge grandchildren's vars to the sample's daughters who became their mothers
  # i.e. sample children's id 'pid_c' equals the sample grandchildren's mother id 'mom'
  left_join(opop, by = c("pid_c" = "mom"), suffix = c("", "_gc"), keep = TRUE) %>% 
  # rename grandchild variables by adding suffix _gc
  rename(mom_gc = mom) %>% 
  rename(pop_gc = pop) %>% 
  rename(fem_gc = fem) %>% 
  rename(dob_year_gc = dob_year) %>% 
  rename(dod_year_gc = dod_year) %>% 
  # keep original variables plus grandchild identifiers and birth and death years 
  select(pid, 
         pid_gc, mom_gc, pop_gc, fem_gc, dob_year_gc, dod_year_gc)


## Build register containing children of MALE children --> sample's grandchildren
gchildren_p <- children %>% 
  # select male childre
  filter(fem_c == 0) %>% 
  # merge grandchildren's vars to the sample's sons who became their fathers
  # i.e. sample children's id 'pid_c' equals the sample grandchildren's father id 'pop'
  left_join(opop, by = c("pid_c" = "pop"), suffix = c("", "_gc"), keep = TRUE) %>% 
  # rename grandchild variables by adding suffix _gc
  rename(mom_gc = mom) %>% 
  rename(pop_gc = pop) %>% 
  rename(fem_gc = fem) %>% 
  rename(dob_year_gc = dob_year) %>% 
  rename(dod_year_gc = dod_year) %>% 
  # keep original variables plus grandchild identifiers and birth and death years 
  select(pid, 
         pid_gc, mom_gc, pop_gc, fem_gc, dob_year_gc, dod_year_gc)

## Build one grandchild register for ALL sample individuals (incl those with no children)
gchildren <- bind_rows(gchildren_f, gchildren_p, nochildren)
# ascending order by pid
gchildren <- arrange(gchildren, pid)

## Build register containing oldest child of sample individuals
oldestgc <- gchildren %>%
  group_by(pid) %>% 
  # add number of all children per sample individual
  add_count(pid, name = "numgkids") %>% 
  # select child that was born first
  slice_min(dob_year_gc, with_ties = FALSE) %>%
  ungroup() %>% 
  # set numgkids to NA if person has no grandchildren (i.e. grandchild's yob is NA)
  mutate(numgkids = ifelse(is.na(dob_year_gc), NA, numgkids)) 

#### ONE FAMILY ROOSTER ####
gp_help <- left_join(sample, lateparent)
gp_help <- left_join(gp_help, oldestc)
gp_help <- left_join(gp_help, oldestgc)

# ...only with info needed for generating the generational placement trajectories
gp_base <- gp_help %>% 
  select(pid, dob_year, dod_year, dod_year_p, dob_year_c, dob_year_gc) %>% 
  # right-censor at age 100 (that's when register is censored) 
  mutate(dod_year_p = ifelse(dod_year_p > cohort+max_age, NA, dod_year_p),
         dob_year_c = ifelse(dob_year_c > cohort+max_age, NA, dob_year_c),
         dob_year_gc = ifelse(dob_year_gc > cohort+max_age, NA, dob_year_gc),
         dod_year = ifelse(dod_year > cohort+max_age, NA, dod_year)) %>% 
  # convert years into ages
  mutate(pdage = dod_year_p - dob_year,
         cage = dob_year_c - dob_year,
         gcage = dob_year_gc - dob_year,
         dage = dod_year - dob_year) 


## Average age at family transitions
avg_par <- mean(gp_base$pdage, na.rm = TRUE)
avg_child <- mean(gp_base$cage, na.rm = TRUE)
avg_gchild <- mean(gp_base$gcage, na.rm = TRUE)
avg_dead <- mean(gp_base$dage, na.rm = TRUE)



#### GP DATA SET ####

# 1) Make long dataframe with one row for each age in the life course
# 2) Fill rows with numbers indicating gp by kin type depending on age at gp-transitions 
# 3) Create "sum" for gp to have one summary indicator

n_sample  <- nrow(sample)
ages <- as.character(c(0:max_age))
lab_ages <- paste0("age", ages)
empty_cols <- data.frame(matrix(NA, nrow = n_sample, ncol = max_age+1))
colnames(empty_cols) <- lab_ages

# 1)
gp_long <- gp_base %>% 
  select(pid, dage, pdage, cage, gcage) %>% 
  # duplicate each sample individual times the age range we consider to create long format
  slice(rep(row_number(), each = max_age+1)) %>% 
  # generate age var and start with 0
  group_by(pid) %>% 
  mutate(age = (row_number()-1)) %>% 
  
  # set ages after parents have died to 2 (1 = at least one parent; 2 = NO parent) 
  # if parental death not yet happened to 1
  mutate(gpp = ifelse(age >= pdage & !is.na(pdage), 2, 1),
         
         # 2)
         # set ages after child(ren) were born to 10 (10 = at least one child; 20 = NO child)
         # if childless and before to 20
         gpc = ifelse(age >= cage & !is.na(cage), 10, 20),
         
         # set ages after grandchild(ren) were born to 100 (100 = at least one grandchild; 200 = NO grandchild)
         # if grandchildless and before to 200
         gpg = ifelse(age >= gcage & !is.na(gcage), 100, 200),
         
         # 3)
         # generate "sum" of generational placements
         gp = gpp + gpc + gpg,
         
         # add status "dead" to gp
         gp = ifelse(age >= dage & !is.na(dage), 0, gp)) %>% 
  ungroup


# Reshape to wide dataframe to obtain STS format (STates Sequences; used by TraMineR)
library(tidyr)
gp_wide <- gp_long %>%
  pivot_wider(
    id_cols = c(pid, dage, pdage, cage, gcage),
    names_from = age,
    values_from = gp,
    names_glue = "{.value}_{age}"
  )



# Transform numerics into characters
# mapping dictionary (attention: logic between numerics and characters is inverted (P1xxxx or P0xxxx refers to xx1 or xx2))
mapping <- c("D" = 0,
             "P1C0G0" = 221, "P0C0G0" = 222, 
             "P1C1G0" = 211, "P0C1G0" = 212, 
             "P1C1G1" = 111, "P0C1G1" = 112)

gp <- gp_wide %>% 
  mutate(across(starts_with("gp"), ~case_when(
    . == mapping[1] ~ "D",
    . == mapping[2] ~ "P1C0G0",
    . == mapping[3] ~ "P0C0G0",
    . == mapping[4] ~ "P1C1G0",
    . == mapping[5] ~ "P0C1G0",
    . == mapping[6] ~ "P1C1G1",
    . == mapping[7] ~ "P0C1G1",
    TRUE ~ as.character(.)
  )))

# add other demographic information from sample individuals plus family to gp-dataframe
gp <- left_join(gp, lateparent)
gp <- left_join(gp, oldestc)
gp <- left_join(gp, oldestgc)
gp <- left_join(gp, sample)

save(gp, file = paste0(folder,"/sim_results_", supfile, "_",seed,"_/gp",cohort,max_age,".RData"))

# end cohort selection loop
}

# end simulation round loop
}

#### MERGE MULTIPLE SIMULATION OUTPUTS ####

load(paste0(folder, "/sim_results_", supfile, "_",base_seed,"1_/gp00.RData"))
gp_1 <- gp

#### CONTINUE HERE!!! ####

load(paste0(folder, "/sim_results_", supfile, "_2403083_/gp00.RData"))
gp_3 <- gp
load(paste0(folder, "/sim_results_", supfile, "_2403084_/gp00.RData"))
gp_4 <- gp

gp <- rbind(gp_2, gp_3, gp_4)

save(gp, file = paste0(folder,"/sim_results_", supfile, "_",seed,"4_/gp00.RData"))

load(paste0(folder, "/sim_results_", supfile, "_2403084_/gp00.RData"))

# create folder to store graphs
graph.folder <- paste0(folder,"/sim_results_", supfile, "_",seed,"_/graphs/",cohort,"_",max_age,"/")# Check if the folder already exists
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}


#### SEQUENCE ANALYSIS ####
library(foreign)
library(TraMineR)
library(TraMineRextras)
library(WeightedCluster)
library(RColorBrewer)
library(bookdown)

# Define characteristics of sequences

# Alphabet, ie all possible states
gpalpha <- c("P1C0G0", "P0C0G0", 
             "P1C1G0", "P0C1G0", 
             "P1C1G1", "P0C1G1")

# Labels of states
gplabels <- c("child", "no ancestors/descentants",
              "child and parent", "parent",
              "child, parent, and grandparent", 
              "parent and grandparent")

# Short labels for states
gpstates <-c("G1G2", "G2",
             "G1G2G3", "G2G3",
             "G1G2G3G4", "G2G3G4")

# Color palette
cblind <- brewer.pal(6, "RdBu")


# Define sequence object using the above characteristics
seq <- seqdef(gp, 6:106, 
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates,
              missing = "D", right = "DEL")

table(gp$gp_99)
# Plot state distribution of life course
#seqIplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE)
png(file = paste0(graph.folder, "seqD_full.png"),
    width=964, height=556)
seqdplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - ",cohort+max_age,", ",cohort," birth cohort), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", with.missing = T)
dev.off()


png(file = paste0(graph.folder, "seqI_full.png"),
    width=964, height=556)
seqIplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - ",cohort+max_age,", ",cohort," birth cohort), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
                       missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqi100_full.png"),
    width=964, height=556)
seqiplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - ",cohort+max_age,", ",cohort," birth cohort), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", idxs = 1:100)
dev.off()


# Meantime in general pop
seqmeant(seq)
png(file = paste0(graph.folder, "meant_full.png"),
    width=964, height=556)
seqmtplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - ",cohort+max_age,", ",cohort," birth cohort), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop))
dev.off()



#### CLUSTER ANALYSIS ####

# 1. Generate dissimilarity matrix

# Chi2 distance
# Chi2 distance with number of periods *K* set to length of sequence 
# --> makes Chi2 to a "position-wise" Chi2-measure, sensitive to timing (similar to HAM)
chi <- seqdist(seq, method = "CHI2", step = max(seqlength(seq)))

# OM distance with transition rate based costs
#omt <- seqdist(seq, method = "OM", indel = 1, sm = "TRATE")

# OM distance with constant costs
#omc <- seqdist(seq, method = "OM", indel = 1, sm = "CONSTANT")


# 2. Clustering

# 2a. Hierarchical clustering using WARD (1-10 clusters)
chi_ward <- hclust(as.dist(chi), method = "ward.D")
#omt_ward <- hclust(as.dist(omt), method = "ward.D")
#omc_ward <- hclust(as.dist(omc), method = "ward.D")

chi_ward10 <- as.clustrange(chi_ward, diss = chi, ncluster = 10)
chiWard.qual <- chi_ward10
plot(chiWard.qual, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2)

# omt_ward10 <- as.clustrange(omt_ward, diss = chi, ncluster = 10)
# omtWard.qual <- omt_ward10
# plot(omtWard.qual, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2)
# 
# omc_ward10 <- as.clustrange(omc_ward, diss = chi, ncluster = 10)
# omcWard.qual <- omc_ward10
# plot(omcWard.qual, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2)

# 2b. Partitioning around medoids
# PAM + Ward (as starting point)

chi_pam10 <- wcKMedRange(chi,
                         kvals = 2:10,
                         initialclust = chi_ward)

# Plot cluster quality
png(file = paste0(graph.folder, "clustqual.png"),
    width=964, height=556)
plot(chi_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2, 
     main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                   "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
dev.off()
summary(chi_pam10, max.rank = 3)

# omt_pam10 <- wcKMedRange(omt,
#                          kvals = 2:10,
#                          initialclust = omt_ward)
# plot(omt_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2)
# summary(omt_pam10, max.rank = 3)
# 
# omc_pam10 <- wcKMedRange(omc,
#                          kvals = 2:10,
#                          initialclust = omc_ward)
# plot(omc_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2)
# summary(omc_pam10, max.rank = 3)


#### D plots ####
# State distribution plots by clusters
# CHI2

w = 2000
h = 1250

png(file = paste0(graph.folder, "seqD_4.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster4, 
         border = NA, ltext = c(gpstates), with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_5.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_6.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
        border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
        main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                      "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
        missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_7.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

#### F plots ####
png(file = paste0(graph.folder, "seqF20_4.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster4,
         border = NA, ltext = c(gpstates),  with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_5.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_6.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_7.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

# OM trate
# seqdplot(seq, group = omt_pam10$clustering$cluster4,
#          border = NA, ltext = c(gpstates), main = "OM trate PAM: 4 Clusters")
# 
# seqdplot(seq, group = omt_pam10$clustering$cluster5,
#          border = NA, ltext = c(gpstates), main = "OM trate PAM: 5 Clusters")
# 
# seqdplot(seq, group = omt_pam10$clustering$cluster6,
#          border = NA, ltext = c(gpstates), main = "OM trate PAM: 6 Clusters")

# OM constant
# seqdplot(seq, group = omc_pam10$clustering$cluster4,
#          border = NA, ltext = c(gpstates), main = "OM constant PAM: 4 Clusters")
# 
# seqdplot(seq, group = omc_pam10$clustering$cluster5,
#          border = NA, ltext = c(gpstates), main = "OM constant PAM: 5 Clusters")
# 
# seqdplot(seq, group = omc_pam10$clustering$cluster6,
#          border = NA, ltext = c(gpstates), main = "OM constant PAM: 6 Clusters")





# Mean time spent in each state by cluster
by(seq, chi_pam10$clustering$cluster5, seqmeant)
png(file = paste0(graph.folder, "mean_plot_6.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster6, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 6 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "mean_plot_5.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster5, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 5 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
          missing.color = "#f7f7f7")
dev.off()

#### LABELLED GRAPH FOR OPTIMAL CLUSTER SOLUTION ####

# We extract X clusters and re-label them from 1 to X to replace the medoid identifiers

# identify medoids sorted by frequency
mc <- chi_pam10$clustering$cluster5
med <- as.data.frame(sort(table(mc), decreasing = TRUE))
med1 <- as.character(med[1,1])
med2 <- as.character(med[2,1])
med3 <- as.character(med[3,1])
med4 <- as.character(med[4,1])
med5 <- as.character(med[5,1])
# med6 <- as.character(med[6,1])
# med7 <- as.character(med[7,1])
# create factor containing medioids incl labels
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5),
                    as.character("1","2","3","4","5"),
                    labels = c("Cluster 1 -\n 3-gen family",
                               "Cluster 2 -\n Childless",
                               "Cluster 3 -\n 4-gen family",
                               "Cluster 4 -\n 3-gen (via 2-gen) family",
                               "Cluster 5 -\n 2-gen family"))

# attach to dataframe to use as weights in plots
gp$chi <- mc.factor

# generate new sequence object 
# seqchi <- seqdef(gp, 6:106,
#                  labels = gplabels,  
#                  cnames = ages, 
#                  tick.last = TRUE, 
#                  xtstep = 5, 
#                  cpal = cblind, 
#                  alphabet = gpalpha, 
#                  states = gpstates,
#                  missing = "D", right = "DEL")

# different plots with labels
png(file = paste0(graph.folder, "seqD_5_lab.png"),
    width=w, height=h)
seqdplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
dev.off()

png(file = paste0(graph.folder, "seqI_5_lab.png"),
    width=w, height=h)
seqIplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqf20_5_lab.png"),
    width=w, height=h)
seqfplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "mean_plot_5_lab.png"),
    width=w, height=h)
seqmtplot(seq, group = group.p(gp$chi), border = NA,
          ltext = c(gpstates), 
          main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                        "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "seqr_5_lab.png"),
    width=w, height=h)
seqrplot(seq, group = group.p(gp$chi), border = NA,
          ltext = c(gpstates), 
          main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                        "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE, diss = chi)
dev.off()

# seqfplot(seqchi, group = group.p(gp$chi), idxs = 1:50,
#          ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
#          ylab = NA, yaxis = FALSE, border = NA)
# 
# seqiplot(seqchi, group = group.p(gp$chi), idxs = 1:500,
#          ltext = gpstates, use.layout = TRUE, cex.legend = 1.2, space = 0,
#          ylab = NA, yaxis = FALSE)

# cross-sectional entropy plot
# seqHtplot(seqchi, group = group.p(gp$chi),
#          ltext = gpstates)

#### TABLES ####
library(gtsummary)
library(vtable)

# merge parent alive until end of observation, number of kids and number of grandkids to data
# clusters <- left_join(gp, lateparent)
# clusters <- left_join(clusters, oldestc)
# clusters <- left_join(clusters, oldestgc)
# clusters <- left_join(clusters, sample)


descr <- gp %>% 
  select(pid, dage, pdage, cage, gcage, chi, dead_p, numkids, numgkids, fem, mstat) %>% 
  # set (grand)parent status to 0 if number of (grand)kids missing, 1 otherwise
  mutate(isparent = ifelse(is.na(numkids), 0, 1),
         isgparent = ifelse(is.na(numgkids), 0, 1)) %>% 
  # convert marital status var to factor
  mutate(mstat = factor(mstat, levels = c(1:4),
            labels = c("Single",
                       "Divorced",
                       "Widowed",
                       "Married"))) 

# Save dataframe for comparing with other simulations later
# save(descr, file = paste0(folder, "/sim_results_", supfile, "_", seed, "_/sumfile_",seed,".RData"))
# dev.off()

# Descriptive table using sumtable

st(descr, 
   vars = c("isparent", "numkids", "cage", 
            "isgparent", "numgkids", "gcage", 
            "dead_p", "pdage",
            "fem", "mstat", "dage"),
   summ=c("mean(x)",
          "sd(x)"),
   summ.names = c("Mean",
                  "Std.Dev"),
   group = "chi", 
   group.test =TRUE,
   labels = c(paste0("Parent at age ", max_age), "Number of children", "Age at birth of first child",
              paste0("Grandparent at age ", max_age), "Number of grandchildren", "Age at birth of first grandchild",
              paste0("Both parents dead at age ", max_age), "Age at death of second parent", 
              "Female", "Marital status", "Age at own death"),
   title = paste0("Summary Statistics (Simulated data, ", het, " heterogeneous fertility, ", bint,  ", opop size = ", seed, ")"),
   out = "csv",
   file = paste0(graph.folder, "sumtable_5"))


