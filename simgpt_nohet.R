
# install.packages("devtools")
# install.packages("tidyverse")
library(tidyverse) #For data wrangling
# install.packages("ggh4x")
library(ggh4x) #for extended facet plotting
# install.packages("data.table") #for large datasets
library(data.table)

# remove.packages("rsocsim")
# devtools::install_github("MPIDR/rsocsim")
library(rsocsim)

## set wd
setwd("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT")

# Load functions to write SOCSIM rate files
source("functions.R")

# # convert HFD data to SOCSIM format
# write_socsim_rates_HFD(Country = "NOR") 
# 
# #convert HMD data to SOCSIM format
# write_socsim_rates_HMD(Country = "NOR")

#### SIMULATION ####

## Adjust to specifications you set in sup-file!

het <- "" # set to no or leave empty, depending on setting in sup
bint <- "birth int. 9"
alpha <- 1
beta <- 0.3

## CREATE INITIAL POPULATION AND MARRIAGE FILES
# Set size of initial population
size_opop <-  5000

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
for (i in c(1, 1)) {
# Specify the folder where the supervisory and the rate files are. 
# If the R session is running through the project, you can use the following command. 
folder <- getwd()

# Type the name of the supervisory file  stored in the above folder:
supfile <- "socsim_NOR.sup"

#### Choose a seed number (today's date) for the simulation: ####

  seed <- paste0("240228",i)

# Run a single SOCSIM-simulation with a given folder and the provided supervisory file
# using the "future" process method
rsocsim::socsim(folder, supfile, seed, process_method = "future")


## IMPORT OUTPUTS TO R
# Read the opop file using the read_opop function
opop <- rsocsim::read_opop(folder = getwd(), supfile = supfile, 
                           seed = seed, suffix = "",  fn = NULL)


# Read the omar file using the read_opop function
omar <- rsocsim::read_omar(folder = getwd(), supfile = supfile, 
                           seed = seed, suffix = "",  fn = NULL)



## ESTIMATE AGE-SPECIFIC RATES FROM SIMULATED DATA
# Estimate age-specific fertility rates
asfr_sim <- rsocsim::estimate_fertility_rates(opop = opop,
                                              final_sim_year = 2022, #[Jan-Dec]
                                              year_min = 1847, # Closed [
                                              year_max = 2022, # Open )
                                              year_group = 5, 
                                              age_min_fert = 15, # Closed [
                                              age_max_fert = 50, # Open )
                                              age_group = 5) #[,)

# Estimate age-specific mortality rates
asmr_sim <- rsocsim::estimate_mortality_rates(opop = opop, 
                                              final_sim_year = 2022, # [Jan-Dec]
                                              year_min = 1847, # Closed
                                              year_max = 2022, # Open )
                                              year_group = 5,
                                              age_max_mort = 110, # Open )
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
HFD <- read.table(file = "Input/NORasfrRR.txt", 
                  as.is=T, header=T, skip=2, stringsAsFactors=F)

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
  mutate(Source = "HFD/HMD", 
         Rate = "ASFR")

# Wrangle SOCSIM data
SocsimF <- asfr_sim %>% 
  rename(ASFR = socsim) %>% 
  mutate(Source = "SOCSIM",
         Rate = "ASFR")


# MORTALITY
# Extract year and age breaks used in the estimate_mortality_rates() to apply the same values to HMD data
# Year breaks. Extract all the unique numbers from the intervals 
year_breaks_mort <- unique(as.numeric(str_extract_all(asmr_sim$year, "\\d+", simplify = T)))

# Year range to filter HMD data
year_range_mort <- min(year_breaks_mort):max(year_breaks_mort-1)

# Age breaks of mortality rates. Extract all the unique numbers from the intervals 
age_breaks_mort <- unique(as.numeric(str_extract_all(asmr_sim$age, "\\d+", simplify = T)))


# Read the same HMD female life table 1x1 file used to write the input mortality files
ltf <- read.table(file="Input/fltper_1x1.txt",
                  as.is=T, header=T, skip=2, stringsAsFactors=F)

# Read the same HMD male life table 1x1 file used to write the input mortality files
ltm <- read.table(file="Input/mltper_1x1.txt",
                  as.is=T, header=T, skip=2, stringsAsFactors=F)

# Wrangle HMD life tables
HMD <- ltf %>%
  select(Year, Age, mx) %>%
  mutate(Sex = "Female") %>% 
  bind_rows(ltm %>% 
              select(Year, Age, mx) %>%  
              mutate(Sex = "Male")) %>% 
  mutate(Age = ifelse(Age == "110+", "110", Age),
         Age = as.numeric(Age)) %>% 
  filter(Year %in% year_range_mort) %>% 
  mutate(year = cut(Year, breaks = year_breaks_mort, 
                    include.lowest = F, right = F, ordered_results = T),
         age = cut(Age, breaks = age_breaks_mort, 
                   include.lowest = F, right = F, ordered_results = T)) %>%
  filter(!is.na(age)) %>% 
  group_by(year, Sex, age) %>% 
  summarise(mx = mean(mx)) %>%
  ungroup() %>% 
  mutate(Source = "HFD/HMD",
         Rate = "ASMR")

# Wrangle SOCSIM data
SocsimM <- asmr_sim %>% 
  rename(mx = socsim) %>% 
  mutate(Sex = ifelse(sex == "male", "Male", "Female"),
         Source = "SOCSIM",
         Rate = "ASMR") %>% 
  select(year, Sex, age,  mx, Source, Rate)



## PLOT SIMULATED AND INPUT DATA TOGETHER

# create folder to store graphs
graph.folder <- paste0(folder,"/sim_results_", supfile, "_",seed,"_/graphs/")# Check if the folder already exists
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}

## Plotting ASFR and ASMR (for females) from HFD/HMD vs SOCSIM 
yrs_plot <- c("[1847,1852)", "[1877,1882)", "[1907,1912)", "[1937,1942)", "[1967,1972)", "[1997,2002)", "[2017,2022)") 

# Get the age levels to define them before plotting and avoid wrong order
age_levels <- levels(asmr_sim$age)

png(file = paste0(graph.folder,"rates.png"),
    width = 964, height = 556)
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
  scale_linetype_manual(values = c("HFD/HMD" = "solid","SOCSIM" = "11")) +
  facetted_pos_scales(y = list(ASFR = scale_y_continuous(),
                               ASMR =  scale_y_continuous(trans = "log10")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title = paste0("Age-Specific Fertility and Mortality rates in Norway (1846-2022) 
       \n retrieved from HFD, HFC, HMD and a SOCSIM simulation 
       \n ", het," heterogeneous fertility, ", bint,  ", opop size = ", size_opop, 
                      "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"), 
       x = "Age") + 
  theme_bw()
dev.off()


## GETTING ESTIMATES -- IDENTIFY KIN!

#Use the final simulation year (January-December)
#and the last simulated month to convert our monthly dates into yearly ones
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

#Parameters specific to this simulation: will need to be changed
last_month <- max(opop$dob) # Last simulated month
final_sim_year <- 2022

#Cleaning our population file
opop <- opop %>% 
  #Fixing dates of death for individuals still living at the end
  mutate(dod = if_else(dod == 0, last_month, dod)) %>%
  #Dates of birth and death in years for both individual and parents
  mutate(dob_year = asYr(dob, last_month=last_month, final_sim_year=final_sim_year),
         dod_year = asYr(dod, last_month=last_month, final_sim_year=final_sim_year))


# Saving opop data frame for usage in R on my Mac
save(opop, file = paste0(folder,"/sim_results_", supfile, "_",seed,"_/opop.RData"))



#### PARENTS ####
## Build register containing death dates for parents for sample individuals (born 1953 and alive in 2019)

# load simulated register data
load(paste0(folder, "/sim_results_", supfile, "_", seed, "_/opop.RData"))

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

sample <- opop %>% 
  # keep if born in 1953
  filter(dob_year==1953)  %>% 
  # keep if alive in 2019 or later
  filter(dod_year>2019)

mothers <- sample %>% 
  # join to mother identifier (mom) in sample (left) data from opop (right) 
  left_join(opop, by = c("mom" = "pid"), suffix = c("", "_p"), keep = TRUE) %>% 
  # keep original variables plus parental death years 
  select(pid, mom, pop, pid_p, dod_year_p) %>% 
  mutate(mother = 1)

fathers <- sample %>% 
  # join to father identifier (mom) in sample (left) data from opop (right) 
  left_join(opop, by = c("pop" = "pid"), suffix = c("", "_p"), keep = TRUE) %>% 
  # keep original variables plus parental death years 
  select(pid, mom, pop, pid_p, dod_year_p) %>% 
  mutate(mother = 0)

## Build one parent register for ALL sample individuals
parents <- bind_rows(mothers, fathers)
# ascending order by pid
parents <- arrange(parents, pid)

## Select parent who died last
lateparent <- parents %>% 
  group_by(pid) %>% 
  slice_max(dod_year_p, with_ties = FALSE) %>% 
  ungroup %>% 
  mutate(dead_p = ifelse(dod_year_p > 2019, 0, 1))

#### CHILDREN ####
## Build register containing children of FEMALE sample individuals
children_f <- sample %>% 
  # select females
  filter(fem == 1) %>% 
  # merge sample vars (pid) to sample via sample individuals who became mothers 
  left_join(opop, by = c("pid" = "mom"), suffix = c("", "_c"), keep = TRUE) %>% 
  # keep original variables plus child birth and death years 
  select(pid,  
         pid_c, fem_c, mom_c, pop_c, dob_year_c, dod_year_c)
      
## Build register containing children of MALE sample individuals
children_p <- sample %>% 
  # select males
  filter(fem == 0) %>% 
  # merge sample vars (pid) to sample via sample individuals who became fathers 
  left_join(opop, by = c("pid" = "pop"), suffix = c("", "_c"), keep = TRUE) %>% 
  # keep original variables plus child birth and death years 
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
  # add number of all children per sample individual
  add_count(pid, name = "numkids") %>% 
  # select child that was born first
  slice_min(dob_year_c) %>%
  ungroup() %>% 
  mutate(numkids = ifelse(is.na(dob_year_c), NA, numkids)) 


#### GRANDCHILDREN ####
## Build register containing children of FEMALE children --> sample's grandchildren
gchildren_f <- children %>% 
  # select females
  filter(fem_c == 1) %>% 
  # merge sample vars (pid_c) to sample via children who became mothers 
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


## Build register containing children of MALE sample individuals
gchildren_p <- children %>% 
  # select females
  filter(fem_c == 0) %>% 
  # merge sample vars (pid_c) to sample via children who became fathers 
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
  mutate(numgkids = ifelse(is.na(dob_year_gc), NA, numgkids)) 

#### ONE FAMILY ROOSTER ####
gp_wide <- left_join(sample, lateparent)
gp_wide <- left_join(gp_wide, oldestc)
gp_wide <- left_join(gp_wide, oldestgc)

# ...with most basic info
gp_base <- gp_wide %>% 
  select(pid, dob_year, dod_year_p, dob_year_c, dob_year_gc) %>% 
  # right-censor at year 2019
  mutate(dod_year_p = ifelse(dod_year_p > 2019, NA, dod_year_p),
         dob_year_c = ifelse(dob_year_c > 2019, NA, dob_year_c),
         dob_year_gc = ifelse(dob_year_gc > 2019, NA, dob_year_gc)) %>% 
  # convert years into ages
  mutate(dage = dod_year_p - dob_year,
         cage = dob_year_c - dob_year,
         gcage = dob_year_gc - dob_year)



## Average age at family transitions
avg_par <- mean(gp_base$dage, na.rm = TRUE)
avg_child <- mean(gp_base$cage, na.rm = TRUE)
avg_gchild <- mean(gp_base$gcage, na.rm = TRUE)


## Same logic as in STATA

n_sample  <- nrow(sample)
max_age <- 68
ages <- as.character(c(1:max_age-1))
lab_ages <- paste0("age", ages)
empty_cols <- data.frame(matrix(NA, nrow = n_sample, ncol = max_age))
colnames(empty_cols) <- lab_ages


gp_long <- gp_base %>% 
  select(pid, dage, cage, gcage) %>% 
  # duplicate each sample individual times the age range we consider to create long format
  slice(rep(row_number(), each = max_age)) %>% 
  # generate age var
  group_by(pid) %>% 
  mutate(age = row_number()) %>% 
  
  # set ages after parents have died to 2 (1 = at least one parent; 2 = NO parent) 
  # if parental death not yet happened to 1
  mutate(gpp = ifelse(age >= dage & !is.na(dage), 2, 1),
         
         # set ages after child(ren) were born to 10 (10 = at least one child; 20 = NO child)
         # if childless and before to 20
         gpc = ifelse(age >= cage & !is.na(cage), 10, 20),
         
         # set ages after grandchild(ren) were born to 100 (100 = at least one grandchild; 200 = NO grandchild)
         # if grandchildless and before to 200
         gpg = ifelse(age >= gcage & !is.na(gcage), 100, 200),
         
         # generate "sum" of generational placements
         gp = gpp + gpc + gpg) %>% 
  ungroup

library(tidyr)
gp_wide <- gp_long %>%
  pivot_wider(
    id_cols = c(pid, dage, cage, gcage),
    names_from = age,
    values_from = gp,
    names_glue = "{.value}_{age}"
  )

# transform numerics into characters
# mapping dictionary (attention: logic between numerics and characters is inverted (P1xxxx or P0xxxx refers to xx1 or xx2))
mapping <- c("P1C0G0" = 221, "P0C0G0" = 222, 
             "P1C1G0" = 211, "P0C1G0" = 212, 
             "P1C1G1" = 111, "P0C1G1" = 112)

gp <- gp_wide %>% 
  mutate(across(starts_with("gp"), ~case_when(
    . == mapping[1] ~ "P1C0G0",
    . == mapping[2] ~ "P0C0G0",
    . == mapping[3] ~ "P1C1G0",
    . == mapping[4] ~ "P0C1G0",
    . == mapping[5] ~ "P1C1G1",
    . == mapping[6] ~ "P0C1G1",
    TRUE ~ as.character(.)
  )))

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
cblind <- brewer.pal(6, "Paired")

# Aggregate identical sequences to reduce computing time and memory used
# see Studer, 2013 (WeightedCluster Library Manual, Annex 1)

# ac <- wcAggregateCases(gp[, 5:72])
# ac

# Define sequence object
# with agg weights
# seq <- seqdef(gp[ac$aggIndex, 5:72], 
#               labels = gplabels,  
#               cnames = ages, 
#               tick.last = TRUE, 
#               xtstep = 5, 
#               cpal = cblind, 
#               alphabet = gpalpha, 
#               states = gpstates)

seq <- seqdef(gp, 5:72, 
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates)

# Plots
#seqIplot(seq, border = NA, ltext = c(gpstates), with.legend = "right")
png(file = paste0(graph.folder, "seqD_full.png"),
    width=964, height=556)
seqdplot(seq, border = NA, ltext = c(gpstates), with.legend = "right", 
         main = paste0("Simulated Data \n (1846 - 2022), HFC & HMD, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
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

# add weights ", weights = ac$aggWeights" here if using ac
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
                         # add ", weights = ac$aggWeights" if using ac
                         initialclust = chi_ward)

png(file = paste0(graph.folder, "clustqual.png"),
    width=964, height=556)
plot(chi_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2, 
     main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                   "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
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

# State distribution plots by clusters
# CHI2
png(file = paste0(graph.folder, "seqD_4.png"),
    width=964, height=556)
seqdplot(seq, group = chi_pam10$clustering$cluster4,
         border = NA, ltext = c(gpstates), 
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
dev.off()

png(file = paste0(graph.folder, "seqD_5.png"),
    width=964, height=556)
seqdplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),  
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
dev.off()

png(file = paste0(graph.folder, "seqD_6.png"),
    width=964, height=556)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
        border = NA, ltext = c(gpstates),  
        main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                      "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
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
by(seq, chi_pam10$clustering$cluster6, seqmeant)
png(file = paste0(graph.folder, "mean_plot.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster6, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 6 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
dev.off()

#### GRAPH ####
# We extract X clusters and re-label them from 1 to X to replace the medoid identifiers
# and attach the vector with the clusters to the main dataframe gp

# identify medoids sorted by frequency
# and disaggregate data

# add "[ac$disaggIndex]" behind first row to disaggregate
mc <- chi_pam10$clustering$cluster6
med <- as.data.frame(sort(table(mc), increasing = TRUE))
med1 <- as.character(med[1,1])
med2 <- as.character(med[2,1])
med3 <- as.character(med[3,1])
med4 <- as.character(med[4,1])
med5 <- as.character(med[5,1])
med6 <- as.character(med[6,1])
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5, med6),
                    as.character("1","2","3","4","5","6"),
                    labels = c("Cluster 1 -\n Childless, early parental death",
                               "Cluster 2 -\n Childless, late parental death",
                               "Cluster 3 -\n 2-gen family",
                               "Cluster 4 -\n 4-gen family)",
                               "Cluster 5 -\n Early 3-gen family",
                               "Cluster 6 -\n Late 3-gen family"))

gp$chi <- mc.factor

seqchi <- seqdef(gp, 5:72,
                 labels = gplabels,  
                 cnames = ages, 
                 tick.last = TRUE, 
                 xtstep = 5, 
                 cpal = cblind, 
                 alphabet = gpalpha, 
                 states = gpstates)

png(file = paste0(graph.folder, "seqD_6_lab.png"),
    width=964, height=556)
seqdplot(seqchi, group = group.p(gp$chi), border = NA,
         ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n alpha = ", alpha, ", beta = ", beta, " (", seed, ")"))
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
clusters <- left_join(gp, lateparent)
clusters <- left_join(clusters, oldestc)
clusters <- left_join(clusters, oldestgc)
clusters <- left_join(clusters, sample)

descr <- clusters %>% 
  select(pid, dage, cage, gcage, chi, dead_p, numkids, numgkids, fem, mstat) %>% 
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
save(descr, file = paste0(folder, "/sim_results_", supfile, "_", seed, "_/sumfile_",seed,".RData"))
dev.off()

# Descriptive table using sumtable

# SUMTABLE EXPORT KLAPPT HIER NOCH NICHT!!
st(descr, 
   vars = c("isparent", "numkids", "cage", 
            "isgparent", "numgkids", "gcage", 
            "dead_p", "dage",
            "fem", "mstat"),
   summ=c("mean(x)",
          "sd(x)"),
   summ.names = c("Mean",
                  "Std.Dev"),
   group = "chi", 
   group.test =TRUE,
   labels = c(paste0("Parent at age ", max_age), "Number of children", "Age at birth of first child",
              paste0("Grandparent at age ", max_age), "Number of grandchildren", "Age at birth of first grandchild",
              paste0("Both parents dead at age ", max_age), "Age at death of second parent", 
              "Female", "Marital status"),
   title = paste0("Summary Statistics (Simulated data, ", het, " heterogeneous fertility, ", bint,  ", opop size = ", seed, ")"),
   out = "csv",
   file = paste0(graph.folder, "sumtable"))

}
