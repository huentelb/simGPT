# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# COMPARISON OF FUTURE SCENARIOS 
# 2000 BIRTH COHORT - AGE RANGE 0-100

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse) # for rownames_to_column
library(janitor)

cohort <- 2000
max_age <- 100
ages <- as.character(c(0:max_age))
lab_ages <- paste0("age", ages)


# Generate folders to store results

# 1. Upper level folder based on simulation base_seed
folder.baseseed <- paste0(folder,"/sim_results_", base_seed,"_/")
if (!dir.exists(folder.baseseed)) {
  # If not, create the new folder
  dir.create(folder.baseseed)
  cat("Folder created:", folder.baseseed, "\n")
} else {
  cat("Folder already exists:", folder.baseseed, "\n")
}

# 2. Lower level folder for cohort-max age-specific output 
graph.folder <- paste0(folder.baseseed, cohort, "_future_scenarios/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}



# Laod gp data of both cohorts and store as separate dataframes with according suffix
load(paste0(folder.baseseed, "gp2000100.RData"))
gpmed <- gp %>% 
  mutate(scenario = 1)

load(paste0(folder, "/sim_results_260306_low/gp2000100.RData"))
gplo <- gp %>% 
  mutate(scenario = 0)

load(paste0(folder, "/sim_results_260305_high/gp2000100.RData"))
gphi <- gp %>% 
  mutate(scenario = 2)



### OVERALL ####

### AGGREGATE-LEVEL COMPARISONS ####

# Means
agg_mean <- gpmed %>% 
  rbind(gphi) %>%
  rbind(gplo) %>% 
  dplyr::select(scenario, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(scenario) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  dplyr::select(scenario, dage:gcage) # exchange 'cohort' for 'dob_year'

# Medians
agg_p50 <- gpmed %>% 
  rbind(gphi) %>%
  rbind(gplo) %>%  
  dplyr::select(scenario, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(scenario) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  dplyr::select(scenario, dage:gcage) %>% # exchange 'cohort' for 'dob_year'
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -scenario)

# SD
agg_sd <- gpmed %>% 
  rbind(gphi) %>%
  rbind(gplo) %>%  
  dplyr::select(scenario, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(scenario) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  dplyr::select(scenario, dage:gcage) %>% # exchange 'cohort' for 'dob_year'
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -scenario)

# Combine means, medians and SD into one df
agg <- agg_mean %>%
  left_join(agg_sd, by = "scenario") %>%
  left_join(agg_p50, by = "scenario") %>%
  dplyr::select (
    scenario,
    starts_with("dage"),
    starts_with("dead_p"),
    starts_with("pdage"),
    starts_with("isparent"),
    starts_with("numkids"),
    starts_with("cage"),
    starts_with("isgparent"),
    starts_with("numgkids"),
    starts_with("gcage")
  )


# Swap rows and columns of indic_mean
tab_a2_agg <- agg %>%
  mutate(across(-scenario, ~ round(.x, 2))) %>%
  pivot_longer(-scenario, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = scenario, values_from = value) 


# Number of observations
tab_a2_n <- gpmed %>% 
  rbind(gphi) %>%
  rbind(gplo) %>% 
  count(scenario) %>% 
  pivot_wider(names_from = scenario, values_from = n) %>% 
  as.data.frame

rownames(tab_a2_n) <- "N"

tab_a2_n <- rownames_to_column(tab_a2_n, var = "variable")




### DEFINE SEQUENCES ####
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
seqmed <- seqdef(gpmed, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
                labels = gplabels,  
                cnames = ages, 
                tick.last = TRUE, 
                xtstep = 5, 
                cpal = cblind, 
                alphabet = gpalpha, 
                states = gpstates,
                missing = "D", right = "DEL")

seqlo <- seqdef(gplo, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
                 labels = gplabels,  
                 cnames = ages, 
                 tick.last = TRUE, 
                 xtstep = 5, 
                 cpal = cblind, 
                 alphabet = gpalpha, 
                 states = gpstates,
                 missing = "D", right = "DEL")

seqhi <- seqdef(gphi, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
                 labels = gplabels,  
                 cnames = ages, 
                 tick.last = TRUE, 
                 xtstep = 5, 
                 cpal = cblind, 
                 alphabet = gpalpha, 
                 states = gpstates,
                 missing = "D", right = "DEL")



### COMPARISON OF SEQUENCES ####


#### State Distribution Plot #####

gp <- gplo %>% 
  rbind(gpmed, gphi)


### Define joint sequence including 'dead' as explicit state for state distribution plot

gpalpha_D <- c("P1C0G0", "P0C0G0", 
             "P1C1G0", "P0C1G0", 
             "P1C1G1", "P0C1G1", "D")

# Labels of states
gplabels_D <- c("child", "no ancestors/descentants",
              "child and parent", "parent",
              "child, parent, and grandparent", 
              "parent and grandparent", "dead")

# Short labels for states
gpstates_D <-c("G1G2", "G2",
             "G1G2G3", "G2G3",
             "G1G2G3G4", "G2G3G4", "Dead")

cblind_D <- c(cblind, "#f7f7f7")


seq <- seqdef(gp, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
                labels = gplabels_D,  
                cnames = ages, 
                tick.last = TRUE, 
                xtstep = 5, 
                cpal = cblind_D, 
                alphabet = gpalpha_D, 
                states = gpstates_D)

png(file = paste0(graph.folder, "seqD_scenarios.png"),
    width=800, height=375)
seqdplot(seq, group = gp$scenario, 
         ltext = c(gpstates_D), legend.prop = .15,
         use.layout = TRUE, with.legend = "auto",
         rows = 1, cols = 3,
         xlab = "Age",
         main = c("Low Fertility", "Medium Fertility (main)", "High Fertility"),
         border = NA,
         cex.lab = 1.3,
         cex.main = 1.6,
         cex.axis = 1.2)
dev.off()

pdf(file = paste0(graph.folder, "seqD_scenarios.pdf"),
    width=10, height=4.5)
seqdplot(seq, group = gp$scenario, 
         ltext = c(gpstates_D), legend.prop = .175,
         use.layout = TRUE, with.legend = "auto",
         rows = 1, cols = 3,
         xlab = "Age",
         main = c("Low Fertility", "Medium Fertility (main)", "High Fertility"),
         border = NA,
         cex.lab = 1.3,
         cex.main = 1.4,
         cex.axis = 1.2)
dev.off()


#### Sequence indicators ####

# "taking nonvisited states into account -> variance as a predictability indicator of the state duration"
# "ignoring nonvisited states -> measuring the variance of the observed spells." (Ritschard 2023: 2045)

# in our case: number of spells = number of visited states; whenever a new spell begins, a new state is visited
# as one cannot go back to a state that they left (parents are not reborn; in our setting, (grand)children do not die)
# -> recurrence meaningless

indicmed <- seqindic(seqmed, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index

indiclo <- seqindic(seqlo, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index

indichi <- seqindic(seqhi, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index

# store mean across full sample
indic_mean_med <- indicmed %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(scenario = 1) %>% 
  dplyr::select(-Visited) %>%
  dplyr::select(scenario, Lgth:Cplx)

indic_mean_lo <- indiclo %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(scenario = 0) %>% 
  dplyr::select(-Visited) %>%
  dplyr::select(scenario, Lgth:Cplx)

indic_mean_hi <- indichi %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(scenario = 2) %>% 
  dplyr::select(-Visited) %>%
  dplyr::select(scenario, Lgth:Cplx)

indic_mean <- indic_mean_lo %>% 
  rbind(indic_mean_med) %>%  
  rbind(indic_mean_hi)


# Swap rows and columns of indic_mean
tab_a2_ind <- indic_mean %>%
  pivot_longer(-scenario, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = scenario, values_from = value) 



#### BIC & LRT ####

# Compare occurrence of states
comp_occ_01 <- round(seqCompare(seqlo, seqdata2 = seqmed, stat = "all", 
                             method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5),
                  2)

comp_occ_02 <- round(seqCompare(seqlo, seqdata2 = seqhi, stat = "all", 
                                method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5),
                     2)

comp_occ_12 <- round(seqCompare(seqmed, seqdata2 = seqhi, stat = "all", 
                                method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5),
                     2)

# Compare timing of statesgp$
comp_time_01 <- round(seqCompare(seqlo, seqdata2 = seqmed, stat = "all", 
                              method = "CHI2", step = 101),
                   2)

comp_time_02 <- round(seqCompare(seqlo, seqdata2 = seqhi, stat = "all", 
                                 method = "CHI2", step = 101),
                      2)

comp_time_12 <- round(seqCompare(seqmed, seqdata2 = seqhi, stat = "all", 
                                 method = "CHI2", step = 101),
                      2)

# Compare duration in statesgp$
comp_dur_01 <- round(seqCompare(seqlo, seqdata2 = seqmed, stat = "all", 
                             method = "OMstran", otto = 0.5, sm = "INDELSLOG"),
                  2)

comp_dur_02 <- round(seqCompare(seqlo, seqdata2 = seqhi, stat = "all", 
                                method = "OMstran", otto = 0.5, sm = "INDELSLOG"),
                     2)

comp_dur_12 <- round(seqCompare(seqmed, seqdata2 = seqhi, stat = "all", 
                                method = "OMstran", otto = 0.5, sm = "INDELSLOG"),
                     2)


# Combine all BIC and LRT into one df
bic <- as.data.frame(comp_occ_01, row.names = "Occurence (low vs med)") %>% 
  rbind(as.data.frame(comp_time_01, row.names = "Timing (low vs med)"),
        as.data.frame(comp_dur_01, row.names = "Duration (low vs med)"),
        as.data.frame(comp_occ_02, row.names = "Occurence (low vs high)"),
        as.data.frame(comp_time_02, row.names = "Timing (low vs high)"),
        as.data.frame(comp_dur_02, row.names = "Duration (low vs high)"),
        as.data.frame(comp_occ_12, row.names = "Occurence (med vs high)"), 
        as.data.frame(comp_time_12, row.names = "Timing (med vs high)"),
        as.data.frame(comp_dur_12, row.names = "Duration (med vs high)")) %>% 
  select(3,1,2) # select and order columns (BIC, LRT, p-value) 

bic <- rownames_to_column(bic, var = "Measure")



#### AVERAGE GP MEAN TIME ####
gp_mt_lo <- round(seqmeant(seqlo),2) #, serr = T for calculating Var, SD, SE
gp_mt_med <- round(seqmeant(seqmed),2)
gp_mt_hi <- round(seqmeant(seqhi),2)


colnames(gp_mt_lo) <- 0
colnames(gp_mt_med) <- 1
colnames(gp_mt_hi ) <- 2

tab_a2_mt <- gp_mt_lo %>% 
  cbind(gp_mt_med) %>% 
  cbind(gp_mt_hi) %>% 
  as.data.frame() 

tab_a2_mt <- rownames_to_column(tab_a2_mt, var = "variable")



#### CREATE TAB A2 ####
tab_a2 <- tab_a2_agg %>% 
  rbind(tab_a2_n, tab_a2_mt, tab_a2_ind) %>%
  as.data.frame()


library(gt)
library(flextable)

set_flextable_defaults(
  font.size = 11,
  border.color = 'black',
  line_spacing = 1.3,
)

ft_a2 <- flextable(tab_a2) %>% 
  add_header_row(colwidths = c(1,3), values = c(" ", "Scenario")) %>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft_a2

save_as_docx("Table A2" = ft_a2, path = paste0(graph.folder, "Tab_A2.docx"), align = "left")

ft_bic <- flextable(bic)%>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft_bic

save_as_docx("BIC_A2" = ft_bic, path = paste0(graph.folder, "BIC_A2.docx"), align = "left")

