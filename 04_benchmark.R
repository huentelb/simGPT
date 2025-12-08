# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# BENCHMARK MICROSIMULATED AGAINST EMPIRICAL REGISTER DATA
# 1960 BIRTH COHORT - AGE RANGE 0-67 

# ONLY RUNS ON SERVER WITH NORWEGIAN REGISTER DATA

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### SET-UP #####

library(tidyverse) 
library(dplyr)
library(foreign)
library(gtable)

cohort <- 1960
max_age <- 59
ages <- as.character(c(0:max_age))



# Generate folders to store results on Norwegian server

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

##### LOAD MICROSIMULATED DATA ####


# You need to manually upload the synthetic gp data to the Norwegian server 
# Load the data from the folder you stored it in

# load the gp-data for from the simulation
load("N:\\durable\\Data20\\Project_JW_FamWealth\\simgpt\\data\\gp196059.RData")

# store as gp-dataframe with suffix m (microsimulation)
# and re-order and select columns to match dataframe of empirical data
gpm <- gp %>% 
  select(pid, gp_0:gp_59, dage:gcage, numkids, numgkids, dead_p, isparent, isgparent) %>% 
  rename_with(~ gsub("^gp_(\\d+)$", "gp\\1", .), starts_with("gp")) %>% # remove _ 
  mutate(source = "Microsimulation",
         pid = paste0("P", pid)) # add suffix P to pid 


##### LOAD EMPIRICAL REGISTER DATA #####

gp <- read.dta("N:\\durable\\Data20\\Project_JW_FamWealth\\simgpt\\data\\gp_wide.dta", 
               convert.factors = FALSE)

# store as gp-dataframe with suffix e (empirical)
gpe <- gp %>% 
  mutate(pid = id) %>% 
  filter(dage > 0 | is.na(dage)) %>% # exclude individuals who did not survive until first birthday
  select(pid, gp0:gp59, dage, pdage, cage, gcage, numkids, numgkids, dead_p,
         isparent, isgparent) %>% 
  mutate(source = "Register")




### AGGREGATE-LEVEL COMPARISONS ####

# Means
agg_mean <- gpe %>% 
  rbind(gpm) %>% 
  select(source, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(source) %>% 
  summarise_all(mean, na.rm = TRUE)
  

# Medians (for cont vars)
agg_p50 <- gpe %>% 
  rbind(gpm) %>% 
  select(source, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(source) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -source) 

# SD (for cont vars)
agg_sd <- gpe %>% 
  rbind(gpm) %>% 
  select(source, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(source) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -source) 


# Combine means, medians and SD into one df
agg <- agg_mean %>% 
  left_join(agg_sd, by = "source") %>% 
  left_join(agg_p50, by = "source") %>% 
  select(
    source,
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
tab1_agg <- agg %>%
  mutate(across(-source, ~ round(.x, 2))) %>%
  pivot_longer(-source, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = source, values_from = value) %>%
  column_to_rownames("variable") # Moves "variable" column to row names




# Add significance test for difference
gp_boot <- gp60 %>% 
  rbind(gp00) %>% 
  dplyr::select(dob_year, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  mutate(cohort = dob_year) 

# bootstrapped means
library(tidyverse)
library(mosaic)


b <- 1000 # number of samples for bootstrapping

# Bootstrapping MEANS
for (i in (c("dage", "dead_p", "pdage", "isparent", "numkids", "cage", "isgparent", "numgkids", "gcage"))) {
  
  set.seed(123456) # put into loop to get same samples across all indicators
  formula <- as.formula(paste(i, "~ cohort"))
  
  gp_boot_mean <- do(b)*mean(formula, na.rm = TRUE, data = mosaic::resample(gp_boot))
  gp_boot_mean <- gp_boot_mean %>%
    mutate(diff = X1960-X2000)
  
  # add CI to Table 1
  tab1_agg[i,3] <- round(tab1_agg[i,1]-tab1_agg[i,2],2) #mean(gp_boot_mean$diff)
  tab1_agg[i,4:5] <- round(confint(gp_boot_mean$diff, level = 0.95),2) #same as quantile(gp_boot_mean$diff, c(0.025, 0.975))
  
  
}

# Bootstrapping SD
for (i in (c("dage", "pdage", "numkids", "cage", "numgkids", "gcage"))) {
  
  set.seed(123456)
  formula <- as.formula(paste(i, "~ cohort"))
  
  gp_boot_sd <- do(b)*sd(formula, na.rm = TRUE, data = mosaic::resample(gp_boot))
  gp_boot_sd <- gp_boot_sd %>%
    mutate(diff = X1960-X2000)
  
  # add CI to Table 1
  tab1_agg[paste0(i,"_sd"),3] <- round(tab1_agg[paste0(i,"_sd"),1]-tab1_agg[paste0(i,"_sd"),2],2) #mean(gp_boot_mean$diff)
  tab1_agg[paste0(i,"_sd"),4:5] <- round(confint(gp_boot_sd$diff, level = 0.95),2)
  
}

### HERE IS SOMETHING STILL OFF! ####
# Bootstrapping MEDIANS
for (i in (c("dage", "pdage", "numkids", "cage", "numgkids", "gcage"))) {
  
  set.seed(123456) 
  formula <- as.formula(paste(i, "~ cohort"))
  
  gp_boot_p50 <- do(b)*median(formula, na.rm = TRUE, data = mosaic::resample(gp_boot))
  gp_boot_p50 <- gp_boot_p50 %>%
    mutate(diff = X1960-X2000)
  
  # add CI to Table 1
  tab1_agg[paste0(i,"_p50"),3] <- round(tab1_agg[paste0(i,"_p50"),1]-tab1_agg[paste0(i,"_p50"),2],2) #mean(gp_boot_mean$diff)
  tab1_agg[paste0(i,"_p50"),4:5] <- round(confint(gp_boot_p50$diff, level = 0.95),2)
  
}



# Number of observations
tab1_n <- gpe %>%
  rbind(gpm) %>% 
  count(source) %>% 
  pivot_wider(names_from = source, values_from = n) %>% 
  as.data.frame()

rownames(tab1_n) <- "N"



#### SEQUENCE ANALYSIS ####

library(foreign)
library(TraMineR)
library(TraMineRextras)
library(WeightedCluster)
library(RColorBrewer)
library(bookdown)


### DEFINE SEQUENCES ####
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
# (without ac weights for seq indicators and BIC/LRT comparisons)

seqm <- seqdef(gpm, 2:paste0(max_age+2),
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")


seqe <- seqdef(gpe, 2:paste0(max_age+2),
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates,
              missing = "D", right = "DEL")



### COMPARISON OF SEQUENCES ####

##### Sequence indicators ####

# "taking nonvisited states into account -> variance as a predictability indicator of the state duration"
# "ignoring nonvisited states -> measuring the variance of the observed spells." (Ritschard 2023: 2045)

# in our case: number of spells = number of visited states; whenever a new spell begins, a new state is visited
# as one cannot go back to a state that they left (parents are not reborn; in our setting, (grand)children do not die)
# -> recurrence meaningless

indicm <- seqindic(seqm, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index

indice <- seqindic(seqe, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index


# store mean across full sample
indic_meane <- indice %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(source = "Register") %>% 
  select(-Visited) %>% 
  select(source, Lgth:Cplx)

indic_meanm <- indicm %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(source = "Microsimulation") %>% 
  select(-Visited) %>% 
  select(source, Lgth:Cplx)

indic_mean <- indic_meane %>% 
  rbind(indic_meanm) 



# Swap rows and columns of indic_mean
tab1_ind <- indic_mean %>%
  pivot_longer(-source, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = source, values_from = value) %>%
  column_to_rownames("variable")



##### BIC & LRT ####

# Compare occurrence of states
comp_occ <- round(seqCompare(seqm, seqdata2 = seqe, stat = "all", 
                             method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5),
                  2)

# Compare timing of statesgp$
comp_time <- round(seqCompare(seqm, seqdata2 = seqe, stat = "all", 
                              method = "CHI2", step = 101),
                   2)

# Compare duration in statesgp$
comp_dur <- round(seqCompare(seqm, seqdata2 = seqe, stat = "all", 
                             method = "OMstran", otto = 0.5, sm = "INDELSLOG"),
                  2)

# Combine all BIC and LRT into one df
bic <- as.data.frame(comp_occ, row.names = "Occurence") %>% 
  rbind(as.data.frame(comp_time, row.names = "Timing"),
        as.data.frame(comp_dur, row.names = "Duration")) %>% 
  select(3,1,2) # select and order columns (BIC, LRT, p-value) 

write.csv(bic, paste0(graph.folder, "tab1_bic.csv"))



##### AVERAGE GP MEAN TIME ####
gp_mte <- round(seqmeant(seqm),2) #, serr = T for calculating Var, SD, SE
gp_mtm <- round(seqmeant(seqe),2)

colnames(gp_mte) <- "Register"
colnames(gp_mtm) <- "Microsimulation"

tab1_mt <- gp_mte %>% 
  cbind(gp_mtm) %>% 
  as.data.frame() 


##### CREATE TAB 1 ####
tab1 <- tab1_agg %>% 
  rbind(tab1_mt, tab1_ind, tab1_n) %>% 
  select(Register, Microsimulation)

tab1 <- rownames_to_column(tab1, var = "Measure") 

write.csv(tab1, paste0(graph.folder, "tab1.csv"))





### Define sequence objects using AC weights ####
acm <- wcAggregateCases(gpm[, 2:paste0(max_age+2)])
acm

seqm <- seqdef(gpm[acm$aggIndex, 2:paste0(max_age+2)],
               weights = acm$aggWeights,
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")


ace <- wcAggregateCases(gpe[, 2:paste0(max_age+2)])
ace

seqe <- seqdef(gpe[ace$aggIndex, 2:paste0(max_age+2)],
               weights = ace$aggWeights,
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")


#### CLUSTER ANALYSIS EMPIRICAL ####

help.folder <- graph.folder

graph.folder <- paste0(help.folder, "cluster_emp/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}


# 1. Generate dissimilarity matrix
seq <- seqe
ac <- ace

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

chi_ward10 <- as.clustrange(chi_ward, diss = chi, ncluster = 10, weights = ac$aggWeights)
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
                         initialclust = chi_ward,
                         weights = ac$aggWeights)
saveRDS(chi_pam10, file = paste0(graph.folder, "chipam10.RData"))

chi_pam10 <- readRDS(paste0(graph.folder, "chipam10.RData"))




# Plot cluster quality
png(file = paste0(graph.folder, "clustqual.png"),
    width=964, height=556)
plot(chi_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2, 
     main = paste0("Empricial data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                   "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"))
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



##### D plots ####
# State distribution plots by clusters
# CHI2

w = 2000
h = 1250

png(file = paste0(graph.folder, "seqD_4.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster4, 
         border = NA, ltext = c(gpstates), with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_5.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_6.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_7.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

##### F plots ####
png(file = paste0(graph.folder, "seqF20_4.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster4,
         border = NA, ltext = c(gpstates),  with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_5.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_6.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_7.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n emp., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
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
by(seq, chi_pam10$clustering$cluster6, seqmeant)
png(file = paste0(graph.folder, "mean_plot_6.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster6, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 6 Clusters emp. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "mean_plot_5.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster5, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 5 Clusters emp. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7")
dev.off()




##### LABELLED GRAPH FOR OPTIMAL CLUSTER SOLUTION ####

# We extract X clusters and re-label them from 1 to X to replace the medoid identifiers

# identify medoids sorted by frequency
mc <- chi_pam10$clustering$cluster6[ac$disaggIndex]
med <- as.data.frame(sort(table(mc), decreasing = TRUE))
med1 <- as.character(med[1,1])
med2 <- as.character(med[2,1])
med3 <- as.character(med[3,1])
med4 <- as.character(med[4,1])
med5 <- as.character(med[5,1])
med6 <- as.character(med[6,1])

# store size of clusters for each cluster to add to titles
propmed <- as.data.frame(sort(prop.table(table(mc)), decreasing = TRUE))
propmed1 <- round(propmed[1,2], digits = 2)*100
propmed2 <- round(propmed[2,2], digits = 2)*100
propmed3 <- round(propmed[3,2], digits = 2)*100
propmed4 <- round(propmed[4,2], digits = 2)*100
propmed5 <- round(propmed[5,2], digits = 2)*100
propmed6 <- round(propmed[6,2], digits = 2)*100


# create factor containing medoids incl labels
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5, med6),
                    as.character(c("1","2","3","4","5","6")))


# store labels as values for later use
l1 <- as.character(paste0("Cluster 1 -\n Late 3-gen family (", propmed1, "%)"))
l2 <- as.character(paste0("Cluster 2 -\n 4-gen family (", propmed2, "%)"))
l3 <- as.character(paste0("Cluster 3 -\n 2-gen family (", propmed3, "%)"))
l4 <- as.character(paste0("Cluster 4 -\n Non-parent late (", propmed4, "%)"))
l5 <- as.character(paste0("Cluster 5 -\n Early 3-gen (", propmed5, "%)"))
l6 <- as.character(paste0("Cluster 6 -\n Non-parent early (", propmed6, "%)"))

# attach to dataframe to use as weights in plots
gpe$chi <- factor(mc.factor,
                 labels = c(l1,
                            l2,
                            l3,
                            l4,
                            l5,
                            l6))

# generate new sequence object (without aggregation)
seqchi <- seqdef(gpe, 2:paste0(max_age+2),
                 labels = gplabels,
                 cnames = ages,
                 tick.last = TRUE,
                 xtstep = 5,
                 cpal = cblind,
                 alphabet = gpalpha,
                 states = gpstates,
                 missing = "D", right = "DEL")

# different plots with labels
png(file = paste0(graph.folder, "seqD_6_lab.png"),
    width=w, height=h)
seqdplot(seqchi, group = gpe$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2)
dev.off()

png(file = paste0(graph.folder, "seqI_6_lab.png"),
    width=w, height=h)
seqIplot(seqchi, group = gpe$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqF100_6_lab.png"),
    width=w, height=h)
seqfplot(seqchi, group = gpe$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7", idxs = 1:100)
dev.off()

by(seqchi, gpe$chi, seqmeant)
png(file = paste0(graph.folder, "mean_plot_6_lab.png"),
    width=w, height=h)
seqmtplot(seqchi, group = gpe$chi, border = NA,
          ltext = c(gpstates), 
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "seqr_6_lab.png"),
    width=w, height=h)
seqrplot(seqchi, group = gpe$chi, border = NA,
         ltext = c(gpstates), 
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




##### Relative frequency plot (empirical) ####

w <- 1000
h <- 625

# Check different parameters against theoretically most reasonable: CHI2
# Sorting on the first MDS factor extracted from a dissimilarity matrix built using the CHI-square

# Sample random n=5,000 cases & define sequence object (analysis does not work with full dataset)
set.seed(2407)
testgp <- sample_n(gpe, 5000)
testseq <- seqdef(testgp, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
                  labels = gplabels,  
                  cnames = ages, 
                  tick.last = TRUE, 
                  xtstep = 5, 
                  cpal = cblind, 
                  alphabet = gpalpha, 
                  states = gpstates,
                  missing = "D", right = "DEL")


# 1) CHI2 distance
testchi <- seqdist(testseq, method = "CHI2", step = max(seqlength(testseq)))

# Select medoids based on distance
srfchi <- seqrf(testseq,
                diss = testchi,
                sortv = "mds",
                grp.meth = "first")

# RF plot: 
# Plot all k = 100 medoids + average distance of repr. sequences to medoid
png(file = paste0(graph.folder, "seqrf.png"),
    width=w, height=h)
plot(srfchi, which.plot = "both", main = paste0("Random sample from ", cohort, " birth cohort (n = 5,000)")) 
dev.off()
# summary(srfchi)

# For sequence index plot in order of RF plot: 
# a. Assign representing medoid to each sequence (medoid_id)
testgp <- testgp %>% 
  mutate(medoid_id = srfchi[["rf"]][["kmedoid.index"]])

# b. Sort medoid_id according to order of rfplot
testgp <- testgp %>% 
  mutate(medoid_id = factor(medoid_id, levels = srfchi[["rf"]][["medoids"]]))


# c. seqIplot sorted by medoid_id
seqIplot(testseq, border = NA,
         ltext = c(gpstates), 
         missing.color = "#f7f7f7", with.legend = FALSE,
         sortv = testgp$medoid_id)


##### RFplot by clusters (empirical) ####
w <- 7
h <- 7

# generate one RF plot per cluster
# store gp dataframes per cluster
c1 <- gpe %>% 
  filter(chi == l1)

c2 <- gpe %>% 
  filter(chi == l2)

c3 <- gpe %>% 
  filter(chi == l3)

c4 <- gpe %>% 
  filter(chi == l4)

c5 <- gpe %>% 
  filter(chi == l5)

c6 <- gpe %>% 
  filter(chi == l6)

# Cluster 1
seq1 <- seqdef(c1, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi1 <- seqdist(seq1, method = "CHI2", step = max(seqlength(seq1)))

# Select medoids based on distance
srfchi1 <- seqrf(seq1,
                 diss = chi1,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c1.pdf"),
    width=w, height=h)
plot(srfchi1, which.plot = "both", main = l1)
dev.off()

# Cluster 2
seq2 <- seqdef(c2, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi2 <- seqdist(seq2, method = "CHI2", step = max(seqlength(seq2)))

# Select medoids based on distance
srfchi2 <- seqrf(seq2,
                 diss = chi2,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c2.pdf"),
    width=w, height=h)
plot(srfchi2, which.plot = "both", main = l2)
dev.off()


# Cluster 3
seq3 <- seqdef(c3, 2:paste0(max_age+2), # for max_age 300 to column 306, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi3 <- seqdist(seq3, method = "CHI2", step = max(seqlength(seq3)))

# Select medoids based on distance
srfchi3 <- seqrf(seq3,
                 diss = chi3,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c3.pdf"),
    width=w, height=h)
plot(srfchi3, which.plot = "both", main = l3)
dev.off()


# Cluster 4
seq4 <- seqdef(c4, 2:paste0(max_age+2), # for max_age 400 to column 406, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi4 <- seqdist(seq4, method = "CHI2", step = max(seqlength(seq4)))

# Select medoids based on distance
srfchi4 <- seqrf(seq4,
                 diss = chi4,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c4.pdf"),
    width=w, height=h)
plot(srfchi4, which.plot = "both", main = l4)
dev.off()

# Cluster 5
seq5 <- seqdef(c5, 2:paste0(max_age+2), # for max_age 500 to column 506, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi5 <- seqdist(seq5, method = "CHI2", step = max(seqlength(seq5)))

# Select medoids based on distance
srfchi5 <- seqrf(seq5,
                 diss = chi5,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c5.pdf"),
    width=w, height=h)
plot(srfchi5, which.plot = "both", main = l5)
dev.off()


# Cluster 6
seq6 <- seqdef(c6, 2:paste0(max_age+2), # for max_age 100 to column 101, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi6 <- seqdist(seq6, method = "CHI2", step = max(seqlength(seq6)))

# Select medoids based on distance
srfchi6 <- seqrf(seq6,
                 diss = chi6,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c6.pdf"),
    width=w, height=h)
plot(srfchi6, which.plot = "both", main = l6)
dev.off()


# Combine all per-cluster rfplots into one graph
w <- 750
h <- 600

# png(file = paste0(graph.folder, "seqrf_cluster6.png"),
#     width=w, height=h)
pdf(paste0(graph.folder, "seqrf_cluster6.pdf"), width = 8, height = 6)  # Open PDF device
original_par <- par(no.readonly = TRUE) # store original current parameter
par(mfrow = c(3, 2), # 3 rows, 2 columns
    mar = c(3.5, 2, 3 , 2), # margins of each plot
    mgp = c(2, 1, 0)) # margins around axis title, axis labels, and axis line
plot(srfchi1, which.plot = "medoids", skipar = TRUE, main = l1, cex.main = 1, info = "none")
plot(srfchi2, which.plot = "medoids", skipar = TRUE, main = l2, cex.main = 1, info = "none")
plot(srfchi3, which.plot = "medoids", skipar = TRUE, main = l3, cex.main = 1, info = "none")
plot(srfchi4, which.plot = "medoids", skipar = TRUE, main = l4, cex.main = 1, info = "none")
plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1, info = "none", xlab = "Age")
plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1, info = "none", xlab = "Age")
dev.off()
par(original_par) # reset layout


pdf(paste0(graph.folder, "seqrf_both_cluster6.pdf"), 
    width = 8, height = 9)
par(mfrow = c(3, 4), # 3 rows, 4 columns
    mar = c(3.5, 2, 3 , 2), # margins of each plot
    mgp = c(2, 1, 0)) # margins around axis title, axis labels, and axis line
plot(srfchi1, which.plot = "medoids", skipar = TRUE, main = l1, cex.main = 1.1, info = "none")
plot(srfchi1, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi2, which.plot = "medoids", skipar = TRUE, main = l2, cex.main = 1.1, info = "none")
plot(srfchi2, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi3, which.plot = "medoids", skipar = TRUE, main = l3, cex.main = 1.1, info = "none")
plot(srfchi3, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi4, which.plot = "medoids", skipar = TRUE, main = l4, cex.main = 1.1, info = "none")
plot(srfchi4, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi5, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi6, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)
dev.off()
par(original_par) # reset layout



# 2) OM distance with transition rate based costs
omt <- seqdist(seq, method = "OM", indel = 1, sm = "TRATE")

# Select medoids based on distance
srfomt <- seqrf(seq,
                diss = omt,
                sortv = "mds",
                grp.meth = "prop",
                weights = ac$aggWeights)

# RF plot: 
# Plot all k = 100 medoids + average distance of repr. sequences to medoid
png(file = paste0(graph.folder, "seqrf_omt_e.png"),
    width=w, height=h)
plot(srfomt, which.plot = "both")
dev.off()
summary(srfomt)

# For sequence index plot in order of RF plot: 
# a. Assign representing medoid to each sequence (medoid_id)
testgp <- testgp %>% 
  mutate(medoid_id = srfomt[["rf"]][["kmedoid.index"]])

# b. Sort medoid_id according to order of rfplot
testgp <- testgp %>% 
  mutate(medoid_id = factor(medoid_id, levels = srfomt[["rf"]][["medoids"]]))

# c. seqIplot sorted by medoid_id
seqIplot(testseq, border = NA,
         ltext = c(gpstates), 
         missing.color = "#f7f7f7", with.legend = FALSE,
         sortv = testgp$medoid_id)


# RFplot by clusters
seqrfplot(testseq, group = group.p(testgp$chi), 
          diss = omt,
          sortv = "mds", 
          ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
          ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE)

# seqIplot by clusters sorted by medoid_id
seqIplot(testseq, group = group.p(testgp$chi),
         ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
         ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE,
         sortv = testgp$medoid_id)



##### Complexity by cluster (empirical) ####

indic1 <- seqindic(seq1, indic=c("cplx"), with.missing=F) # complexity index
indic2 <- seqindic(seq2, indic=c("cplx"), with.missing=F) 
indic3 <- seqindic(seq3, indic=c("cplx"), with.missing=F) 
indic4 <- seqindic(seq4, indic=c("cplx"), with.missing=F) 
indic5 <- seqindic(seq5, indic=c("cplx"), with.missing=F) 
indic6 <- seqindic(seq6, indic=c("cplx"), with.missing=F) 


# store means across full sample and rowbind into one dataframe
indic_mean_ec <- indic1 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  rbind(indic2 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic3 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic4 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic5 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic6 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)))

# Swap rows and columns of indic_mean
tab_ind_clusters_e <- as.data.frame(t(indic_mean_ec))


# Store as docx
knitr::kable(tab_ind_clusters_e, "simple")



#### CLUSTER ANALYSIS MICROSIMULATED ####


graph.folder <- paste0(help.folder, "cluster_sim/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}

# 1. Generate dissimilarity matrix
seq <- seqm
ac <- acm

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

chi_ward10 <- as.clustrange(chi_ward, diss = chi, ncluster = 10, weights = ac$aggWeights)
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
                         initialclust = chi_ward,
                         weights = ac$aggWeights)
saveRDS(chi_pam10, file = paste0(graph.folder, "chipam10.RData"))

chi_pam10 <- readRDS(paste0(graph.folder, "chipam10.RData"))




# Plot cluster quality
png(file = paste0(graph.folder, "clustqual.png"),
    width=964, height=556)
plot(chi_pam10, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore", lwd = 2, 
     main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                   "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"))
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



##### D plots ####
# State distribution plots by clusters
# CHI2

w = 2000
h = 1250

png(file = paste0(graph.folder, "seqD_4.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster4, 
         border = NA, ltext = c(gpstates), with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_5.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_6.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_7.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

##### F plots ####
png(file = paste0(graph.folder, "seqF20_4.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster4,
         border = NA, ltext = c(gpstates),  with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_5.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_6.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "seqF20_7.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
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
by(seq, chi_pam10$clustering$cluster6, seqmeant)
png(file = paste0(graph.folder, "mean_plot_6.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster6, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 6 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "mean_plot_5.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster5, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 5 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7")
dev.off()




##### LABELLED GRAPH FOR OPTIMAL CLUSTER SOLUTION ####

# We extract X clusters and re-label them from 1 to X to replace the medoid identifiers

# identify medoids sorted by frequency
mc <- chi_pam10$clustering$cluster6[ac$disaggIndex]
med <- as.data.frame(sort(table(mc), decreasing = TRUE))
med1 <- as.character(med[1,1])
med2 <- as.character(med[2,1])
med3 <- as.character(med[3,1])
med4 <- as.character(med[5,1]) # swap C4 and C5 to match order or emp data
med5 <- as.character(med[4,1])
med6 <- as.character(med[6,1])

# store size of clusters for each cluster to add to titles
propmed <- as.data.frame(sort(prop.table(table(mc)), decreasing = TRUE))
propmed1 <- round(propmed[1,2], digits = 2)*100
propmed2 <- round(propmed[2,2], digits = 2)*100
propmed3 <- round(propmed[3,2], digits = 2)*100
propmed4 <- round(propmed[5,2], digits = 2)*100 # swap C4 and C5 to match order or emp data
propmed5 <- round(propmed[4,2], digits = 2)*100
propmed6 <- round(propmed[6,2], digits = 2)*100


# create factor containing medoids incl labels
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5, med6),
                    as.character(c("1","2","3","4","5","6")))


# store labels as values for later use
l1 <- as.character(paste0("Cluster 1 -\n Late 3-gen family (", propmed1, "%)"))
l2 <- as.character(paste0("Cluster 2 -\n 4-gen family (", propmed2, "%)"))
l3 <- as.character(paste0("Cluster 3 -\n 2-gen family (", propmed3, "%)"))
l4 <- as.character(paste0("Cluster 4 -\n Non-parent late (", propmed4, "%)"))
l5 <- as.character(paste0("Cluster 5 -\n Early 3-gen (", propmed5, "%)"))
l6 <- as.character(paste0("Cluster 6 -\n Non-parent early (", propmed6, "%)"))


# attach to dataframe to use as weights in plots
gpm$chi <- factor(mc.factor,
                  labels = c(l1,
                             l2,
                             l3,
                             l4,
                             l5,
                             l6))

# generate new sequence object (without aggregation)
seqchi <- seqdef(gpm, 2:paste0(max_age+2),
                 labels = gplabels,
                 cnames = ages,
                 tick.last = TRUE,
                 xtstep = 5,
                 cpal = cblind,
                 alphabet = gpalpha,
                 states = gpstates,
                 missing = "D", right = "DEL")

# different plots with labels
png(file = paste0(graph.folder, "seqD_6_lab_m.png"),
    width=w, height=h)
seqdplot(seqchi, group = gpm$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2)
dev.off()

png(file = paste0(graph.folder, "seqI_6_lab_m.png"),
    width=w, height=h)
seqIplot(seqchi, group = gpm$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqf20_6_lab_m.png"),
    width=w, height=h)
seqfplot(seqchi, group = gpm$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

by(seqchi, gpm$chi, seqmeant)
png(file = paste0(graph.folder, "mean_plot_6_lab_m.png"),
    width=w, height=h)
seqmtplot(seqchi, group = gpm$chi, border = NA,
          ltext = c(gpstates), 
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()



##### Relative frequency plot (microsim) ####

w <- 1000
h <- 625

# Check different parameters against theoretically most reasonable: CHI2
# Sorting on the first MDS factor extracted from a dissimilarity matrix built using the CHI-square

# Sample random n=5,000 cases & define sequence object (analysis does not work with full dataset)
set.seed(2407)
testgp <- sample_n(gpm, 5000)
testseq <- seqdef(testgp, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
                  labels = gplabels,  
                  cnames = ages, 
                  tick.last = TRUE, 
                  xtstep = 5, 
                  cpal = cblind, 
                  alphabet = gpalpha, 
                  states = gpstates,
                  missing = "D", right = "DEL")


# 1) CHI2 distance
testchi <- seqdist(testseq, method = "CHI2", step = max(seqlength(testseq)))

# Select medoids based on distance
srfchi <- seqrf(testseq,
                diss = testchi,
                sortv = "mds",
                grp.meth = "first")

# RF plot: 
# Plot all k = 100 medoids + average distance of repr. sequences to medoid
png(file = paste0(graph.folder, "seqrf.png"),
    width=w, height=h)
plot(srfchi, which.plot = "both", main = paste0("Random sample from ", cohort, " birth cohort (n = 5,000)")) 
dev.off()
# summary(srfchi)

# For sequence index plot in order of RF plot: 
# a. Assign representing medoid to each sequence (medoid_id)
testgp <- testgp %>% 
  mutate(medoid_id = srfchi[["rf"]][["kmedoid.index"]])

# b. Sort medoid_id according to order of rfplot
testgp <- testgp %>% 
  mutate(medoid_id = factor(medoid_id, levels = srfchi[["rf"]][["medoids"]]))


# c. seqIplot sorted by medoid_id
seqIplot(testseq, border = NA,
         ltext = c(gpstates), 
         missing.color = "#f7f7f7", with.legend = FALSE,
         sortv = testgp$medoid_id)


#### RFplot by clusters (microsim) ####
w <- 7
h <- 7

# generate one RF plot per cluster
# store gp dataframes per cluster
c1 <- gpm %>% 
  filter(chi == l1)

c2 <- gpm %>% 
  filter(chi == l2)

c3 <- gpm %>% 
  filter(chi == l3)

c4 <- gpm %>% 
  filter(chi == l4)

c5 <- gpm %>% 
  filter(chi == l5)

c6 <- gpm %>% 
  filter(chi == l6)

# Cluster 1
seq1 <- seqdef(c1, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi1 <- seqdist(seq1, method = "CHI2", step = max(seqlength(seq1)))

# Select medoids based on distance
srfchi1 <- seqrf(seq1,
                 diss = chi1,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c1.pdf"),
    width=w, height=h)
plot(srfchi1, which.plot = "both", main = l1)
dev.off()

# Cluster 2
seq2 <- seqdef(c2, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi2 <- seqdist(seq2, method = "CHI2", step = max(seqlength(seq2)))

# Select medoids based on distance
srfchi2 <- seqrf(seq2,
                 diss = chi2,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c2.pdf"),
    width=w, height=h)
plot(srfchi2, which.plot = "both", main = l2)
dev.off()


# Cluster 3
seq3 <- seqdef(c3, 2:paste0(max_age+2), # for max_age 300 to column 306, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi3 <- seqdist(seq3, method = "CHI2", step = max(seqlength(seq3)))

# Select medoids based on distance
srfchi3 <- seqrf(seq3,
                 diss = chi3,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c3.pdf"),
    width=w, height=h)
plot(srfchi3, which.plot = "both", main = l3)
dev.off()


# Cluster 4
seq4 <- seqdef(c4, 2:paste0(max_age+2), # for max_age 400 to column 406, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi4 <- seqdist(seq4, method = "CHI2", step = max(seqlength(seq4)))

# Select medoids based on distance
srfchi4 <- seqrf(seq4,
                 diss = chi4,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c4.pdf"),
    width=w, height=h)
plot(srfchi4, which.plot = "both", main = l4)
dev.off()

# Cluster 5
seq5 <- seqdef(c5, 2:paste0(max_age+2), # for max_age 500 to column 506, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi5 <- seqdist(seq5, method = "CHI2", step = max(seqlength(seq5)))

# Select medoids based on distance
srfchi5 <- seqrf(seq5,
                 diss = chi5,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c5.pdf"),
    width=w, height=h)
plot(srfchi5, which.plot = "both", main = l5)
dev.off()


# Cluster 6
seq6 <- seqdef(c6, 2:paste0(max_age+2), # for max_age 100 to column 101, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi6 <- seqdist(seq6, method = "CHI2", step = max(seqlength(seq6)))

# Select medoids based on distance
srfchi6 <- seqrf(seq6,
                 diss = chi6,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c6.pdf"),
    width=w, height=h)
plot(srfchi6, which.plot = "both", main = l6)
dev.off()


# Combine all per-cluster rfplots into one graph
w <- 750
h <- 600

# png(file = paste0(graph.folder, "seqrf_cluster6.png"),
#     width=w, height=h)
pdf(paste0(graph.folder, "seqrf_cluster6.pdf"), width = 8, height = 6)  # Open PDF device
original_par <- par(no.readonly = TRUE) # store original current parameter
par(mfrow = c(3, 2), # 3 rows, 2 columns
    mar = c(3.5, 2, 3 , 2), # margins of each plot
    mgp = c(2, 1, 0)) # margins around axis title, axis labels, and axis line
plot(srfchi1, which.plot = "medoids", skipar = TRUE, main = l1, cex.main = 1, info = "none")
plot(srfchi2, which.plot = "medoids", skipar = TRUE, main = l2, cex.main = 1, info = "none")
plot(srfchi3, which.plot = "medoids", skipar = TRUE, main = l3, cex.main = 1, info = "none")
plot(srfchi4, which.plot = "medoids", skipar = TRUE, main = l4, cex.main = 1, info = "none")
plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1, info = "none", xlab = "Age")
plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1, info = "none", xlab = "Age")
dev.off()
par(mfrow = c(1, 1)) # reset layout


pdf(paste0(graph.folder, "seqrf_both_cluster6.pdf"), 
    width = 8, height = 9)
par(mfrow = c(3, 4), # 3 rows, 4 columns
    mar = c(3.5, 2, 3 , 2), # margins of each plot
    mgp = c(2, 1, 0)) # margins around axis title, axis labels, and axis line
plot(srfchi1, which.plot = "medoids", skipar = TRUE, main = l1, cex.main = 1.1, info = "none")
plot(srfchi1, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi2, which.plot = "medoids", skipar = TRUE, main = l2, cex.main = 1.1, info = "none")
plot(srfchi2, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi3, which.plot = "medoids", skipar = TRUE, main = l3, cex.main = 1.1, info = "none")
plot(srfchi3, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi4, which.plot = "medoids", skipar = TRUE, main = l4, cex.main = 1.1, info = "none")
plot(srfchi4, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi5, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi6, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)
dev.off()
par(original_par) # reset layout



# 2) OM distance with transition rate based costs
omt <- seqdist(seq, method = "OM", indel = 1, sm = "TRATE")

# Select medoids based on distance
srfomt <- seqrf(seq,
                diss = omt,
                sortv = "mds",
                grp.meth = "prop",
                weights = ac$aggWeights)

# RF plot: 
# Plot all k = 100 medoids + average distance of repr. sequences to medoid
png(file = paste0(graph.folder, "seqrf_omt_e.png"),
    width=w, height=h)
plot(srfomt, which.plot = "both")
dev.off()
summary(srfomt)

# For sequence index plot in order of RF plot: 
# a. Assign representing medoid to each sequence (medoid_id)
testgp <- testgp %>% 
  mutate(medoid_id = srfomt[["rf"]][["kmedoid.index"]])

# b. Sort medoid_id according to order of rfplot
testgp <- testgp %>% 
  mutate(medoid_id = factor(medoid_id, levels = srfomt[["rf"]][["medoids"]]))

# c. seqIplot sorted by medoid_id
seqIplot(testseq, border = NA,
         ltext = c(gpstates), 
         missing.color = "#f7f7f7", with.legend = FALSE,
         sortv = testgp$medoid_id)


# RFplot by clusters
seqrfplot(testseq, group = group.p(testgp$chi), 
          diss = omt,
          sortv = "mds", 
          ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
          ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE)

# seqIplot by clusters sorted by medoid_id
seqIplot(testseq, group = group.p(testgp$chi),
         ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
         ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE,
         sortv = testgp$medoid_id)




##### Complexity by cluster (microsimulated) ####

indic1 <- seqindic(seq1, indic=c("cplx"), with.missing=F) # complexity index
indic2 <- seqindic(seq2, indic=c("cplx"), with.missing=F) 
indic3 <- seqindic(seq3, indic=c("cplx"), with.missing=F) 
indic4 <- seqindic(seq4, indic=c("cplx"), with.missing=F) 
indic5 <- seqindic(seq5, indic=c("cplx"), with.missing=F) 
indic6 <- seqindic(seq6, indic=c("cplx"), with.missing=F) 


# store means across full sample and rowbind into one dataframe
indic_mean_mc <- indic1 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  rbind(indic2 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic3 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic4 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic5 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)),
        indic6 %>% summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)))

# Swap rows and columns of indic_mean
tab_ind_clusters_m <- as.data.frame(t(indic_mean_mc))


# Store as docx
knitr::kable(tab_ind_clusters_m, "simple")




#### DESCRIPTION OF CLUSTERS ####

# MEANS
agg <- gpe %>% 
  rbind(gpm) %>% 
  select(source, chi, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(chi, source) %>%
  summarise_all(mean, na.rm = TRUE) #%>% 
# select(cohort, dage:gcage) # exchange 'cohort' for 'dob_year'

# Means
agg_mean <- gpe %>% 
  rbind(gpm) %>% 
  select(source, chi, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(chi, source) %>% 
  summarise_all(mean, na.rm = TRUE)

# Medians (for cont vars)
agg_p50 <- gpe %>% 
  rbind(gpm) %>% 
  select(source, chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(chi, source) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -c(source, chi))

# SD (for cont vars)
agg_sd <- gpe %>% 
  rbind(gpm) %>% 
  select(source, chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(chi, source) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -c(source, chi)) 


# Combine means, medians and SD into one df
agg <- agg_mean %>% 
  left_join(agg_sd, by = c("source", "chi")) %>% 
  left_join(agg_p50, by = c("source", "chi")) %>% 
  select(
    source,
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


# # Swap rows and columns 
# tab1_agg <- agg %>%
#   mutate(across(-c(chi, source), ~ round(.x, 2))) %>%
#   pivot_longer(-c(chi, source), names_to = "variable", values_to = "value") %>%
#   pivot_wider(names_from = source, values_from = value) %>%
#   column_to_rownames("variable") # Moves "variable" column to row names



# Rounds values of indicators
agg_round <- agg %>%
  ungroup() %>% 
  mutate(across(-c(chi, source), ~ round(.x, 2))) 

# Swaps rows and columns
tab2_agg <- agg_round %>%
  arrange(desc(source), chi) %>% # arrange alphabetically
  pivot_longer(cols = -c(chi, source), names_to = "variable", values_to = "value") %>%
  unite("group", chi, source, sep = "_") %>%  # Combine chi and source into a single column
  pivot_wider(names_from = group, values_from = value) %>% 
  column_to_rownames("variable") %>% 
  select(
    starts_with("Cluster 1"),
    starts_with("Cluster 2"),
    starts_with("Cluster 3"),
    starts_with("Cluster 4"),
    starts_with("Cluster 5"),
    starts_with("Cluster 6")
  )



# Number of observations
tab2_n <- gpe %>%
  rbind(gpm) %>% 
  count(chi, source) %>% 
  group_by(source) %>%
  mutate(prop = n / sum(n)*100) %>% 
  pivot_longer(cols = -c(chi, source), names_to = "variable", values_to = "value") %>% 
  mutate(across(c(value), ~ round(.x, ))) %>% 
  unite("group", chi, source, sep = "_") %>%  # Combine chi and source into a single column
  pivot_wider(names_from = group, values_from = value) %>%
  column_to_rownames("variable") %>% 
  as.data.frame() %>% 
  select(
    starts_with("Cluster 1"),
    starts_with("Cluster 2"),
    starts_with("Cluster 3"),
    starts_with("Cluster 4"),
    starts_with("Cluster 5"),
    starts_with("Cluster 6")
  )



##### CREATE TAB 2 ####
tab2 <- tab2_agg %>% 
  rbind(tab2_n) 

tab2 <- rownames_to_column(tab2, var = "Measure") 

write.csv(tab2, paste0(help.folder, "tab2.csv"))


### last line ###