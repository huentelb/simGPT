# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# BENCHMARK MICROSIMULATED AGAINST EMPIRICAL REGISTER DATA
# 1953 BIRTH COHORT - AGE RANGE 0-67 

# ONLY RUNS ON SERVER WITH NORWEGIAN REGISTER DATA

# huenteler@demogr.mpg.de
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### SET-UP #####

library(tidyverse) 
library(foreign)

cohort <- 1953
max_age <- 66
ages <- as.character(c(0:max_age))


# You need to manually upload the synthetic gp data to the Norwegian server 
# Load the data from the folder you stored it in

load(paste0(folder.baseseed, "gp195366.RData"))

# store as gp-dataframe with suffix m (microsimulation)
# and re-order and select columns to match dataframe of empirical data
gpm <- gp %>% 
  select(pid, gp_0:gp_66, dage:gcage, numkids, numgkids, dead_p)



# Generate folders to store results on Norwegian server

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




#### EMPIRICAL REGISTER DATA #####

# load 

# store as gp-dataframe with suffix e (empirical)
gpe <- gp

# add new variable "group" to the df for later comparison with simulated data
gpe <- cbind(gpe, group = 2)



##### MERGE EMPIRICAL and MICROSIMULATED DATA ####
gp <- gpe %>% 
  rbind(gpm)


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
seq <- seqdef(gp, 2:paste0(max_age+2), # for max_age 100 to column 106, for max_age 66 to column 72
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates,
              missing = "D", right = "DEL")


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

### DESCRIPTION OF SEQUENCES ####
indic <- seqindic(seq, indic=c("lgth", "visited", "trans", "entr", "turb2n", "cplx"), with.missing=F)

indic_mean <- indic %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))


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
seqdplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

#### F plots ####
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
by(seq, chi_pam10$clustering$cluster5, seqmeant)
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
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"))
dev.off()

png(file = paste0(graph.folder, "seqI_5_lab.png"),
    width=w, height=h)
seqIplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqf20_5_lab.png"),
    width=w, height=h)
seqfplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:20)
dev.off()

png(file = paste0(graph.folder, "mean_plot_5_lab.png"),
    width=w, height=h)
seqmtplot(seq, group = group.p(gp$chi), border = NA,
          ltext = c(gpstates), 
          main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                        "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "seqr_5_lab.png"),
    width=w, height=h)
seqrplot(seq, group = group.p(gp$chi), border = NA,
         ltext = c(gpstates), 
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
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



#### DESCRIPTION OF CLUSTERS ####
indic <- seqindic(seq, indic=c("lgth", "visited", "trans", "entr", "turb2n", "cplx"), with.missing=F)

indic$cluster <- gp$cluster

indic_mean_cl <- indic %>%
  group_by(cluster) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))


#### COMPARISON OF CLUSTERS ####
gpf <- gp %>% 
  filter(fem == 1)
gpm <- gp %>% 
  filter(fem == 0)
seqf <- seqdef(gpf, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates,
              missing = "D", right = "DEL")
seqm <- seqdef(gpm, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

gp <- gp %>% 
  filter(cohort != "1 1944-48")
  
# Compare occurrence of states
comp_occ <- seqCompare(seq,  group = gp$cohort, stat = "all", method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5)

# Compare timing of statesgp$
comp_time <- seqCompare(seq, group = gp$cohort, stat = "all", method = "CHI2", step = max(seqlength(seq)))

# Compare duration in statesgp$
comp_dur <- seqCompare(seq,  group = gp$cohort, stat = "all", method = "OMstran", otto = 0.5, sm = "INDELSLOG")


seq$fem <- gp$fem
seqiplot(seq, border = NA, group = gp$fem,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("GPT by Sex"),
         missing.color = "#f7f7f7", idxs = 1:100)

seqdplot(seq, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Seq2"),
         missing.color = "#f7f7f7")


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
# save(descr, file = paste0(folder, "/sim_results_", supfile, "_", base_seed, "_/sumfile_",base_seed,".RData"))
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
   title = paste0("Summary Statistics (Simulated data, ", het, " heterogeneous fertility, ", bint,  ", opop size = ", base_seed, ")"),
   out = "csv",
   file = paste0(graph.folder, "sumtable_5"))

### last line ###