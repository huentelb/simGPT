# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SEQUENCE & CLUSTER ANALYSIS 
# 1953 BIRTH COHORT - AGE RANGE 0-100

# huenteler@demogr.mpg.de
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse) 
library(flextable)

# 1. Upper level folder based on simulation base_seed
folder.baseseed <- paste0(folder,"/sim_results_", supfile, "_",base_seed,"_/")
if (!dir.exists(folder.baseseed)) {
  # If not, create the new folder
  dir.create(folder.baseseed)
  cat("Folder created:", folder.baseseed, "\n")
} else {
  cat("Folder already exists:", folder.baseseed, "\n")
}


# load the gp-data for from each simulation from each birth cohorts
for (c in c(1960, 2000)) {
  
  for (i in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) {
    
    load(paste0(folder, "/sim_results_", supfile, "_",base_seed,i,"_/gp", c, max_age, ".RData"))
    assign(paste0("gp_", i), gp)  
    }
    
  # merge the simulations into one gp-dataframe for each birth cohort
  gp <- rbind(gp_1, gp_2, gp_3, gp_4, gp_5, gp_6, gp_7, gp_8, gp_9, gp_10)
    
    
  # store combined gp dataframe in base_seed folder
  save(gp, file = paste0(folder.baseseed, "gp", c, max_age, ".RData"))
  
}



cohort <- 1960
max_age <- 100
ages <- as.character(c(0:max_age))

# Generate folders to store results

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
graph.folder <- paste0(folder.baseseed, cohort, "_", max_age, "/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}



#### OPEN MERGED SIMULATION FILE FOR 2000 COHORT ####

load(paste0(folder.baseseed, "gp", cohort, max_age, ".RData"))




#### SEQUENCE ANALYSIS ####
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
ac <- wcAggregateCases(gp[,6:paste0(max_age+6)])
ac

seq <- seqdef(gp[ac$aggIndex, 6:paste0(max_age+6)], # for max_age 100 to column 106, for max_age 66 to column 72
              weights = ac$aggWeights,
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


png(file = paste0(graph.folder, "seqr100_full.png"),
    width=964, height=556)
seqfplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - ",cohort+max_age,", ",cohort," birth cohort), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", idxs = 1:250)
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

chi_ward10 <- as.clustrange(chi_ward, diss = chi, ncluster = 10, weigths = ac$aggWeights)
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
                         kvals = 2:8,
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

tab_q <- round(summary(chi_pam10, max.rank = 3), 2)
tab_q


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
seqdplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

#### F plots ####
png(file = paste0(graph.folder, "seqF50_4.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster4,
         border = NA, ltext = c(gpstates),  with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 4 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_5.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_6.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_7.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_8.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster8,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 8 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()




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
mc <- chi_pam10$clustering$cluster6[ac$disaggIndex]
med <- as.data.frame(sort(table(mc), decreasing = TRUE))
med1 <- as.character(med[1,1])
med2 <- as.character(med[2,1])
med3 <- as.character(med[3,1])
med4 <- as.character(med[4,1])
med5 <- as.character(med[5,1])
med6 <- as.character(med[6,1])
# med7 <- as.character(med[7,1])
# med8 <- as.character(med[8,1])

# create factor containing medioids incl labels
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5, med6),
                    as.character(c("1","2","3","4","5","6")))
                    # labels = c("Cluster 1 -\n 3-gen family",
                    #            "Cluster 2 -\n Childless",
                    #            "Cluster 3 -\n 4-gen family",
                    #            "Cluster 4 -\n 3-gen (via 2-gen) family",
                    #            "Cluster 5 -\n 2-gen family",
                    #            "Cluster 6",
                    #            "Cluster 7",
                    #            "Cluster 8"))

# attach to dataframe to use as weights in plots
gp$chi <- mc.factor


# save new dataframe for later comparison
save(gp, file = paste0(folder.baseseed, "gp", cohort, max_age, "_chi.RData"))


# generate new sequence object 
seq <- seqdef(gp, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
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
seqdplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"))
dev.off()

png(file = paste0(graph.folder, "seqI_6_lab.png"),
    width=w, height=h)
seqIplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqF50_6_lab.png"),
    width=w, height=h)
seqfplot(seq, group = group.p(gp$chi), border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "mean_plot_6_lab.png"),
    width=w, height=h)
seqmtplot(seq, group = group.p(gp$chi), border = NA,
          ltext = c(gpstates), 
          main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                        "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "seqr_6_lab.png"),
    width=w, height=h)
seqrplot(seq, group = group.p(gp$chi), border = NA,
         ltext = c(gpstates), 
         main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                       "\n ", cohort," birth cohort; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", with.legend = FALSE, diss = chi)
dev.off()


### Relative frequency plot ####

# Check different parameters against theoretically most reasonable: CHI2
# Sorting on the first MDS factor extracted from a dissimilarity matrix built using the CHI-square

# Sample random n=5,000 cases & define sequence object (analysis does not work with full dataset)
set.seed(6789)
testgp <- sample_n(gp, 5000)
testseq <- seqdef(testgp, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
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
plot(srfchi, which.plot = "both")
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


# RFplot by clusters
seqrfplot(testseq, group = group.p(testgp$chi), 
          diss = testchi,
          sortv = "mds", 
          ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
          ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE)

# seqIplot by clusters sorted by medoid_id
seqIplot(testseq, group = group.p(testgp$chi),
         ltext = gpstates, use.layout = TRUE, cex.legend = 1.2,
         ylab = NA, yaxis = FALSE, border = NA, with.legend = FALSE,
         sortv = testgp$medoid_id)



# 2) OM distance with transition rate based costs
omt <- seqdist(testseq, method = "OM", indel = 1, sm = "TRATE")

# Select medoids based on distance
srfomt <- seqrf(testseq,
               diss = omt,
               sortv = "mds",
               grp.meth = "first")

# RF plot: 
# Plot all k = 100 medoids + average distance of repr. sequences to medoid
plot(srfomt, which.plot = "both")


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







#### DESCRIPTION OF CLUSTERS (preliminary) ####

indic <- seqindic(seq, indic=c("lgth", "visited", "trans", "entr", "turb2n", "cplx"), with.missing=F)

indic$cluster <- gp$cluster

indic_mean_cl <- indic %>%
  group_by(cluster) %>% 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

### last line ###

