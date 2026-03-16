# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SEQUENCE & CLUSTER ANALYSIS on JOINT SAMPLE
# 1960 and 2000 BIRTH COHORTS - AGE RANGE 0-100

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse) 
library(flextable)


max_age <- 100
ages <- as.character(c(0:max_age))



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
graph.folder <- paste0(folder.baseseed, "joint_", max_age, "/")
if (!dir.exists(graph.folder)) {
  # If not, create the new folder
  dir.create(graph.folder)
  cat("Folder created:", graph.folder, "\n")
} else {
  cat("Folder already exists:", graph.folder, "\n")
}



#### OPEN MERGED SIMULATION FILES ####

# Laod gp data of both cohorts and store as separate dataframes with according suffix
load(paste0(folder.baseseed, "gp1960100.RData"))
gp60 <- gp
load(paste0(folder.baseseed, "gp2000100.RData"))
gp00 <- gp

gp <- gp60 %>% 
  rbind(gp00)

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
         main = paste0("Simulated Data \n (1846 - 2100, 1960 + 2000 birth cohorts), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", with.missing = T)
dev.off()


png(file = paste0(graph.folder, "seqr100_full.png"),
    width=964, height=556)
seqfplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - 2100, 1960 + 2000 birth cohorts), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", idxs = 1:250)
dev.off()

png(file = paste0(graph.folder, "seqi100_full.png"),
    width=964, height=556)
seqiplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
         main = paste0("Simulated Data \n (1846 - 2100, 1960 + 2000 birth cohorts), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop), 
         missing.color = "#f7f7f7", idxs = 1:100)
dev.off()


# Meantime in general pop
seqmeant(seq)
png(file = paste0(graph.folder, "meant_full.png"),
    width=964, height=556)
seqmtplot(seq, border = NA, ltext = c(gpstates), with.legend = FALSE, 
          main = paste0("Simulated Data \n (1846 - 2100, 1960 + 2000 birth cohorts), ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop))
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
saveRDS(chi_ward10, file = paste0(graph.folder, "chiward10.RData"))

chi_ward10 <- readRDS(paste0(graph.folder, "chiward10.RData"))

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
     legendpos="topright",
     main = paste0("Simulated data, ", het, " fertility heterogeneity, ", bint,  ", opop size = ", size_opop, 
                   "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"))
dev.off()

tab_q <- round(summary(chi_pam10, max.rank = 3), 2)
tab_q


#### D plots ####
# State distribution plots by clusters
# CHI2

w = 2000
h = 1250

png(file = paste0(graph.folder, "seqD_5.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_6.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_7.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqD_8.png"),
    width=w, height=h)
seqdplot(seq, group = chi_pam10$clustering$cluster8,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 8 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7")
dev.off()

#### F plots ####
png(file = paste0(graph.folder, "seqF50_5.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster5,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 5 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_6.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster6,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 6 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_7.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster7,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 7 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()

png(file = paste0(graph.folder, "seqF50_8.png"),
    width=w, height=h)
seqfplot(seq, group = chi_pam10$clustering$cluster8,
         border = NA, ltext = c(gpstates),   with.legend = FALSE, cex.axis = 2,
         main = paste0("Chi PAM: 8 Clusters \n sim., ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                       "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
         missing.color = "#f7f7f7", idxs = 1:50)
dev.off()




# Mean time spent in each state by cluster
by(seq, chi_pam10$clustering$cluster8, seqmeant)
png(file = paste0(graph.folder, "mean_plot_8.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster8, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 8 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()

png(file = paste0(graph.folder, "mean_plot_7.png"),
    width=964, height=556)
seqmtplot(seq, group = chi_pam10$clustering$cluster7, border = NA,
          ltext = c(gpstates), main = paste0("Chi Ward: 7 Clusters sim. ", het, " het. fert., ", bint,  ", opop size = ", size_opop, 
                                             "\nboth birth cohorts; alpha = ", alpha, ", beta = ", beta, " (", base_seed, ")"),
          missing.color = "#f7f7f7")
dev.off()




#### LABELLED GRAPH FOR OPTIMAL CLUSTER SOLUTION ####

# We extract X clusters and re-label them from 1 to X to replace the medoid identifiers

# identify medoids sorted according to separate cluster solution of 1960 cohort
mc <- chi_pam10$clustering$cluster7[ac$disaggIndex]
med <- as.data.frame(sort(table(mc), decreasing = TRUE))
med1 <- as.character(med[1,1]) # 3-gen 
med2 <- as.character(med[3,1]) # 4-gen
med3 <- as.character(med[2,1]) # 3-gen via 2-gen
med4 <- as.character(med[5,1]) # 2-gen fuzzy
med5 <- as.character(med[4,1]) # Non-parent
med6 <- as.character(med[7,1]) # Non-parent + early death
med7 <- as.character(med[6,1]) # new cluster

# store size of clusters for each cluster to add to titles
propmed <- as.data.frame(sort(prop.table(table(mc)), decreasing = TRUE))
propmed1 <- round(propmed[1,2], digits = 2)*100 # adjust order of clusters also here!
propmed2 <- round(propmed[3,2], digits = 2)*100
propmed3 <- round(propmed[2,2], digits = 2)*100
propmed4 <- round(propmed[5,2], digits = 2)*100
propmed5 <- round(propmed[4,2], digits = 2)*100
propmed6 <- round(propmed[7,2], digits = 2)*100
propmed7 <- round(propmed[6,2], digits = 2)*100

# create factor containing medoids incl labels
mc.factor <- factor(mc, levels = c(med1, med2, med3, med4, med5, med6, med7),
                    as.character(c("1","2","3","4","5","6","7")))


# store labels as values for later use
l1 <- as.character(paste0("Cluster 1 -\n 3-gen family (", propmed1, "%)"))
l2 <- as.character(paste0("Cluster 2 -\n 4-gen family (", propmed2, "%)"))
l3 <- as.character(paste0("Cluster 3 -\n 3-gen (via 2-gen) family (", propmed3, "%)"))
l4 <- as.character(paste0("Cluster 4 -\n 2-gen family/fuzzy (", propmed4, "%)"))
l5 <- as.character(paste0("Cluster 5 -\n Non-parent (", propmed5, "%)"))
l6 <- as.character(paste0("Cluster 6 -\n Non-parent + early death (", propmed6, "%)"))
l7 <- as.character(paste0("Cluster 7 -\n 3-gen family + early death (", propmed7, "%)"))

# attach to dataframe to use as weights in plots
gp$chi <- factor(mc.factor,
                 labels = c(l1,
                            l2,
                            l3,
                            l4,
                            l5,
                            l6,
                            l7))


# save new dataframe for later comparison
save(gp, file = paste0(folder.baseseed, "gp_joint_chi.RData"))

load(paste0(folder.baseseed, "gp_joint_chi.RData"))

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

w = 2000
h = 1250

# different plots with labels
png(file = paste0(graph.folder, "seqD_7_lab.png"),
    width=w, height=h)
seqdplot(seq, group = gp$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2)
dev.off()

png(file = paste0(graph.folder, "seqI_7_lab.png"),
    width=w, height=h)
seqIplot(seq, group = gp$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7")
dev.off()

png(file = paste0(graph.folder, "seqF100_7_lab.png"),
    width=w, height=h)
seqfplot(seq, group = gp$chi, border = NA,
         ltext = gpstates, with.legend = FALSE, cex.axis = 2,
         missing.color = "#f7f7f7", idxs = 1:100)
dev.off()

by(seq, gp$chi, seqmeant)
png(file = paste0(graph.folder, "mean_plot_7_lab.png"),
    width=w, height=h)
seqmtplot(seq, group = gp$chi, border = NA,
          ltext = c(gpstates), 
          missing.color = "#f7f7f7", with.legend = FALSE)
dev.off()
# 
# png(file = paste0(graph.folder, "seqr_8_lab.png"),
#     width=w, height=h)
# seqrplot(seq, group = gp$chi, border = NA,
#          ltext = c(gpstates), 
#          missing.color = "#f7f7f7", with.legend = FALSE, diss = chi)
# dev.off()


### Relative frequency plot ####

w <- 1000
h <- 625

# Check different parameters against theoretically most reasonable: CHI2
# Sorting on the first MDS factor extracted from a dissimilarity matrix built using the CHI-square

# Sample random n=5,000 cases & define sequence object (analysis does not work with full dataset)
set.seed(2407)
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
png(file = paste0(graph.folder, "seqrf.png"),
    width=w, height=h)
plot(srfchi, which.plot = "both", main = paste0("Random sample from joint birth cohorts (n = 5,000)")) 
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


#### RFplot by clusters ####
w <- 7
h <- 7

# generate one RF plot per cluster
# store gp dataframes per cluster
c1 <- gp %>% 
  filter(chi == l1)

c2 <- gp %>% 
  filter(chi == l2)

c3 <- gp %>% 
  filter(chi == l3)

c4 <- gp %>% 
  filter(chi == l4)

c5 <- gp %>% 
  filter(chi == l5)

c6 <- gp %>% 
  filter(chi == l6)

c7 <- gp %>% 
  filter(chi == l7)


# Cluster 1
seq1 <- seqdef(c1, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
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
seq2 <- seqdef(c2, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
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
seq3 <- seqdef(c3, 6:paste0(max_age+6), # for max_age 300 to column 306, for max_age 66 to column 72
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
seq4 <- seqdef(c4, 6:paste0(max_age+6), # for max_age 400 to column 406, for max_age 66 to column 72
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
seq5 <- seqdef(c5, 6:paste0(max_age+6), # for max_age 500 to column 506, for max_age 66 to column 72
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
seq6 <- seqdef(c6, 6:paste0(max_age+6), # for max_age 100 to column 101, for max_age 66 to column 72
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


# Cluster 7
seq7 <- seqdef(c7, 6:paste0(max_age+6), # for max_age 100 to column 101, for max_age 66 to column 72
               labels = gplabels,  
               cnames = ages, 
               tick.last = TRUE, 
               xtstep = 5, 
               cpal = cblind, 
               alphabet = gpalpha, 
               states = gpstates,
               missing = "D", right = "DEL")

# CHI2 distance
chi7 <- seqdist(seq7, method = "CHI2", step = max(seqlength(seq7)))

# Select medoids based on distance
srfchi7 <- seqrf(seq7,
                 diss = chi7,
                 sortv = "mds",
                 grp.meth = "first")

pdf(file = paste0(graph.folder, "seqrf_c7.pdf"),
    width=w, height=h)
plot(srfchi7, which.plot = "both", main = l7)
dev.off()






# Combine all per-cluster rfplots into one graph
w <- 750
h <- 600

# png(file = paste0(graph.folder, "seqrf_cluster6.png"),
#     width=w, height=h)
pdf(paste0(graph.folder, "seqrf_cluster7.pdf"), 
    width = 8, height = 6)

original_par <- par(no.readonly = TRUE) # store original current parameter

par(mfrow = c(4, 2), # 3 rows, 2 columns
    mar = c(3.5, 2, 3 , 2), # margins of each plot
    mgp = c(2, 1, 0)) # margins around axis title, axis labels, and axis line
plot(srfchi1, which.plot = "medoids", skipar = TRUE, main = l1, cex.main = 1, info = "none")
plot(srfchi2, which.plot = "medoids", skipar = TRUE, main = l2, cex.main = 1, info = "none")
plot(srfchi3, which.plot = "medoids", skipar = TRUE, main = l3, cex.main = 1, info = "none")
plot(srfchi4, which.plot = "medoids", skipar = TRUE, main = l4, cex.main = 1, info = "none")
plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1, info = "none")
plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1, info = "none", xlab = "Age")
plot(srfchi7, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1, info = "none", xlab = "Age")
dev.off()
par(original_par) # reset layout


pdf(paste0(graph.folder, "seqrf_both_cluster7.pdf"), 
    width = 8, height = 9)
par(mfrow = c(4, 4), # 3 rows, 4 columns
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

plot(srfchi5, which.plot = "medoids", skipar = TRUE, main = l5, cex.main = 1.1, info = "none")
plot(srfchi5, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi6, which.plot = "medoids", skipar = TRUE, main = l6, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi6, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

plot(srfchi7, which.plot = "medoids", skipar = TRUE, main = l7, cex.main = 1.1, info = "none", xlab = "Age")
plot(srfchi7, which.plot = "diss.to.med", skipar = TRUE, cex.main = 1)

dev.off()
par(original_par) # reset layout





### Regression by cohort ####
# Package for multinomial logistic regression
library(nnet)
# To transpose output results
library(broom)
library(gtsummary)
library(ggeffects)
library(marginaleffects)

# Set categorical variables as factors
gp_reg <- gp %>%
  mutate(
    cohort = factor(dob_year, 
                       labels = c("1960",
                                  "2000")),
    chi = factor(chi,
                 labels = c("C1 - 3-gen",
                            "C2 - 4-gen",
                            "C3 - 3- via 2-gen",
                            "C4 - 2-gen/fuzzy",
                            "C5 - Non-parent",
                            "C6 - Non-parent, early death",
                            "C7 - 3-gen, early death")),
  ) %>% 
  dplyr::select(chi, cohort)

fit_basic <- multinom(chi ~ cohort, data = gp_reg)
tidy(fit_basic, conf.int=TRUE)
tbl_regression(fit_basic, exp = TRUE)

ggaverage(fit_basic, terms = "cohort", ci_level = 0.95) %>%
  plot() 

ggeffect(fit_basic, terms = "cohort", ci_level = 0.95) %>%
  plot() 

contrast <- avg_comparisons(fit_basic)

pdf(paste0(graph.folder, "contrast.pdf"), 
    width = 8, height = 6)
ggplot(data = contrast,
       aes(x = estimate, 
           y = group)) +
  geom_vline(xintercept = 0, color = "darkred") +
  geom_point() +
  geom_errorbar(aes(xmin = conf.low,
                    xmax = conf.high),
                width = 0.15) +
  labs(
    x = "Contrast: 2000 - 1960 birth cohort",
    y = ""
  )+
  scale_y_discrete(limits = rev(levels(contrast$group))) + # reverse order of displayed clusters
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", size = 11),  
    axis.title = element_text(color = "black", size = 11)   
  )
dev.off()

# Save results
pprob_cohort <- ggeffect(fit_basic, terms = "cohort")

# Build plot
ggplot(data = pprob_cohort,
       aes(x = x, y = predicted,
           color = response.level, group = response.level)) +
  # geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high,
                    color = response.level,
                    group = response.level),
                width = .05) +
  scale_color_brewer(palette = "Dark2",
                     name = "",
                     labels = c(l1, l2, l3, l4, l5, l6, l7)) +
  labs(
    x = "Cohort",
    y = "Probability"
  ) +
  # Set the theme
  theme_minimal() +
  theme(
    legend.position = "none", 
    axis.title = element_text(size = 14) # increase axis title size
  ) +
  facet_wrap("response.level",
             nrow = 4)
  

# Rename the levels of response.level to "C1", "C2", etc.
pprob_cohort$response.level <- factor(pprob_cohort$response.level,
                                      levels = unique(pprob_cohort$response.level),
                                      labels = paste0("Cluster ", 1:length(unique(pprob_cohort$response.level))))

pdf(paste0(graph.folder, "predprob.pdf"), 
    width = 8, height = 6)
ggplot(data = pprob_cohort,
       aes(x = x, y = predicted)) +
  geom_bar(
    stat = "identity",
    fill = "grey",
    color = "black",
    position = position_dodge()
  ) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = .05,
    color = "black",
    position = position_dodge(.9)
  ) +
  geom_text(
    aes(label = round(predicted, 2)),
    position = position_dodge(.9),
    vjust = 2.5,
    size = 3
  ) +
  labs(
    x = "Cohort",
    y = "Probability"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap(
    ~ response.level,
    nrow = 1
  )
dev.off()



#### DESCRIPTION OF CLUSTERS ####

# Means
agg_mean <- gp %>% 
  select(chi, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(chi) %>% 
  summarise_all(mean, na.rm = TRUE)

# Medians (for cont vars)
agg_p50 <- gp %>% 
  select(chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(chi) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -c(chi))

# SD (for cont vars)
agg_sd <- gp %>% 
  select(chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(chi) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -c(chi)) 


# Combine means, medians and SD into one df
agg <- agg_mean %>% 
  left_join(agg_sd, by = c("chi")) %>% 
  left_join(agg_p50, by = c("chi")) %>% 
  select(
    chi, 
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



tab_a9_agg <- agg %>%
  mutate(across(-chi, ~ round(.x, 2))) %>%
  pivot_longer(-chi, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = chi, values_from = value) 


# Number of observations
tab_a9_n <- gp %>%
  group_by(chi) %>% 
  count() %>% 
  pivot_wider(names_from = chi, values_from = n) %>% 
  as.data.frame

rownames(tab_a9_n) <- "N"
tab_a9_n <- rownames_to_column(tab_a9_n, var = "variable")



##### CREATE TAB a9 ####
tab_a9 <- tab_a9_agg %>% 
  rbind(tab_a9_n) 

write.csv(tab_a9, paste0(graph.folder, "tab_a9.csv"))



### last line ###

