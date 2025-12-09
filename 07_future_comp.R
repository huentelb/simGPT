# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# COMPARISON OF SEQUENCES ACROSS COHORTS 
# 1960 & 2000 BIRTH COHORT - AGE RANGE 0-100

# Code written by Bettina Hünteler
# huenteler@demogr.mpg.de
# https://github.com/huentelb/simGPT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse) # for rownames_to_column
library(janitor)


max_age <- 100
ages <- as.character(c(0:max_age))
lab_ages <- paste0("age", ages)

# Folder based on simulation base_seed
folder.baseseed <- paste0(folder,"/sim_results_", supfile, "_",base_seed,"_/")
if (!dir.exists(folder.baseseed)) {
  # If not, create the new folder
  dir.create(folder.baseseed)
  cat("Folder created:", folder.baseseed, "\n")
} else {
  cat("Folder already exists:", folder.baseseed, "\n")
}



# Laod gp data of both cohorts and store as separate dataframes with according suffix
load(paste0(folder.baseseed, "gp1960100.RData"))
gp60 <- gp
load(paste0(folder.baseseed, "gp2000100.RData"))
gp00 <- gp


### OVERALL ####

### AGGREGATE-LEVEL COMPARISONS ####

# Means
agg_mean <- gp60 %>% 
  rbind(gp00) %>% 
  dplyr::select(dob_year, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(dob_year) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  mutate(cohort = dob_year) %>% 
  dplyr::select(cohort, dage:gcage) # exchange 'cohort' for 'dob_year'

# Medians
agg_p50 <- gp60 %>% 
  rbind(gp00) %>% 
  dplyr::select(dob_year, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(dob_year) %>% 
  summarise_all(median, na.rm = TRUE) %>% 
  mutate(cohort = dob_year) %>% 
  dplyr::select(cohort, dage:gcage) %>% # exchange 'cohort' for 'dob_year'
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -cohort)

# SD
agg_sd <- gp60 %>% 
  rbind(gp00) %>% 
  dplyr::select(dob_year, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  group_by(dob_year) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  mutate(cohort = dob_year) %>% 
  dplyr::select(cohort, dage:gcage) %>% # exchange 'cohort' for 'dob_year'
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -cohort)

# Combine means, medians and SD into one df
agg <- agg_mean %>%
  left_join(agg_sd, by = "cohort") %>%
  left_join(agg_p50, by = "cohort") %>%
  dplyr::select (
    cohort,
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
  mutate(across(-cohort, ~ round(.x, 2))) %>%
  pivot_longer(-cohort, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  column_to_rownames("variable") %>% # Moves "variable" column to row names
  mutate(diff = NA, # columns for diff and CI upper and lower bound
         lb = NA,
         ub = NA)



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
tab1_n <- gp60 %>%
  rbind(gp00) %>% 
  mutate(cohort = dob_year) %>% 
  count(cohort) %>% 
  pivot_wider(names_from = cohort, values_from = n)

rownames(tab1_n) <- "N"






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
seq60 <- seqdef(gp60, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
              labels = gplabels,  
              cnames = ages, 
              tick.last = TRUE, 
              xtstep = 5, 
              cpal = cblind, 
              alphabet = gpalpha, 
              states = gpstates,
              missing = "D", right = "DEL")

seq00 <- seqdef(gp00, 6:paste0(max_age+6), # for max_age 100 to column 106, for max_age 66 to column 72
                labels = gplabels,  
                cnames = ages, 
                tick.last = TRUE, 
                xtstep = 5, 
                cpal = cblind, 
                alphabet = gpalpha, 
                states = gpstates,
                missing = "D", right = "DEL")



### COMPARISON OF SEQUENCES ####

#### Sequence indicators ####

# "taking nonvisited states into account -> variance as a predictability indicator of the state duration"
# "ignoring nonvisited states -> measuring the variance of the observed spells." (Ritschard 2023: 2045)

# in our case: number of spells = number of visited states; whenever a new spell begins, a new state is visited
# as one cannot go back to a state that they left (parents are not reborn; in our setting, (grand)children do not die)
# -> recurrence meaningless

indic60 <- seqindic(seq60, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                               "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                               "cplx"), with.missing=F) # complexity index

indic00 <- seqindic(seq00, indic=c("lgth", "visited", "visitp", "transp", # length, states visited + prop, no. transitions
                                   "entr", "meand", "dustd", # mean spell duration + SD (accounting for non-visited states: , "meand2", "dustd2")
                                   "cplx"), with.missing=F) # complexity index


# store mean across full sample
indic_mean60 <- indic60 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),3)) %>% 
  mutate(cohort = 1960) %>% 
  dplyr::select(-Visited) %>%
  dplyr::select(cohort, Lgth:Cplx)

indic_mean00 <- indic00 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),3)) %>% 
  mutate(cohort = 2000) %>% 
  dplyr::select(-Visited) %>% 
  dplyr::select(cohort, Lgth:Cplx)

indic_mean <- indic_mean60 %>% 
  rbind(indic_mean00) 

# Swap rows and columns of indic_mean
tab1_ind <- indic_mean %>%
  pivot_longer(-cohort, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  column_to_rownames("variable") %>% 
  mutate(pval = NA)


# Add significance test for difference
indic_t60 <- indic60 %>% 
  mutate(cohort = 1960)
indic_t00 <- indic00 %>% 
  mutate(cohort = 2000)
indic_t <- indic_t60 %>% 
  rbind(indic_t00)


t1 <- t.test(Lgth ~ cohort, data = indic_t)
t2 <- t.test(Visitp ~ cohort, data = indic_t)
t3 <- t.test(Transp ~ cohort, data = indic_t)
t4 <- t.test(Entr ~ cohort, data = indic_t)
t5 <- t.test(MeanD ~ cohort, data = indic_t)
t6 <- t.test(Dustd ~ cohort, data = indic_t)
t7 <- t.test(Cplx ~ cohort, data = indic_t)

tab1_ind[1,3] <- round(t1[["p.value"]], 4)
tab1_ind[2,3] <- round(t2[["p.value"]], 4)
tab1_ind[3,3] <- round(t3[["p.value"]], 4)
tab1_ind[4,3] <- round(t4[["p.value"]], 4)
tab1_ind[5,3] <- round(t5[["p.value"]], 4)
tab1_ind[6,3] <- round(t6[["p.value"]], 4)
tab1_ind[7,3] <- round(t7[["p.value"]], 4)



#### BIC & LRT ####

# Compare occurrence of states
comp_occ <- round(seqCompare(seq60, seqdata2 = seq00, stat = "all", 
                             method = "OMspell", sm = "INDELS", indel = 2, expcost = 0.5),
                  2)

# Compare timing of statesgp$
comp_time <- round(seqCompare(seq60, seqdata2 = seq00, stat = "all", 
                              method = "CHI2", step = 101),
                   2)

# Compare duration in statesgp$
comp_dur <- round(seqCompare(seq60, seqdata2 = seq00, stat = "all", 
                             method = "OMstran", otto = 0.5, sm = "INDELSLOG"),
                  2)

# Combine all BIC and LRT into one df
bic <- as.data.frame(comp_occ, row.names = "Occurence") %>% 
  rbind(as.data.frame(comp_time, row.names = "Timing"),
        as.data.frame(comp_dur, row.names = "Duration")) %>% 
  select(3,1,2) # select and order columns (BIC, LRT, p-value) 

bic <- rownames_to_column(bic, var = "Measure")


# Swap rows and columns of bic?


#### AVERAGE GP MEAN TIME ####
gp_mt60 <- round(seqmeant(seq60),2) #, serr = T for calculating Var, SD, SE
gp_mt00 <- round(seqmeant(seq00),2)

colnames(gp_mt60) <- "1960"
colnames(gp_mt00) <- "2000"

tab1_mt <- gp_mt60 %>% 
  cbind(gp_mt00) %>% 
  as.data.frame() 

#### continue here ####
# Bootstrapping mean time
gp_boot60 <- seq60 %>% 
  mutate(cohort=1960)
gp_boot00 <- seq00 %>% 
  mutate(cohort=2000) 

gp_boot <- gp_boot60 %>% 
  rbind(gp_boot00)

set.seed(123456) # put into loop to get same samples across all indicators
gp_boot_mt <- do(100)*seqmeant(mosaic::resample(seq60))


  gp_boot_mean <- gp_boot_mean %>%
    mutate(diff = X1960-X2000)
  
  # add CI to Table 1
  tab1_agg[i,3] <- round(tab1_agg[i,1]-tab1_agg[i,2],2) #mean(gp_boot_mean$diff)
  tab1_agg[i,4:5] <- round(confint(gp_boot_mean$diff, level = 0.95),2)
  
}




#### CREATE TAB 1 ####
tab1 <- tab1_agg %>% 
  rbind(tab1_n, tab1_mt, tab1_ind)

tab1 <- rownames_to_column(tab1, var = "Measure")

library(gt)
library(flextable)

set_flextable_defaults(
  font.size = 11,
  border.color = 'black',
  line_spacing = 1.3,
)

ft1 <- flextable(tab1) %>% 
  add_header_row(colwidths = c(1,2), values = c(" ", "Cohort")) %>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft1

save_as_docx("Table 1" = ft1, path = paste0(folder.baseseed, "Tab1.docx"), align = "left")

ft_bic <- flextable(bic)%>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft_bic

save_as_docx("BIC" = ft_bic, path = paste0(folder.baseseed, "BIC.docx"), align = "left")





### CLUSTERS ####

# Laod gp data of both cohorts and store as separate dataframes with according suffix
load(paste0(folder.baseseed, "gp1960100_chi.RData"))
gp60 <- gp
load(paste0(folder.baseseed, "gp2000100_chi.RData"))
gp00 <- gp


### AGGREGATE-LEVEL COMPARISONS ####
# Means
agg_mean <- gp60 %>% 
  rbind (gp00) %>%
  dplyr::select (dob_year, chi, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  mutate(cohort = dob_year) %>% 
  group_by(chi, cohort) %>%
  summarise_all (mean, na.rm = TRUE)

# Medians (for cont vars)
agg_p50 <- gp60 %>% 
  rbind (gp00) %>%
  dplyr::select (dob_year, chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  mutate(cohort = dob_year) %>% 
  group_by(chi, cohort) %>%
  summarise_all (median, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_p50"), .cols = -c(cohort, chi))

# SD (for cont vars)
agg_sd <- gp60 %>% 
  rbind (gp00) %>%
  dplyr::select (dob_year, chi, dage, pdage, numkids, cage, numgkids, gcage) %>% 
  mutate(cohort = dob_year) %>% 
  group_by(chi, cohort) %>%
  summarise_all (sd, na.rm = TRUE) %>% 
  rename_with(.fn = ~ paste0(.x, "_sd"), .cols = -c(cohort, chi))

# Combine means, medians and SD into one df
agg <- agg_mean %>%
  left_join(agg_sd, by = c("cohort", "chi")) %>% 
  left_join(agg_p50, by = c("cohort", "chi")) %>% 
  dplyr::select( 
    cohort, 
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

# Round values of indicators
agg_round <- agg %>% 
  ungroup() %>% 
  mutate(across(-c(chi, cohort), ~round(.x, 2)))

# Swap rows and columns
# sorted in descending order (size) within cohorts
tab3_agg <- agg_round %>% 
  arrange(cohort, chi) %>% 
  pivot_longer(cols = -c(chi, cohort), names_to = "variable", values_to = "value") %>% 
  unite("group", chi, cohort, sep = "_") %>% # combine chi and cohort into a single column
  pivot_wider(names_from = group, values_from = value) %>% 
  column_to_rownames("variable") %>% 
  dplyr::select(
    starts_with("Cluster 1"),
    starts_with("Cluster 2"),
    starts_with("Cluster 3"),
    starts_with("Cluster 4"),
    starts_with("Cluster ") # first both non-parent clusters for 1960, then 2000
  )


# Number of observations
tab3_n <- gp60 %>%
  rbind(gp00) %>% 
  mutate(cohort = dob_year) %>% 
  count(cohort, chi) %>% 
  group_by(cohort) %>% 
  mutate(prop = n / sum(n)*100) %>% 
  pivot_longer(cols = -c(chi, cohort), names_to = "variable", values_to = "value") %>% 
  mutate(across(c(value), ~round(.x, ))) %>% 
  unite("group", chi, cohort, sep = "_") %>% 
  pivot_wider(names_from = group, values_from = value) %>% 
  column_to_rownames("variable") %>%
  dplyr::select(
    starts_with("Cluster 1"),
    starts_with("Cluster 2"),
    starts_with("Cluster 3"),
    starts_with("Cluster 4"),
    starts_with("Cluster ") # first both non-parent clusters for 1960, then 2000
  ) %>% 
  as.data.frame()

tab3 <- tab3_agg %>% 
  rbind(tab3_n)

tab3 <- rownames_to_column(tab3, var = "Variable")


set_flextable_defaults(
  font.size = 11,
  border.color = 'black',
  line_spacing = 1.3,
)

ft3 <- flextable(tab3) %>% 
  add_header_row(colwidths = c(1,2,2,2,2,3), values = c(" ", l1, l2, l3, l4, l5)) %>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft3

save_as_docx("Table 3" = ft3, path = paste0(folder.baseseed, "Tab3.docx"), align = "left")


# Mean time in each state per cluster
by(seq60, gp60$chi, seqmeant)
by(seq00, gp00$chi, seqmeant)


### last line ###
