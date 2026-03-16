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
folder.baseseed <- paste0(folder,"/sim_results_",base_seed,"_/")
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
  pivot_wider(names_from = cohort, values_from = value) 



# Add significance test for difference
gp_boot <- gp60 %>% 
  rbind(gp00) %>% 
  dplyr::select(dob_year, dage, dead_p, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  mutate(cohort = dob_year) 

# bootstrapped means
library(tidyverse)
library(purrr)

# build function for bootstrapping mean, sd and median
boot_ci_diff <- function(df, var) {
  
  reps  <- 1000 # number of sampling
  probs <- c(0.5, 0.025, 0.975) # quantiles: mean diff, upper and lower 95% CI
  
  # calculate three different summary statistics
  stats <- list(
    mean   = function(x) mean(x, na.rm = TRUE),
    sd     = function(x) sd(x, na.rm = TRUE),
    p50 = function(x) median(x, na.rm = TRUE)
  )
  
  # one summary statistic per row
  summary <- lapply(names(stats), function(s) {
    
    # store different quantiles
    qs <- replicate(reps, {
      resampled <- df[sample(nrow(df), replace = TRUE), ] # draw random samples with replacement
      
      # split random samples by group (cohort)
      x1 <- resampled[[var]][resampled$dob_year == 1960]
      x2 <- resampled[[var]][resampled$dob_year == 2000]
      
      # calculate difference of the statistics s between groups (x1 and x2) 
      stats[[s]](x1) - stats[[s]](x2)
    }) %>% 
      quantile(probs = probs, na.rm = TRUE) # calculate quantiles as defined above
    
    # how to structure output
    tibble(
      stat = s,
      diff = round(qs[1],2),
      lb  = round(qs[2],2),
      ub = round(qs[3],2)
    )
  })
  
  dplyr::bind_rows(summary)
}

# apply function to following vars
vars <- c("dage", "dead_p", "pdage", "isparent", "numkids", "cage", "isgparent", "numgkids", "gcage")

diffs_boot <- map_dfr( # rowbind output of function
  vars,
  ~ boot_ci_diff(gp_boot, .x) %>% 
    mutate(var = .x, .before = 1)) 

# adjust names from diffs_boot table to match tab1 
diffs_boot_tab1 <- diffs_boot %>% 
  mutate(suffix = ifelse(stat == "mean", "", 
                         ifelse(stat == "sd", "_sd", "_p50")),
         variable = paste0(var, suffix)) %>% 
  dplyr::select(variable, diff, lb, ub)



tab1_agg <- merge(tab1_agg, diffs_boot_tab1, 
                   by = "variable", 
                   sort = FALSE) %>% 
  column_to_rownames("variable") # Moves "variable" column to row names
  


# Number of observations
tab1_n <- gp60 %>%
  rbind(gp00) %>% 
  mutate(cohort = dob_year) %>% 
  count(cohort) %>% 
  pivot_wider(names_from = cohort, values_from = n) %>% 
  mutate(diff = NA, 
         lb = NA, 
         ub = NA) %>% 
  as.data.frame

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
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(cohort = 1960) %>% 
  dplyr::select(-Visited) %>%
  dplyr::select(cohort, Lgth:Cplx)

indic_mean00 <- indic00 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(cohort = 2000) %>% 
  dplyr::select(-Visited) %>% 
  dplyr::select(cohort, Lgth:Cplx)

indic_mean <- indic_mean60 %>% 
  rbind(indic_mean00) 

# Swap rows and columns of indic_mean
tab1_ind <- indic_mean %>%
  pivot_longer(-cohort, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cohort, values_from = value) 


# Bootstrapping for sequence indicators
boot_seqind_diff <- function(seq_1960, seq_2000, var) {
  
  reps  <- 1000 # number of sampling
  probs <- c(0.5, 0.025, 0.975) # quantiles: mean diff, upper and lower 95% CI
  
  boot_diffs <- replicate(reps, {
    
    resampled1 <- seq_1960[sample(nrow(seq_1960), replace = TRUE),]
    resampled2 <- seq_2000[sample(nrow(seq_2000), replace = TRUE),]

    # select one indicator at a time
    x1 <- as.numeric(resampled1[[var]])
    x2 <- as.numeric(resampled2[[var]])
    
    mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)
  })
  
  qs <- quantile(boot_diffs, probs = probs, na.rm = TRUE)
  
    # how to structure output
    tibble(
      diff = round(qs[1],2),
      lb  = round(qs[2],2),
      ub = round(qs[3],2)
    )}

seqind_names <- c("Lgth", "Visitp", "Transp", "Entr", "MeanD", "DustD", "Cplx")

diffs_indic <- map_dfr( # rowbind output of function
  seqind_names,
  ~ boot_seqind_diff(seq_1960 = indic60, 
                     seq_2000 = indic00, .x) %>% 
    mutate(var = .x, .before = 1)) 

tab1_ind <- cbind(tab1_ind, diffs_indic) %>% 
  dplyr::select(-var) %>% 
  column_to_rownames("variable") # Moves "variable" column to row names


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



# build function for bootstrapping mean time
boot_ci_seqmeant <- function(seq_1960, seq_2000) {
  
  reps  <- 1000
  probs <- c(0.5, 0.025, 0.975)
  
  boot_mat <- replicate(reps, {
    
    idx1 <- sample(seq_len(nrow(seq_1960)), replace = TRUE)
    idx2 <- sample(seq_len(nrow(seq_2000)), replace = TRUE)
    
    m1 <- TraMineR::seqmeant(seq_1960[idx1, ])
    m2 <- TraMineR::seqmeant(seq_2000[idx2, ])
    
    as.numeric(m1 - m2)   # Vektor der Länge 6
  })
  
  # Zeilenweise Quantile → Ergebnis: 6 × 3
  qs <- t(apply(boot_mat, 1, quantile, probs = probs, na.rm = TRUE))
  
  tibble(
    state = seq_len(nrow(qs)),
    diff  = round(qs[, 1],2),
    lb   = round(qs[, 2],2),
    ub  = round(qs[, 3],2)
  )
}

diffs_mt <- boot_ci_seqmeant(
  seq_1960 = seq60,
  seq_2000 = seq00
)

tab1_mt <- cbind(tab1_mt, diffs_mt) %>% 
  dplyr::select(-state)
  



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
  add_header_row(colwidths = c(1,2,3), values = c(" ", "Cohort","")) %>% 
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

# Labels from 06_future_gpt00.R -> adjust if changes
l1 <- as.character(paste0("Cluster 1 -\n 3-gen family"))
l2 <- as.character(paste0("Cluster 2 -\n 4-gen family"))
l3 <- as.character(paste0("Cluster 3 -\n 3-gen (via 2-gen) family"))
l4 <- as.character(paste0("Cluster 4 -\n 2-gen family/fuzzy"))
l5 <- as.character(paste0("Cluster 5 -\n Non-parent"))

# Swap rows and columns
# sorted in descending order (size) within cohorts
tab2_agg <- agg_round %>% 
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
tab2_n <- gp60 %>%
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

tab2 <- tab2_agg %>% 
  rbind(tab2_n)

tab2 <- rownames_to_column(tab2, var = "Variable")


set_flextable_defaults(
  font.size = 11,
  border.color = 'black',
  line_spacing = 1.3,
)

ft2 <- flextable(tab2) %>% 
  add_header_row(colwidths = c(1,2,2,2,2,3), values = c(" ", l1, l2, l3, l4, l5)) %>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft2

save_as_docx("Table 2" = ft2, path = paste0(folder.baseseed, "tab2.docx"), align = "left")


# Mean time in each state per cluster
by(seq60, gp60$chi, seqmeant)
by(seq00, gp00$chi, seqmeant)


### last line ###
