# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# COMPARISON OF SEQUENCES ACROSS COHORTS 
# 1960 & 2000 BIRTH COHORT - AGE RANGE 0-100

# huenteler@demogr.mpg.de
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

agg <- gp60 %>% 
  rbind(gp00) %>% 
  select(dob_year, dage, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  group_by(dob_year) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  mutate(cohort = dob_year) %>% 
  select(cohort, dage:gcage) # exchange 'cohort' for 'dob_year'


# Swap rows and columns of indic_mean
tab1_agg <- agg %>%
  mutate(across(-cohort, ~ round(.x, 2))) %>%
  pivot_longer(-cohort, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  column_to_rownames("variable") # Moves "variable" column to row names

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
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(cohort = 1960) %>% 
  select(-Visited) %>% 
  select(cohort, Lgth:Cplx)

indic_mean00 <- indic00 %>%
  summarise(round(across(everything(), \(x) mean(x, na.rm = TRUE)),2)) %>% 
  mutate(cohort = 2000) %>% 
  select(-Visited) %>% 
  select(cohort, Lgth:Cplx)

indic_mean <- indic_mean60 %>% 
  rbind(indic_mean00) 

# Swap rows and columns of indic_mean
tab1_ind <- indic_mean %>%
  pivot_longer(-cohort, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cohort, values_from = value) %>%
  column_to_rownames("variable")



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

agg <- gp60 %>% 
  rbind(gp00) %>% 
  select(dob_year, chi, dage, pdage, isparent, numkids, cage, isgparent, numgkids, gcage) %>% 
  mutate(cohort = dob_year) %>% 
  group_by(cohort, chi) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(cohort, chi, dage:gcage) # exchange 'cohort' for 'dob_year'


# Swap rows and columns 
tab3_agg <- agg %>%
  ungroup() %>% 
  mutate(across(-c(cohort, chi), ~ round(.x, 2))) %>%
  pivot_longer(-c(cohort, chi), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = c(cohort, chi), values_from = value) %>%
  column_to_rownames("variable") # Moves "variable" column to row names

# Number of observations
tab3_n <- gp60 %>%
  rbind(gp00) %>% 
  mutate(cohort = dob_year) %>% 
  count(cohort, chi) %>% 
  pivot_wider(names_from = c(cohort,chi), values_from = n) 

rownames(tab3_n) <- "N"

tab3 <- tab3_agg %>% 
  rbind(tab3_n)

tab3 <- rownames_to_column(tab3, var = "Variable")


set_flextable_defaults(
  font.size = 11,
  border.color = 'black',
  line_spacing = 1.3,
)

ft3 <- flextable(tab3) %>% 
  add_header_row(colwidths = c(1,6,5), values = c(" ", "1960 birth cohort", "2000 birth cohort")) %>% 
  align(align = "center", part = "header") %>% 
  align(j = 1, align = "left", part = "body") # first column left-align
ft3

save_as_docx("Table 3" = ft3, path = paste0(folder.baseseed, "Tab3.docx"), align = "left")





### last line ###
