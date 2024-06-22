### Generate dataset containing ASFR and ASMR for rsocsim ###
### Add years 2023 - 2100 from UNWPP to HMD and HFC/HFD
### Country Norway

# 1. Prepare data from UNWPP (2023 - 2100) to match structure from HMD and HFC/HFD
# 2. Append UNWPP data to already prepared HMD and HFC/HFD (NORasfrRR.txt)

library(dplyr)
library(ggplot2)

#### FERTILITY ####

# Import UNWPP csv file downloaded from https://population.un.org/wpp/  
unwpp_fert <- read.csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2022_Fertility_by_Age1.csv")

# Read HFD ASFR
NORasfrRR <- read.table(file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/NORasfrRR.txt", 
                   as.is=T, header=T, skip=2, stringsAsFactors=F)

# Wrangle data and keep only some years for comparison between both sources
hfd <- 
  NORasfrRR %>% 
  mutate(Age = case_when(Age == "14-" ~ "14",
                         Age == "50+" ~ "50", 
                         TRUE ~ Age),
         Age = as.numeric(Age),
         Source = "HFD") %>% 
  filter(Year > 1949)

# Select years and country and drop irrelevant vars
unwpp_fer_NOR <- 
  unwpp_fert %>% 
  filter(Location == "Norway", 
         Time > 1949,
         Variant == "Medium" | Variant == "Low" | Variant == "Momentum") %>% 
  mutate(ASFR_o = ASFR,
         # convert so it matches ASFR in HFD
         ASFR = ASFR_o / 1000,
         Source = "UNWPP",
         Age = AgeGrp,
         Year = Time) %>% 
  select(Variant, Year, Age, ASFR, Source)

# Compare ASFR between HFD and UNWPP
bind_rows(hfd,  
          unwpp_fer_NOR %>% filter(Variant == "Medium")) %>% 
  filter(Year == 1950 | Year == 1975 | Year == 2000 | Year == 2015) %>% 
  ggplot(aes(Age, ASFR, group = interaction(Year, Source))) +
  geom_line(aes(colour = as.factor(Year), linetype = Source), linewidth = .8) +
  scale_color_discrete(name = "Year")

# Compare UNWPP ASFR between different scenarios
unwpp_fer_NOR %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = interaction(Year, Variant))) +
  geom_line(aes(colour = as.factor(Year), linetype = Variant), linewidth = .8) +
  scale_color_discrete(name = "Period") +
  scale_linetype_discrete()

# Generate txt file that can be appended to existing NORasfrRR.txt (contains ASFR for past (HFD/HFC) already)

## 1. data.frame containing ages 14- and 50+ with ASFR of 0.0000
# Create a vector containing the Years 2023 - 2100
year <- 2023:2100

# Create a vector for Age with "14-" repeated 2100-2022 times and "50+" repeated 2100-2022 times
age <- c(rep(14, 2100-2022), rep(50, 2100-2022))

# Create a vector for ASFR with 0 repeated *2 times (times for each category) as placeholder
asfr <- rep(0, (2100-2022)*2)

# Combine the vectors into a data.frame
limits <- data.frame(Year = year, Age = age, ASFR = asfr)

# 2. Combine UNWPP dataframe with additional ages (medium variant)
unwpp_fer_med <- 
  bind_rows(limits,
            unwpp_fer_NOR %>%
              filter (Year > 2022,
                      Variant == "Medium")) %>% 
  group_by(Year) %>% 
  arrange(Age, .by_group = TRUE) %>% 
  mutate(Age = as.character(Age),
         Age = case_when(Age == "14" ~ "14-",
                         Age == "50" ~ "50+", 
                         TRUE ~ Age)) 
    
# Compare UNWPP ASFR between different years
unwpp_fer_med %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Period") 


#### APPEND TO EXISTING (past) INPUT RATES (NORasfrRR) ####

NORasfRR_med <- 
  bind_rows(NORasfrRR, 
            unwpp_fer_med %>% 
              select(-Source, - Variant))

# Compare ASFR between different years
NORasfRR_med %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1915 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Period") 

# Export table to txt
# tabs to separate, with row names, and no quotation marks around characters
write.table(NORasfRR_med, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/NORasfrRR_med.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)







#### MALE MORTALITY ####

# Import csv file downloaded from https://population.un.org/wpp/  
unwpp_mlt <- read.csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2022_Life_Table_Complete_Medium_Male_2022-2100.csv")

# Read HMD ASMR
mltper <- read.table(file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/mltper_1x1.txt", 
                  as.is=T, header=T, skip=2, stringsAsFactors=F)

unwpp_mlt_NOR <- 
  unwpp_mlt %>% 
  filter(Location == "Norway") %>% 
  mutate(Year = Time,
         Age = as.numeric(AgeGrpStart)) %>% 
  filter(Year > 2022) %>% 
  select(Year, Age, mx, qx, ax, lx, dx, Lx, Tx, ex) 

unwpp_mlt_NOR %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") 


# Add Ages 101-110 and add NA to content rows
# Round numbers as in mltper_1x1?
# year <- rep(2023:2100, each = 10) # rep times ages to cover
# age <- rep(101:110, 2100-2022) # rep times years to cover
# mx <- rep(0, (2100-2022)*10) # rep times ages over years to cover
# qx <- rep(1, (2100-2022)*10)
# ax <- rep(0, (2100-2022)*10)
# lx <- rep(0, (2100-2022)*10)
# dx <- rep(0, (2100-2022)*10)
# Lx <- rep(0, (2100-2022)*10)
# Tx <- rep(0, (2100-2022)*10)
# ex <- rep(0, (2100-2022)*10)
# 
# highage <- data.frame(Year = year, Age = age, mx, qx, ax, lx, dx, Lx, Tx, ex)
# 
# # Combine UNWPP dataframe with additional Ages
# unwpp_mlt_med <- 
#   bind_rows(highage,
#             unwpp_mlt_NOR) %>% 
#               filter (Year > 2022) %>% 
#   group_by(Year) %>% 
#   arrange(Age, .by_group = TRUE) 

# This chunk doesnt work as it messes up the order of Age -- dont know how to solve yet
# mutate(Age = as.character(Age),
#        Age = case_when(Age == "110" ~ "110+",
#                        TRUE ~ Age)) %>% 
  

#### APPEND TO EXISTING INPUT RATES ####

# 1. Delete Ages 101 - 110+ from HMD data
hmd <- 
  mltper %>% 
  mutate(Age = case_when(Age == "110+" ~ "110",
                         TRUE ~ Age),
         Age = as.numeric(Age)) %>% 
  filter(Age <= 100) %>% 
  # leave (observed) mortality RATE mx (central death rate) as is
  # set PROB of dying between last age and next (qx; age- spec. mortality -> probability) to 1 for last age cat
  mutate(qx = ifelse(Age == 100, 1, qx),
         # adjust ax based on formula (72) in HMD documentation (https://mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf)
         ax = ifelse(Age == 100, 1/mx, ax),
         # formula (76) incl (75) -> p = 1-q -> leave as is
         lx = ifelse(Age == 100, lx, lx),
         # for open age cat dx = lx (77)
         dx = ifelse(Age == 100, lx, dx),
         # for open age cat lx * ax
         Lx = ifelse(Age == 100, lx*ax, Lx),
         # for open age cat Tx = Lx
         Tx = ifelse(Age == 100, Lx, Tx),
         # formula (80)
         ex = ifelse(Age == 100, Tx/lx, ex)) 

mltper_1x1_med <- 
  bind_rows(hmd, 
            unwpp_mlt_NOR) 
# mlttest <- 
#   mltper_1x1_med %>% 
#   mutate(Age = as.character(Age),
#          Age = ifelse(Age == "100", "100+", Age))
# 

mltper_1x1_med %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  # Filter rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(breaks = c(1,20,40,60,80,100,110))

# Export table to txt
# tabs to separate, with row names, and no quotation marks around characters
write.table(mltper_1x1_med, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/mltper_1x1_med.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)






#### FEMALE MORTALITY ####

# Import csv file downloaded from https://population.un.org/wpp/  
unwpp_flt <- read.csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2022_Life_Table_Complete_Medium_Female_2022-2100.csv")

# Read HMD ASMR
fltper <- read.table(file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/fltper_1x1.txt", 
                     as.is=T, header=T, skip=2, stringsAsFactors=F)

unwpp_flt_NOR <- 
  unwpp_flt %>% 
  filter(Location == "Norway") %>% 
  mutate(Year = Time,
         Age = as.numeric(AgeGrpStart)) %>% 
  filter(Year > 2022) %>% 
  select(Year, Age, mx, qx, ax, lx, dx, Lx, Tx, ex) 

unwpp_flt_NOR %>% 
  filter(Year == 2022 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") 


# Add Ages 101-110 and add NA to content rows
# Round numbers as in mltper_1x1?
# year <- rep(2023:2100, each = 10) # rep times ages to cover
# age <- rep(101:110, 2100-2022) # rep times years to cover
# mx <- rep(0, (2100-2022)*10) # rep times ages over years to cover
# qx <- rep(1, (2100-2022)*10)
# ax <- rep(0, (2100-2022)*10)
# lx <- rep(0, (2100-2022)*10)
# dx <- rep(0, (2100-2022)*10)
# Lx <- rep(0, (2100-2022)*10)
# Tx <- rep(0, (2100-2022)*10)
# ex <- rep(0, (2100-2022)*10)
# 
# highage <- data.frame(Year = year, Age = age, mx, qx, ax, lx, dx, Lx, Tx, ex)
# 
# # Combine UNWPP dataframe with additional Ages
# unwpp_mlt_med <- 
#   bind_rows(highage,
#             unwpp_mlt_NOR) %>% 
#               filter (Year > 2022) %>% 
#   group_by(Year) %>% 
#   arrange(Age, .by_group = TRUE) 

# This chunk doesnt work as it messes up the order of Age -- dont know how to solve yet
# mutate(Age = as.character(Age),
#        Age = case_when(Age == "110" ~ "110+",
#                        TRUE ~ Age)) %>% 


#### APPEND TO EXISTING INPUT RATES ####

# 1. Delete Ages 101 - 110+ from HMD data
hmd <- 
  fltper %>% 
  mutate(Age = case_when(Age == "110+" ~ "110",
                         TRUE ~ Age),
         Age = as.numeric(Age)) %>% 
  filter(Age <= 100) %>% 
  # leave (observed) mortality RATE mx (central death rate) as is
  # set PROB of dying between last age and next (qx; age- spec. mortality -> probability) to 1 for last age cat
  mutate(qx = ifelse(Age == 100, 1, qx),
         # adjust ax based on formula (72) in HMD documentation (https://mortality.org/File/GetDocument/Public/Docs/MethodsProtocolV6.pdf)
         ax = ifelse(Age == 100, 1/mx, ax),
         # formula (76) incl (75) -> p = 1-q -> leave as is
         lx = ifelse(Age == 100, lx, lx),
         # for open age cat dx = lx (77)
         dx = ifelse(Age == 100, lx, dx),
         # for open age cat lx * ax
         Lx = ifelse(Age == 100, lx*ax, Lx),
         # for open age cat Tx = Lx
         Tx = ifelse(Age == 100, Lx, Tx),
         # formula (80)
         ex = ifelse(Age == 100, Tx/lx, ex)) 

fltper_1x1_med <- 
  bind_rows(hmd, 
            unwpp_flt_NOR)

yrs_plot <- c("[1847,1852)", "[1877,1882)", "[1907,1912)", "[1937,1942)", "[1967,1972)", "[1997,2002)", "[2017,2022)")

fltper_1x1_med %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  # Filter rates of 0, infinite (N_Deaths/0_Pop) and NaN (0_Deaths/0_Pop) values
  filter(mx != 0 & !is.infinite(mx) & !is.nan(mx)) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(breaks = c(1,20,40,60,80,100,110))

# Export table to txt
# tabs to separate, with row names, and no quotation marks around characters
write.table(fltper_1x1_med, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/fltper_1x1_med.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)

