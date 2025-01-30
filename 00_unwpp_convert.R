### Generate dataset containing ASFR and ASMR for rsocsim ###
### Add years 2023 - 2100 from UNWPP to HMD and HFC/HFD
### Country Norway

# 1. Prepare data from UNWPP (2023 - 2100) to match structure from HMD and HFC/HFD
# 2. Append UNWPP data to already prepared HMD and HFC/HFD (NORasfrRR.txt)

library(dplyr)
library(ggplot2)
library(readxl)

#### FERTILITY ####

# Import UNWPP csv file downloaded from https://population.un.org/wpp/  
unwpp_fert <- read.csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2024_Fertility_by_Age1.csv.gz")
# to fill year 2023 (not yet contained in HFD, neither in 2022 'low' and 'momentum' variants)
# unwpp_fert22 <- read.csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2022_Fertility_by_Age1.csv")

# Read HFC ASFR
NORasfrRR <- read.table(file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/NORasfrRR.txt", 
                   as.is=T, header=T, skip=2, stringsAsFactors=F)


# Wrangle data and keep only some years for comparison between both sources
hfc <- 
  NORasfrRR %>% 
  mutate(Age = case_when(Age == "14-" ~ "14",
                         Age == "50+" ~ "50", 
                         TRUE ~ Age),
         Age = as.numeric(Age),
         Source = "HFC") %>% 
  filter(Year > 1949)


# Select years and country and drop irrelevant vars
unwpp_fert_NOR <- 
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
bind_rows(hfc,  
          unwpp_fert_NOR %>% filter(Variant == "Medium")) %>% 
  filter(Year == 1950 | Year == 1975 | Year == 2000 | Year == 2015) %>% 
  ggplot(aes(Age, ASFR, 
             group = interaction(Year, Source))) +
  geom_line(aes(colour = as.factor(Year), 
                linetype = Source), linewidth = .8) +
  scale_color_discrete(name = "Year") +
  labs(title = paste0("Comparing ASFR between HFC and UNWPP (", Sys.Date(), ")"))

# Compare UNWPP ASFR between different scenarios
unwpp_fer_NOR %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, 
             group = interaction(Year, Variant))) +
  geom_line(aes(colour = as.factor(Year), 
                linetype = Variant), linewidth = .8) +
  scale_color_discrete(name = "Year") +
  scale_linetype_discrete() +
  labs(title = paste0("Comparing UNWPP ASFR between different scenarios (", Sys.Date(), ")"))


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

# 2. Combine UNWPP dataframe with additional ages 
# 2.1 medium variant from 2023
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

# 2.2 low variant (only from 2024; 2023 <- Medium = historical)
unwpp_fer_low <- 
  bind_rows(limits,
            unwpp_fer_NOR %>%
              filter (Year == 2023,
                      Variant == "Medium"),
            unwpp_fer_NOR %>%
              filter (Year > 2023,
                      Variant == "Low")) %>% 
  group_by(Year) %>% 
  arrange(Age, .by_group = TRUE) %>% 
  mutate(Age = as.character(Age),
         Age = case_when(Age == "14" ~ "14-",
                         Age == "50" ~ "50+", 
                         TRUE ~ Age)) 

# 2.3 momentum variant (only from 2024; 2023 <- Medium = historical)
unwpp_fer_hi <- 
  bind_rows(limits,
            unwpp_fer_NOR %>%
              filter (Year == 2023,
                      Variant == "Medium"),
            unwpp_fer_NOR %>%
              filter (Year > 2023,
                      Variant == "Momentum")) %>% 
  group_by(Year) %>% 
  arrange(Age, .by_group = TRUE) %>% 
  mutate(Age = as.character(Age),
         Age = case_when(Age == "14" ~ "14-",
                         Age == "50" ~ "50+", 
                         TRUE ~ Age)) 
    
# Compare UNWPP ASFR between different years
med <- unwpp_fer_med %>% 
  filter(Year == 2023 | Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, 
             group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year") +
  labs(title = paste0("Medium"))

low <- unwpp_fer_low %>% 
  filter(Year == 2023 | Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year")  +
  labs(title = paste0("Low"))

hi <- unwpp_fer_hi %>% 
  filter(Year == 2023 | Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year")  +
  labs(title = paste0("High (", Sys.Date(), ")"))

library(patchwork)
low + med + hi +
  plot_layout(guides = 'collect')

#### APPEND TO EXISTING (past) INPUT RATES (NORasfrRR) ####

NORasfRR_med <- 
  bind_rows(NORasfrRR, 
            unwpp_fer_med %>% 
              select(-Source, - Variant))

NORasfRR_low <- 
  bind_rows(NORasfrRR, 
            unwpp_fer_low %>% 
              select(-Source, - Variant))

NORasfRR_hi <- 
  bind_rows(NORasfrRR, 
            unwpp_fer_hi %>% 
              select(-Source, - Variant))

# Compare ASFR between different years
med <- NORasfRR_med %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1915 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year")   +
  labs(title = paste0("Medium"))

low <- NORasfRR_low %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1915 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year")   +
  labs(title = paste0("Low"))

hi <- NORasfRR_hi %>% 
  filter(Year == 1850 | Year == 1900 | Year == 1915 | Year == 1950 | Year == 2000 | Year == 2050 | Year == 2100) %>% 
  ggplot(aes(Age, ASFR, group = Year)) +
  geom_line(aes(colour = as.factor(Year)), linewidth = .8) +
  scale_color_discrete(name = "Year")   +
  labs(title = paste0("High (", Sys.Date(), ")"))

low + med + hi +
  plot_layout(guides = 'collect')

# Export table to txt
# tabs to separate, with row names, and no quotation marks around characters
write.table(NORasfRR_med, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/NORasfrRR_med.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)
write.table(NORasfRR_low, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/NORasfrRR_low.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)
write.table(NORasfRR_hi, file = "/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/simGPT/Input/NORasfrRR_hi.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)






#### MALE MORTALITY ####

# Import csv file downloaded from https://population.un.org/wpp/  
unwpp_mlt <- read_csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2024_Life_Table_Complete_Medium_Male_2024-2100.csv.gz")

# Read HMD ASMR
hmd_m <- read_xlsx("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/HMD_mlt_1x1.xlsx",
                   skip = 2)


unwpp_mlt_NOR <- 
  unwpp_mlt %>% 
  filter(Location == "Norway") %>% 
  mutate(Year = Time,
         Age = as.numeric(AgeGrpStart)) %>% 
  filter(Year >= 2024) %>% 
  select(Year, Age, mx, qx, ax, lx, dx, Lx, Tx, ex) 

unwpp_mlt_NOR %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") 




#### APPEND TO EXISTING INPUT RATES ####

# 1. Delete Ages 101 - 110+ from HMD data
hmd_m_mod <- 
  hmd_m %>% 
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
  bind_rows(hmd_mod, 
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
unwpp_flt <- read_csv("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/WPP2024_Life_Table_Complete_Medium_Female_2024-2100.csv.gz")


# Read HMD ASMR
hmd_f <- read_xlsx("/Users/Bettina/sciebo/projects/GenerationalPlacements/SimGPT/analysis/data/HMD_flt_1x1.xlsx",
                   skip = 2)

unwpp_flt_NOR <- 
  unwpp_flt %>% 
  filter(Location == "Norway") %>% 
  mutate(Year = Time,
         Age = as.numeric(AgeGrpStart)) %>% 
  filter(Year >= 2024) %>% 
  select(Year, Age, mx, qx, ax, lx, dx, Lx, Tx, ex) 


unwpp_flt_NOR %>% 
  filter(Year == 2025 | Year == 2050 | Year == 2075 | Year == 2100) %>% 
  ggplot(aes(Age, mx, group = Year)) +
  geom_line(aes(color = as.factor(Year))) +
  scale_color_discrete(name = "Period") +
  scale_y_continuous(trans = "log10") 




#### APPEND TO EXISTING INPUT RATES ####

# 1. Delete Ages 101 - 110+ from HMD data
hmd_f_mod <- 
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
  bind_rows(hmd_f_mod, 
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

