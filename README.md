# Generational Placement Trajectories in Norway: Combining Empirical and Simulated Data [simGPT]
Code (WIP) to produce output for ongoing *simGPT* project (joint work by Bettina Hünteler and Diego Alburez). 

To reproduce, have the following packages installed:
1. `rsocsim` (find guide to installing in [Diego's workshop repository](https://github.com/alburezg/rsocsim_workshop_paa))
2. `TraMineR`, `TraMineRextras`, and `WeightedCluster`

## This is the output we want to produce
1. Synthetic population register of individuals living in Norway in which we can link parents to their children
2. Generational placement patterns (see [Hünteler, 2022](https://www.sciencedirect.com/science/article/pii/S104026082100054X)) for the birth cohorts 1953 and 2000 from ages 0 to 100 

## This is how we do it

### 1. Input data 
*No need for you to repeat this because we stored all the necessary data for you in the [Input](Input) folder*

1. Fertility data come from the following sources
    - Human Fertility Collection (HFC) for the past periods 1846 – 1966
    - Human Fertility Database (HFD) for the past periods 1967 – 2022
    - United Nations World Population Prospects (UNWPP) for the future periods 2023 – 2100

2. Mortality data come from the following sources
    - Human Mortality Database (HMD) for the past periods 1846 – 2022
    - UNWPP for the future periods 2023 – 2100

The HFC/HFD data contain the period (1x1) age-specific female fertility rates between ages 14 and 50 (open-ended categories) (already combined and stored as [NORasfrRR.txt](Input/NORasfrRR.txt)).

The HMD data contain the period (1x1) life tables by sex ([fltper_1x1.txt](Input/fltper_1x1.txt) for women and [mltper_1x1.txt](Input/mltper_1x1.txt) for men). 

The UNWPP data come in a different format than the HFC/HFD and HMD data. You can find the code to adjust the UNWPP data to the format of the HFC/HFD and HMD data as well as to append it to the time series of these past time series in the [unwpp_convert.R](unwpp_convert.R) file. For the simulations, we chose the medium scenarios for fertility and mortality. You do not need to run it because these data are already stored in the repository, too. 

*The simulations are based on the following data (1846-2100 time series, 1x1), stored in [Input](Input)*
* Fertility: [NORasfrRR_med.txt](Input/NORasfrRR_med.txt)
* Female mortality: [fltper_1x1_med.txt](Input/fltper_1x1_med.txt)
* Male mortality: [mltper_1x1_med.txt](Input/mltper_1x1_med.txt)

### 2. Set up before you start
Before you start with the microsimulation, please adjust the [setup](setup.R). We have automated many of the storing or labelling processes in the code. Therefore, you should check whether the information provided in the setup matches the simulation-specification etc. that you will be using. You should pay attention to the following options you can set:

*Note: If you leave the code as is provided here, you only need to adjust (4) the base seed and (5) your working directory!*
1. Working directory
2. Names of the input files (no need to change, usually)
3. Name of the [supfile](socsim_NOR.sup) (contains settings for microsimulation)
4. Specifications you will use in the supfile (fertility heterogeneity, birth interval, alpha and beta) so that output will contain correct info. (*Important: changing things here will not change anything about the settings of the microsimulation; the information here is only for labelling correctly! If you want to adjust the microsimulation settings, go to the [supfile](socsim_NOR.sup))*
5. Base seed (within the code, we will run i different simulations that will then automatically be starting with the base_seed + i. This will create one output folder for each i simulation)
6. Folder name (usually your working directory)

### 3. Microsimulation (to build synthetic population register)
*The code to run all the remaining steps is included in [simgpt_DEL.R](simgpt_DEL.R).*

1. First, you convert input data into monthly rates using [functions.R](functions.R) which is started from within [simgpt_DEL.R](simgpt_DEL.R).
2. You could adjust the microsimulation settings in the [supfile](socsim_NOR.sup). For our current output we go with the following setting:
    - Birth interval: 9
    - No fertility heterogeneity
      - Therefore, no option to set alpha and beta
    - Run rates from 1846 for 100 years (1200 months) before starting to move forward in time
    - Initial population size = 10,000 (not in supfile but in the main code)


### 4. Build dataset to analyse generational placement trajectories
1. Identify parents and focal's age when the second parent died
2. Identify children and focal's age when the first child was born
3. Identify grandchildren and focal's age when the first grandchild was born
4. Combine in one family rooster combining basic demographic information of focal and their parents, children and grandchildren
5. Generate dataset to analyse generational placement trajectories
    - One column for each age of focal containing their generational placement
    - Focal's generational placement can range from being dead to having all respective kin alive at the same time

### 5. Run sequence and cluster analysis
1. Define individual sequences, i.e. our generational placement trajectories
    - with option "DEL" to delete the positions containing missing values = focal is dead to replace the missing values. See Gabadinho et al. (2010) for more details on the options for handling missing values when defining sequence objects.
3. Calculate pairwise distance between each trajectory
    - Chi<sup>2</sup>-distance measure
4. Cluster them into groups containing trajectories that are similar within and different between the clusters
    - Partitioning around medoids cluster algorithm with Ward as starting point
5. Generate output graphs
    - State distribution plots for full population
    - To compare cluster solutions between different numbers of clusters
      - State distribution plots by cluster 
      - Sequence frequency plots (i most frequent sequences) by cluster
      - Mean time spent in each state by cluster
    - Generate all the above plots with fitted labels for the ideal cluster solution
      - Additionally: representative sequence plot to facilitate labelling
6. Generate descriptive tables to describe composition of clusters
