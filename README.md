# Projecting Generational Placement Trajectories: Empirical and Simulated Populations in Norway
Authored by Bettina Hünteler & Diego Alburez-Gutierrez. For questions contact huenteler@demogr.mpg.de. 

## Objective
We aim to answer two research questions:
1. Can microsimulation be used to validly estimate generational placement trajectories (GPT)?
2. Which typical patterns emerge when considering the full life course (ages 0–100) for the 1960 and 2000 birth cohorts (RQ2)?
3. How do these patterns compare between the cohorts regarding the timing of transitions as well as occurrence and duration of states (RQ3)?

This is the main output we want to produce to answer these questions:
1. Synthetic population register of individuals living in Norway in which we can link parents to their children based on `rsocsim`.
2. Generational placement patterns (see [Hünteler, 2022](https://www.sciencedirect.com/science/article/pii/S104026082100054X)) for two cohorts across two age ranges
    - 1960 birth cohort, age range 0 – 59: for benchmarking against existing historical register data (RQ1),
    - 1960 and 2000 birth cohorts, age range 0 – 100: for projecting generational placement trajectories into the future (RQ2) and comparing them (RQ3). 



## This is how we do it

To reproduce, have the following packages installed:
1. `rsocsim` (find guide to installing and the general idea for and logic of `rsocsim` in [Diego's workshop repository](https://github.com/alburezg/rsocsim_workshop_paa))
2. `TraMineR`, `TraMineRextras`, and `WeightedCluster`


### 1. Input data for microsimulation
*No need for you to repeat this because we stored all the necessary data in the [Input](Input) folder*

1. Fertility data come from the following sources
    - Human Fertility Collection (HFC) for the past periods 1846 – 1966
    - Human Fertility Database (HFD) for the past periods 1967 – 2022
    - United Nations World Population Prospects (UNWPP24) for the past period 2023
    - United Nations World Population Prospects (UNWPP24) for the future periods 2024 – 2100 (medium scenario)

2. Mortality data come from the following sources
    - Human Mortality Database (HMD) for the past periods 1846 – 2023
    - UNWPP24 for the future periods 2024 – 2100

The HFC/HFD data contain the period (1x1) age-specific female fertility rates between ages 14 and 50 (open-ended categories) (already combined and stored as [NORasfrRR.txt](Input/NORasfrRR.txt)).

The HMD data contain the period (1x1) life tables by sex ([fltper_1x1.txt](Input/fltper_1x1.txt) for women and [mltper_1x1.txt](Input/mltper_1x1.txt) for men). 

The UNWPP data come in a different format than the HFC/HFD and HMD data. You can find the code to adjust the UNWPP data to the format of the HFC/HFD and HMD data as well as to append it to the time series of these past time series in the [unwpp_convert.R](unwpp_convert.R) file. You do not need to run it because these data are already stored in the repository, too. 

*The simulations of the main analysis are based on the following data (1846-2100 time series, 1x1), stored in [Input](Input)*
* Fertility (medium scenario): [NORasfrRR_med.txt](Input/NORasfrRR_med.txt)
* Female mortality: [fltper_1x1_med.txt](Input/fltper_1x1_med.txt)
* Male mortality: [mltper_1x1_med.txt](Input/mltper_1x1_med.txt)

* In the appendix, we present results based on simulations with the *high fertility scenario* and the *low fertility scenario*. 


### 2. Set up before you start
Before you start with the microsimulation, please adjust the [setup](setup.R). We have automated many of the storing or labelling processes in the code. Therefore, you should check whether the information provided in the setup matches the simulation-specification etc. that you will be using. You should pay attention to the following options you can set:

*Note: If you leave the code as is provided here, you only need to adjust (5) the base seed and (6) your working directory!*
1. Working directory
2. Names of the input files (no need to change for the main analysis)
3. Name of the [supfile](socsim_NOR.sup) (contains settings for microsimulation)
4. Specifications you will use in the supfile (fertility heterogeneity, birth interval, alpha and beta) so that output will contain correct info. (*Important: changing things here will not change anything about the settings of the microsimulation; the information here is only for labelling correctly! If you want to adjust the microsimulation settings, go to the [supfile](socsim_NOR.sup))*
5. Base seed (within the code, we will run *i* different simulations that will then automatically be starting with setting the *base_seed* + *i*. This will create one output folder for each *i* simulation)
    * For the main analysis we used the base seed *250129* and ten simulation rounds (*i* = 1 to 10). 
6. Folder name (usually your working directory)


### 3. Microsimulation (to build synthetic population register)
*The code to run the following two steps (3. and 4.) is included in [02_simulation.R](02_simulation.R).*

1. First, you convert input data into monthly rates using [functions.R](functions.R) which is started from within [simgpt_DEL.R](simgpt_DEL.R).
2. Second, you run your microsimulation *i* times (set up in a loop). For the main analysis, we ran 10 simulation rounds and merged the resulting synthetic populations to achieve substantial enough cohort sizes. 
    * You can adjust the microsimulation settings in the [supfile](socsim_NOR.sup). For our current output we go with the following setting:
        - Birth interval: 9
        - No fertility heterogeneity
          - Therefore, no option to set alpha and beta
        - Run rates from 1846 for 100 years (1200 months) before starting to move forward in time to set up a population to start from
        - Initial population size = 10,000 (not in supfile but in the main code) of a random age and sex distribution
3. Compare input to simulation rates on the aggregate level. To this end, the code produces two graphs for each simulation (stored in your simulation round's subfolder 'graphs') that compare the age-specific mortality (both sexes) and fertility (only female) rates across different periods. 


### 4. Build dataframes to analyse generational placement trajectories (GPT)
Based on each simulation, you create a dataframe containing the GPT for different subsamples (initiated by a loop varying the birth cohort and age ranges considered):
- *gp196059.RData*: 1960 cohort, age range 0 to 59: for benchmarking against empirical register data
- *gp1960100.RData*: 1960 cohort, age range 0 to 100: for projection into the future
- *gp2000100.RData*: 2000 cohort, age range 0 to 100: for projection into the future

These dataframes are stored in your simulation-subfolders based on *base seed* and simulation round *i*. 

GPT indicate at each age, whether focal is alive and whether they have at least one parent, at least one child, and/or at least one grandchild that is alive. 

The construction of GPT follows the following steps and are applied to each *i* simulation round using a loop:
1. Identify parents and focal's age when the second parent died
2. Identify children and focal's age when the first child was born
3. Identify grandchildren and focal's age when the first grandchild was born
4. Combine this information to define GPT for each age at focal
    - One row for each focal
    - One column for each age of focal containing their generational placement
      - Focal's generational placement can range from being dead to having all respective kin alive at the same time
5. Add GPT to basic demographic information of focal and their parents, children and grandchildren and store as dataframe `gp_i` for each *i* simulation round. 


### 5. Sequence and cluster analysis
The next steps of the analyses all rely on sequence and cluster analysis. We use this, to analyse the GPT and identify a reasonable number of clusters that group together similar GPT and represent typical generational placement patterns. We run various sets of the following analytical steps on our different subsamples from step 4. 

This is the logic:
1. For each subsample, merge the different `gp_i` dataframes for each *i* simulation round into one `gp` dataframe. 
2. Define individual sequences, i.e. our generational placement trajectories
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
      - Additionally: representative sequence plot to reduce overplotting
6. Generate descriptive tables to describe the composition of clusters


##### 5.1 Benchmarking (RQ1)
[03_prep_benchmark.R](03_prep_benchmark.R) and [04_benchmark.R](04_benchmark.R) contain the code to compare the GPT from the synthetic population based on `rsocsim` with the empirical Norwegian register data. Note that you need access to the register data and that the empirical data needs to be prepared on the Norwegian server (GPT are defined using the same logic, but in stata; code not shared here). 

This is the main output these two code files produce:
- Aggregate indicators of demographic events based on both data sources (Table 1)
- Mean time spent in each GPT based on both data sources (Table 1)
- Overall sequence indicators based on both data sources (Table 1)
- BIC differences between sequences based on both data sources (Table 1)
- Cluster characteristics based on both data sources (Table 2)

##### 5.2 Future GPT (RQ2) 
[05_future_gpt60.R](05_future_gpt60.R) and [06_future_gpt00.R](06_future_gpt00.R) analyse the GPT for both birth cohorts (1960 and 2000) from ages 0 to 100, thereby 'projecting' GPT into the future and examining typical GPT for both cohorts.

These two code files produce:
- Graphical representation of typical GPT (Figure 1 and Figure 2)


The [07_future_compare.R](07_future_compare.R) compares the GPT for both cohorts, to investigate change of (typical) GPT over historical time. 

This code file produces:
- Aggregate indicators of demographic events for both cohorts (Table 1)
- Mean time spent in each GPT for both cohorts (Table 1)
- Overall sequence indicators for both cohorts (Table 1)
- BIC differences between sequences contrasting both cohorts (Table 1)
- Cluster characteristics for both cohorts (Table 3)

