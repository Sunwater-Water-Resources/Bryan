# Monte Carlo config file
This is the main config file for the Monte Carlo simulations and uses the keys in the tables below.  

**Table 1: Main keys in the Monte Carlo config file**
| Key | Description |
| ----------- | ----------- |
| ```scheme_config```| A dictionary parameterising the Monte Carlo scheme - see Table 2. The information in the table is used to set up the sampling framework for rainfall depths and storm temporal patterns.|
|```tpt_quantile_analysis```| A dictionary containing the information used to analyse the results of the Monte Carlo analysis using the Total Probability Theorem - see Table 3.  
|  ```confidence_intervals``` | A dictionary containing the information used in the analysis of the Monte Carlo results to determine the  confidence intervals - ***this is a work in progress***.  |

**Table 2: Keys used in the ```scheme_config``` key**
| Key | Description |
| ----------- | ----------- |
| ```lower_aep``` | The lower bound AEP for the range of AEPs to be sampled from the rainfall depths. |
| ```upper_aep``` | The upper bound AEP for the range of AEPs to be sampled from the rainfall depths.  |
| ```number_of_main_divisions``` | The number of stratifications (main divisions) to use for the rainfall sampling. Divisions are spaced equally in standard normal variate scale. |
| ```number_of_sub_divisions``` | Number of samples to make within each stratification (main division). See below for the sampling method. |
| ```sample_method``` | Sampling method to use within each main division. There are three options: **1)** *normally distributed* uses a truncated normal distribution with mean of zero and standard deviation of one to sample standard normal variates within each main division -- this is the recommended method; **2)** *uniformly distributed* uses a uniform distribution for sampling in the standard normal variate scale; **3)** *stratified* uses an equally spaced sampling of standard normal variates across the main division.|
| ```number_of_temporal_patterns``` | The number of storm temporal patterns (usually 10) to sample from. |

**Table 3: Keys used in the ```tpt_quantile_analysis``` key**
| Key | Description |
| ----------- | ----------- |
| ```inflow``` | Range of dam inflow quantiles (m³/s) to estimate the AEP for, with the bounds and divisions of the range specified using the  ```lower```, ```upper```, and ```step``` keys. For example, ```"inflow": {"lower":  1000, "upper":  10000, "step":  1000}```|
| ```level``` | Range of lake level quantiles (m AHD) to estimate the AEP for, with the bounds and divisions of the range specified using the  ```lower```, ```upper```, and ```step``` keys. |
| ```outflow``` | Range of dam outflow quantiles (m³/s) to estimate the AEP for, with the bounds and divisions of the range specified using the  ```lower```, ```upper```, and ```step``` keys. |
