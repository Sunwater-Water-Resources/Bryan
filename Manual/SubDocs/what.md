# What Bryan does
The code has been set up to run URBS and RORB models through either a monte carlo or ensemble scheme. The sampling for the monte carlo scehme includes:
- Stratified **sampling of rainfall** according to the total probability theorem (TPT) method. The user must specify the number of divisions (e.g. 50 divisions) for the stratification along the annual exceedance probability axis (AEP); equally spaced in standard normal variate scale. By default, sampling within each division is done randomly according to a truncated normal distribution. Though a random uniform and equally spaced sampling can also be specified within the Python code. Note that all **AEP** inputs and outputs are in **1 in XXX** format. 
- **Temporal patterns** are sampled as integers between 0 and 9.
- **Storm types** are sampled between ARR point, ARR areal, GSDM, and GTSMR. The sampling is implemented when transitioning from rare to extreme events - randomly selecting between ARR and GSDM/GTSMR. 
- **Initial losses** are radomply scaled according to empirical data in ARR.
- **Preburst** percentiles are sampled to determine the proportion of the catchment average rainfall. Note that this is done after climate change uplift - so preburst is also adjusted for climate change. 
- **Residual preburst**, after subtracting losses, are included using a preburst temporal pattern database. The database was developed using storms extracted from pluvio records for gauges in the dam catchments. 
- **Antecedent lake levels** can be sampled according to an empirical distribution provided by the user. 
- **Embedded bursts** can be filtered.
- **Climate change** impacts on storm intensity, rainfall losses, and temporal patterns are included as per the *Draft update to the Climate Change Considerations chapter in Australian Rainfall and Runoff: A Guide to Flood Estimation* published in late-2023 by a technical working group led by Conrad Wasko at the University of Melbourne. 

For monte carlo simulations, the code can also be used to analyse results by applying the TPT method to create tables of the best estimates of peak inflow, lake level, and outflow. Plots of the results are also produced. For the ensemble simulations the analysis produces tables of the critical inflows, levels, and outflows with associated temporal patterns and durations. Box plots are also produced.  