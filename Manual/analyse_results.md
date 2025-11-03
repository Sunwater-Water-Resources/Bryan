# Result analysis
Once the simulations have been run, the results can be analysed setting the ```Analyse results``` key in the [simulation list](sim_list.md) to *yes*. At this point you may want to set the ```Run models``` key in the [simulation list](sim_list.md) to *no* to avoid rerunning the simulations. 

For the Monte Carlo simulations, the total probability theorem is used to extract the best estimate of the annual exceedance probability for a list of quantiles: inflow, lake level, and outflow. A csv file is created with the outputs as well as plots of the results. 

For the ensemble method, the median (Rank 6) of the inflow, lake level, and outflow is identified for each AEP and storm duration, and then the critical storm duration is identified for the inflow, lake level, and outflow. Again, csv files are created of the outputs and box plots are created of the analysis. 
