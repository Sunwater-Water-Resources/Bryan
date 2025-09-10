![Juniper berries](Manual/stock-photo-juniper-twig-with-berries_3.png)
# Bryan: Sunwater's design flood simulation manager
### A brief guide

**Developers**: Richard Sharpe (Sunwater) and Graigan Panosot (Jacobs).

Bryan is a Python-based platform that has been developed to implement Sunwater's design flood hydrology specification. 

Why Bryan? For a similar reason that Python is called Python - in reference to Monty Python. Early on, Monty Python was a slip in speech when speaking of monte carlo. And *The Life of Brian* is classic movie by the comedy troupe. But Bryan is more pythonic than Brian. So, there it is. 

- [What Bryan does](Manual/what.md)
- [Some pre-requisites](Manual/pre-requisites.md)
- [Installation](Manual/installation.md)
- [Getting started](Manual/getting_started.md)
- [Hydrologic modelling](Manual/hydrologic_modelling.md)
- [The simulations list](Manual/sim_list.md)
- [The config files](Manual/config_files.md)
- [Analysing the model results](Manual/analyse_results.md)

There are also a few Python scripts that can be used for post-processing. For more info ono this [click here](Manual/utilities.md).

Hang on there, Richard... why not just use Storm Injector of the URBS Control Centre for Monte Carlo TPT? Well, Storm Injector had not released Monte Carlo for URBS when we started this project and the URBS Manual did not seem to provide details of all the analyses going on under the hood. We wanted to be sure of those details and test some of our own methods. So, why not just develop our own URBS wrapper. For an example of some of the details we added see:

Criste-Jones, L., Watt, S., Sharpe, R. (2025). Improving the coupling of antecedent storage and storms for design flood estimation. *Hydrology and Water Resources Symposium*, Engineers Australia.
Panosot, G., Watt, S., Sharpe, R. (2025). Development and application of pre-burst temporal pattern ensembles for design flood hydrology. *Hydrology and Water Resources Symposium*, Engineers Australia.
Sharpe, R.G. (2024). Dealing with vertex perplexity in parabolic interpolation of design rainfall. *ANCOLD conference*, November 2024

Well, what about the juniper bushes over there?
