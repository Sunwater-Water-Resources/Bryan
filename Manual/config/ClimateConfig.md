# Climate change config file
This is the config file for managing climate change adjustments according to the draft update to ARR released in late-2023. The main aspects of this file that will need to be checked are the catchment-specific Natural Resource Management Regions cluster and the Köppen-Geiger climate classification. It is unlikely other components of this config file will need to be changed often; perhaps when new IPCC reports and data are released or the ARR/Sunwater guidelines change. The keys used in this config file are dicussed in the tables below.  

**Table 1: Required keys**
| Key | Description |
| ----------- | ----------- |
|```NRM cluster```|The Natural Resource Management Regions cluster for the study catchment. A shapefile of these clusters is provided in the *climate_change/GIS folder*. |
|```KG classification```| The Köppen-Geiger climate classification for the study catchment. A raster of these clusters is provided in the *climate_change/GIS folder*, as well as a PDF world map for the legend. |
|```loss rates file```| Filename for the csv file with scaling factors for the initial and continuing loss rate adjustments, as specified in the draft ARR update. These data are used to compute the loss uplift factors for a given temperature rise. |
|```temporal pattern scaling```| Filename for the csv file with scaling factors for the D50 of temporal pattern adjustments, as specified in the draft ARR update. These data are used to compute the D50 shift for a given temperature rise. |
|```weighting files```| A dictionary pointing to the filepaths of weighting factors used in the probability-weighted sampling of temporal patterns to shift the average D50 of temporal patterns forward. The weightings are specified for the ```ARR areal```, ```ARR point```, ```GSDM```, and ```GTSMR``` temporal patterns.|
|```cap rainfall scaling```|Apply a cap to the percentage increase in rainfall per °C. For example, for a short duration like 1 hours, use an 8% per °C increase in rainfall instead of the 15% per °C prescribed in ARR using the following key pair: ```cap rainfall scaling: 8```. |

**Table 2: Optional keys if using the *SSP* method in the simulations list. However, this method is not in accordance with the Sunwater Climate Change Policy (Sunwater uses the *GWL* method)**.
| Key | Description |
| ----------- | ----------- |
|```base period``` | This key is only needed if using the SSP method with the ```SSP``` and ```Year``` keys provided in the simulation list. This is the median period of the rainfall records used in the rainfall frequency analysis that underlies the design rainfall estimates. For the 2016 revision of the IFDs provided by the Bureau, this is estimated at 1961 to 1990. Insert the period as a list; e.g. ```[1961, 1990]```|
|```files``` | This key is only needed if using the SSP method with the ```SSP``` and ```Year``` keys provided in the simulation list. This is a dictionary of files pointing to time-varying temperature data for both historic records (```historical```) and future projections (```SSP1-1.9```, ```SSP1-2.6```, ```SSP2-4.5```, ```SSP3-7.0```, ```SSP5-1.9```). The items in the dictionary are filepaths to the data (csv files) |

