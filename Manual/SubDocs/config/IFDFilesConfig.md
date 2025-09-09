# IFD config file
This config file contains all the rainfall depth information and uses the keys listed in the table below.  
| Config file | Description |
| ----------- | ----------- |
| ```folder``` | Relative path to folder where IFD data is located.|
|```col_suffix```| The suffix in the headers used in the IFD tables (csv files). This will be stripped from the tables to extract the AEPs. | 
|```durations```| A list of dictionaries with information for each storm duration. Each dictionary has the following keys: ```duration``` to specify the duration in hours, and ```filename``` to specify the filename of the associated rainfall depth table (csv file). ++Note that the subcatchments order in the IFD files must be in the same order as the URBS subcatchment order.++|
|```PMP_depths```| The filename of the table of PMP depths (csv file)|
|```PMP_scaling```| The filename to the table of PMP spatial scaling factors for each sub-catchment (csv file)|
|```AEP_of_PMP```| To specify the AEP of the PMP. This was not automated because it was thought that there may be cases where specific AEPs will be used. However, in most cases: AEP = 1 in 10,000,000 if area < 100 km²; AEP = 1 in 10,000 if area > 100,000 km²; else AEP = 1 in 1 / 10^\[log(area) - 9].|
|```extreme_rainfall_interpolation_method```| The interpolation method to use for extreme rainfall. There are three options: **1)** The parabolic method provided in ARR is ```SiriwardenaWeinmann1998```; or **2)** a variation of this, using a normal scale for the AEP axis, is ```HillAndOthers2000```; **3)** It is recommended that ```GEV``` is used as this will deal with vertex perplexity that sometimes arises with the parabolic methods. |
|```pmp_exceedance_sampling```|(Optional) This key is unlikely to be required. It provides the user an option to deviate from the standard sampling method whereby rainfall for AEPs rarer than the AEP of the PMP are extrapolated from the curve fitted to the extreme rainfalls (```extrapolate```). While  not recommended, the rainfall sampled beyond the AEP of the PMP can be capped to match the PMP depth by setting this key to: ```cap```.|
