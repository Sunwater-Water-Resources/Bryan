# Lake config file
This is the config file for managing antecedent lake levels and uses the keys in the table below.  

**Table 1: Config file keys**
| Config file | Description |
| ----------- | ----------- |
| ```exceedance_layer_info``` | Here, a list (```[]```) or a single entry of probability distribution information for sampling the antecedent lake volume is provided. Using a list enables several probability distributions to be applied, with different relationships across different ranges of AEPs (provided as standard normal variates). If only one distribution is needed, the list will contain only one dictionary and the square brackets can be excluded. See Table 2 below for more information on the keys used in the dictionaries in the list. This file is only needed in Monte Carlo simulations where the ```ADV``` key in the simulation list has been specified as *varying*. |
|```correlation_layer_info```| Here, a list ```[]``` or a single entry of correlation information for sampling the antecedent lake volume is provided. Using a list enables several correlations to be applied across different ranges of AEPs (provided as standard normal variates). For example, if there is a correlation, it is usually only for the more extreme AEPs. Thus, the correlation can be set to zero for the more frequent AEPs. See Table 3 below for more information on the keys used in the dictionaries in the list. This key is only needed in Monte Carlo simulations where the ```ADV``` key in the simulation list has been specified as *varying* and a correlation between rainfall and antecedent lake volume has been identified.|
|```volume_cap```|(Optional) This key is used to handle the capping of the sampled ADV to FSV. By default, if the sampled ADV happens to be larger than the FSV, Bryan caps the ADV to the FSV. However, historical records might indicate that lake levels are at times above FSL prior to a large storm. Thus, this capping can be turned off by setting this key to ```none```.|

**Table 2: ```exceedance_layer_info``` keys**
| Key | Description |
| ----------- | ----------- |
|```lower_z```| Lower bound of the AEP range for this probability band using standard normal variate scale. Use -99 if this band applies to the left-hand tail. |
|```upper_z```|Upper bound of the AEP range for this probability band using standard normal variate scale. Use 99 if this band applies to the right-hand tail.|
|```type```| Can be one of three types of distributions: **1)** ```uniform``` applies a fixed antecedent lake level; **2)** ```sigmoid``` applies a sigmoidal curve based on a logistic function; **3)** The ```empirical``` method applies a user-defined distribution from a csv file -- see the ```filename``` key.|
|```value_ML```|The volume to use (in ML) when applying the ```uniform``` method.|
|```coefficients```|A dictionary of coefficients to use when applying the ```sigmoid``` method only. The coefficients include: ```k``` is the slope of the curve; ```Vf``` is the lower volume (ML) of the curve (floor);```Vc``` is the upper volume (ML) of the curve (ceiling); ```z0``` is the middle position of the curve in the horizontal (standard normal variate); ```H``` is a scaling factor for the curve height. |
|```filename``` | Used for the ```empirical``` method only. The name of the csv file containing the user-defined (empirical) probability relationship. The csv file is assumed to be located in the same folder as the lake config file and must contain a column labeled ```z```, which is the standard normal variate, and ```volume```, which is the antecedent dam volume in ML. Note that Bryan interpolates the sampled *z* using a log transform of the volume. Therefore, zero or negative volumes will cause the code to crash. |

**Table 3: ```correlation_layer_info``` keys**
| Key | Description |
| ----------- | ----------- |
|```lower_z```| Lower bound of the AEP range for this probability band using standard normal variate scale. Use -99 if this band applies to the left-hand tail. |
|```upper_z```|Upper bound of the AEP range for this probability band using standard normal variate scale. Use 99 if this band applies to the right-hand tail.|
|```correlation```|Pearson correlation coefficient for lake standard normal variate verus rainfall standard normal variate for the current AEP band.|
