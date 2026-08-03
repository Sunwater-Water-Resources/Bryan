# Lake config file
This is the config file for managing antecedent lake levels and uses the keys in the table below.  

It is used by Monte Carlo simulations where the ```ADV``` key in the [simulation list](../sim_list.md) is set to *varying*, and by reservoir routing simulations where the ```ADV source``` key resamples the ADV from the ```lake_z``` column of an existing mcdf file.  

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
|```coefficients```|A dictionary of coefficients to use when applying the ```sigmoid``` method only. The coefficients include: ```k``` is the slope of the curve; ```Vf``` is the lower volume (ML) of the curve (floor);```Vc``` is the upper volume (ML) of the curve (ceiling); ```z0``` is the middle position of the curve in the horizontal (standard normal variate); ```H``` is a scaling factor for the curve height, which may be omitted for a standard logistic curve; ```A``` (optional, rarely needed) is the numerator of the sigmoid. **Read the section below before writing or reading a sigmoid curve**, because the value of ```H``` changes what ```Vc``` means and how ```A``` is derived from it. |
|```filename``` | Used for the ```empirical``` method only. The name of the csv file containing the user-defined (empirical) probability relationship. The csv file is assumed to be located in the same folder as the lake config file and must contain a column labeled ```z```, which is the standard normal variate, and ```volume```, which is the antecedent dam volume in ML. Note that Bryan interpolates the sampled *z* using a log transform of the volume. Therefore, zero or negative volumes will cause the code to crash. |

## The ```sigmoid``` curve and the ```A``` coefficient

Bryan evaluates the sigmoid antecedent lake volume curve as

```
log10 V(z) = A / (H + exp(-k (z - z0))) + log10(Vf)
```

so the curve runs from ```Vf``` at low *z* up to ```10 ** (A / H + log10(Vf))``` at high *z*.

**Two conventions for ```A``` are in circulation**, and ```A``` is not normally written in
the config file. **```H``` is what tells them apart**, so Bryan works ```A``` out from
```Vf``` and ```Vc``` -- which the user has already provided -- according to ```H```:

**Table 2a: the two sigmoid conventions**
| | Older models (```H``` > 1) | Newer models (```H``` = 1, or omitted) |
| --- | --- | --- |
| ```A``` | ```log10(Vc)``` | ```log10(Vc) - log10(Vf)``` |
| ```H``` | hand-tuned, typically 4 to 6 | ```1``` -- may be left out of the file altogether |
| ```Vc``` | **not the ceiling** -- it is ```10 ** A```. The true ceiling is a few percent away | **is** the ceiling, exactly |
| ```z0``` | wherever the hand-tuning put it | the exact midpoint between ```Vf``` and ```Vc``` in log space |
| Shape | ceiling and midpoint placed by hand | standard logistic, fitted |

Older models come from the original ADV assessment workbooks, where ```H``` -- the "height
adjustment factor" -- was turned by hand until the curve passed through a chosen midpoint
volume. That makes the upper asymptote equal ```Vc``` only when
```H = log10(Vc) / (log10(Vc) - log10(Vf))```, which is what the hand-tuning was converging
on. It is usually close but rarely exact.

Newer models fit a standard logistic between a chosen floor and ceiling, so ```H``` is 1 and
```Vc``` is genuinely the ceiling -- normally set to the full supply volume, so sampled
storage asymptotes to full supply and cannot pass it. Since ```H``` does nothing in this
case, it can be left out of the file.

> **Nothing needs changing to keep an older model running.** Its ```H``` is above 1, so it is
> read exactly as it always was and the results do not move.

> **```H``` = 1 can never be a hand-tuned curve**, which is what makes this safe to decide
> automatically. In the older convention ```H``` = 1 puts the upper asymptote at
> ```Vc * Vf```, which only means anything if the floor is 1 ML. Any realistic floor gives an
> ```H``` well above 1 -- a 10,000 ML floor under a 140,000 ML ceiling gives 4.6.

Bryan prints the resolved ```A```, which convention it came from, and both asymptotes of the
curve whenever it loads a sigmoid layer. **Check that line when picking up an unfamiliar
model**, because getting the convention wrong fails silently: the run completes, but the
antecedent volumes are absurd, most realisations start above the top of the storage curve,
and only the lake level frequency curve comes out blank -- the inflow curve still looks
perfectly sensible. Bryan also warns if the ceiling lands more than a factor of two above
```Vc```.

So a newer-style curve is written:

```json
{
  "exceedance_layer_info": [
    {
      "lower_z": -99,
      "upper_z": 99,
      "type": "sigmoid",
      "coefficients": {
        "k": 1.5931,
        "Vf": 12000.0,
        "Vc": 129041.0,
        "H": 1.0,
        "z0": -0.5646
      }
    }
  ],
  "volume_cap": "none"
}
```

giving ```A = log10(129041) - log10(12000) = 1.031546``` and a curve running from 12,000 ML
up to exactly 129,041 ML.

### Setting ```A``` by hand

```A``` may also be given in the ```coefficients``` explicitly, which overrides the rule
above. This is only needed for a curve that is neither convention; it is not required for
either of the cases in Table 2a, and normally it should be left out.

Note that ```volume_cap``` interacts with this. The ```ceiling``` option caps the sampled
volume at ```Vc```, which is a no-op for a newer-style curve (```Vc``` is already the
asymptote) but does bite on an older-style one, where the true ceiling sits above ```Vc```.

**Table 3: ```correlation_layer_info``` keys**
| Key | Description |
| ----------- | ----------- |
|```lower_z```| Lower bound of the AEP range for this probability band using standard normal variate scale. Use -99 if this band applies to the left-hand tail. |
|```upper_z```|Upper bound of the AEP range for this probability band using standard normal variate scale. Use 99 if this band applies to the right-hand tail.|
|```correlation```|Pearson correlation coefficient for lake standard normal variate versus rainfall standard normal variate for the current AEP band.|
