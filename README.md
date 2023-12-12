# Detective Cat (DetCat): Detecting Categorical Outliers in Relational Datasets

Poor data quality, commonly caused by data glitches such as outliers, significantly affects any data analytics task, leading to inaccurate decisions and poor predictions of the machine learning models. While detecting outliers in numerical data has been extensively studied, few attempts were made to solve the problem of detecting categorical outliers. In this repository, we provide DetCat for detecting categorical outliers in relational datasets by utilizing the syntactic structure of the values. For a given attribute, DetCat identifies a set of patterns that represents the majority of the values and declares the data values that cannot be generated by those patterns as outliers. Moreover, DetCat considers the trade-off between the expressiveness and the compactness of the generated patterns and introduces a new similarity measure between the patterns in order to reduce the false predictions.

## Running the tool
The tool is implemented in Python and uses dash library for the graphical interface. To run the tool, proceed with the following instructions:
* get the code first:
```
# clone using https
git clone https://github.com/qahtanaa/DetectiveCat.git

cd DetectiveCat
```
* Install the required libraries
```
pip install -r requirements.txt
```
#### Running from the command line

```
python detective_cat table_name min_coverage max_len min_cov_per_pattern [default|custom]
```

#### Running the GUI

```
python appDC.py 
```

## Parameters

Detcat has four main parameters: 

* minimum coverage of the dominant patterns ($\xi$): we add more patterns to the set of dominant patterns until the ratio of the generated values using the dominant patterns to the total number of values is greater than $\xi$; the rest of the patterns are flagged as potential erroneous patterns.
* The max allowed length for values ($L$): no exhaustive search for patterns in values with length greater than $L$; instead, we use only two patterns (the pattern of the original characters and the fully generalized pattern).
* The minimum coverage per dominant pattern ($\gamma$): for a pattern $p_i$ to be flagged as dominant pattern, the ratio of the generated values by $p_i$ to the total number of values in the attribute must be greater than $\gamma$; Parameters $\xi$ and $\gamma$ are considered based on the fact that outliers are rare events and should be minority in any given dataset.
* The similarity measure: we use the customized Levenshtein distance to measure the distance between patterns. However, the original Levenshtein distance can also be tested in filtering the false positives/negatives. In the command line, this can be set using `default` for original Levenshtein distance and `custom` for the customized Levenshtein distance. 

  
## Testing the tool

Download the [hospital dataset](https://github.com/BigDaMa/raha/tree/master/datasets/hospital) and use the following parameters:
Min coverage = 99
Max Length = 15
Min Cov. per Pattern = 0.1
Distance Function = Cust. Lev. distance
