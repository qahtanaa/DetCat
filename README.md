# DetectiveCat

Poor data quality, commonly caused by data glitches such as outliers, significantly affects any data analytics task, leading to inaccurate decisions and poor predictions of machine learning models. While detecting outliers in numerical data has been extensively studied, few attempts were made to solve the problem of detecting categorical outliers. Current industry standards are often based on organizationally defined, manual rules instead of being systemic. In this project, we developed a tool that can detect categorical outliers in relational datasets by identifying the syntactic structure of the values or the patterns in a given
column and detecting data elements that do not conform with the majority of the values. Recent studies aimed to solve the problem by developing frameworks that
extract regular expressions, however, this project developed more granular techniques on pattern identification. It focused on the trade-off between the expressiveness and
compactness of the generated patterns to obtain minimally descriptive patterns for a given attribute. A similarity metric, defining conformity of a value to a pattern
has been proposed and the tool was compared to existing manual and automated categorical error detection methods. 

In this document, we will present the steps necessary to use the system. In its current state, the tool is presented in the jupyter notebook "Detective_Cat". 
1. Load the jupyter notebook.
2. Enter the name of the dataset, along with the separator (default as comma), as well as the column in the main function:
       ![image](https://github.com/qahtanaa/DetectiveCat/assets/93333370/11c1b860-ef02-4def-84ad-41dd2ce97e33)
    ![image](https://github.com/qahtanaa/DetectiveCat/assets/93333370/8a807e5b-a110-471f-822f-62025a965e3e)
3. Execute the main() function cell. Enter 'Y' to run the tool with its default settings. Enter 'N' to input custom settings for the following values:
       A. limlength= As noted in Chapter 4, to improve the efficiency of the system, longer strings generate a limited set of patterns. The variable "limlength" determines at which length of strings the limited set is               generated. 
       B. maxlen= As noted in Chapter 4, very long strings only generate two patterns: their literal and fully abstracted representations. The variable "maxlen" determines wat which length of strings the binary patterns             are generated.
       C. genlevelcostmultiplier= As noted in Chapter 4, a score is assigned to each pattern based on their frequency, and number of characters abstracted (generalization level). The generalization level is multiplied by            a cost value. When set >1, the cost multiplier favors more expressive patterns. When set to <1, more generalized patterns will receive a higher final score. 
       D. minpatmultiplier= As noted in Chapter 4, patterns must cross a coverage percentage threshold in order to be aggregated and considered as valid values. This threshold is defined as the average coverage percentage           for all the patterns, multiplied by 0.7, as suggested by Jolliffe. This multiplier may be edited by the user. Values higher than 0.7 will increase the minimum threshold, potentially causing less                            patterns to be identified as dominant. 
       E. acc_threshold= As noted in Chapter 4, patterns will be aggregated until a certain percentage of coverage (explained variance) is acheived. This accumulation threshold is set at 0.95 (95%), but can be edited by             the user. Higher values mean the tool will potentially identify more patterns as dominant. A lower ceiling means the tool will potentially identify less patterns as dominant. Consider how likely you think                  outliers may be in your dataset when setting this value (higher ceiling means less outliers, lower means possibly more outliers are present).
       F. lencomparison= When dealing with a large number of patterns, the false positive and false negative detector engine will be unable to fully execute, given the number of crossproduct, pairwise computations. This             variable limits the number of pairwise comparisons done by only considering values whose length is less than the "lencomparison" number.
       G. pvalue= To determine whether a pattern is a false positive or false negative, the similarity scores are translated into z-scores, for which a one tailed statistical test is done to classify potential FP/FN                 values. The pvalue is set at 5% (0.05), however this can be edited by the user. A higher value will mean more patterns will be identified as potential FP/FN, while a lower will mean less. 
       H. smallnfpvalue= In the case that there are not enough pairwise computations for a statistical test to be used, the tool will determine whether a value is a false positive if its distance with a dominant pattern             is less than a value. The default value is 1.5, or the cost of one insertion of the same class character, but can be edited.
       I. smallnfnvalue= In the case that there are not enough pariwise computations for a statistical test to be used, the tool will determine whether a value is a false negative if its distance with a dominant pattern             is more than a value. The default value is 8.0, or the maximum substitution cost within the generalization tree, but can be edited. 
4. Results for the column are outpur as a dictionary that contains valid values, outlier values, possible false positives, and possible false negatives. A csv is also output with the same data. Additional information is      provided with the number of total patterns generated, the aggregation coverage results for the column, and the final pattern list. 



it's not even directed to me??? tf im not watching this. 5 times so 5 ppl wished u good luck?,,, i told u you're overprepared and too intense. so inspector is more of like just manual work. repetition. 

'this is for you' --> starts unbuttoning, then stops... 
