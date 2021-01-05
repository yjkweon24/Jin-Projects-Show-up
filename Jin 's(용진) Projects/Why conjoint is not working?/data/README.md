Data Files
----------
* profiles.tsv
* comparsions.tsv
* selections.tsv

File Format
-----------
profiles.tsv
    This file contains all the different profiles of that will be presented. The column names should be 'Profile' and the paramters. For instance, the tea unit test would have the column names 'Profile', 'Milk', and 'Sugar.' Each of the rows will simply be the profile number and the attribute for that profile. For instance, for the second profile of the computers unit test, the row values would be '2', 'Low', 'Low', 'Med.'

comparisons.tsv
    This file contains all the different comparisons that will be presented. The column names should be 'Comparisons' followed by 'Profile{}'.format(i) for i in [1, # of profiles]. These column names are not important as they will be dropped but allows the data to be more consistent through the different files. The row values are simply the comparison number and the profile number for each profile presented.

selections.tsv
    This file contains all the selections chosen. The column names are "Comparison" and "Individual{}".format(i) for i in [1, # of individuals]. These column names, similar to the column names in comparison.tsv, will be dropped. However, the consistency allows for easy of understanding across different data sets. The row values are simply the comparison number and the profile number **based off of the profile number that was in the column name in comparisons.tsv.** Thus, the values should not exceed the number of profiles that are presented in each comparison.
