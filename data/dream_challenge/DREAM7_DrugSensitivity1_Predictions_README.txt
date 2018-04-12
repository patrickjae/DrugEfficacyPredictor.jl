~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_Predictions.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the template file that participants will complete and is the
submission format for this challenge.  The file is organized as follows:

---------------------------------------------------------
| DrugAnonID   |      Drug1      |     Drug2      | ... |
---------------------------------------------------------
|  Cell Line1  |  rank of Cell   | rank of Cell   | ... |
|    	       |  Line1 on Drug1 | Line1 on Drug2 |     |  
---------------------------------------------------------
---------------------------------------------------------
|  Cell Line2  |  rank of Cell   | rank of Cell   | ... |
|    	       |  Line2 on Drug1 | Line2 on Drug2 |     |  
---------------------------------------------------------
...						        |
---------------------------------------------------------

Please note that the lowest ranks (1,2,3,..) correspond to most
efficacious compounds (lower GI50) to inhibit growth. The highest
ranks (..., 29, 30, 31) correspond to the weakest compounds in
inhibiting growth.


The anonymized drugs  are listed as columns and the cell lines are listed
as the rows.  Particpants are to replace the placeholder rank with
the predicted rank of Cell Line X on Drug Y.  Keep the header names
and row names. Save the predictions as a comma-separated-value (.csv)
file.  Rename the file with your <TeamName> as follows:
       DREAM7_DrugSensitivity1_Predictions_<TeamName>.csv

Note that participants are given drug response values for the cell
lines in the training set and thus the rank order of these cell lines.
The test set cell lines are to be rank ordered in relation to the cell
lines in the training set.  As there are 35 cell lines in the training
set and 18 cell lines in the test set, for submission, there will be a
rank ordered list of 53 cell lines for every drug.


-- Created by Jim Costello for the DREAM7 Drug Sensitivity challenge on 6/12/12.
