~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_Drug_Response_Training.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are a total of 53 cell lines and 31 drugs in the DREAM7 drug
sensitivity challenge 1 dataset.  DREAM7 participants are supplied
with drug response data on 35 cell lines and all 31 drugs.
Participants are challenged to predict the rank order of drug response
on the remaining 18 cell lines for all 31 drugs.

This file contains the training drug response data in the form of GI50
concentrations (Growth Inhibition of 50%).  The GI50 is a measure to
assess the efficacy of therapeutic compounds in a cell culture. To
estimate the GI50 a series of assays were performed as previously
described (Kuo, et al 2009).  Briefly, cells were treated for 72 hours
with a set of 9 doses of each compound in 1:5 serial dilution.  Cell
viability was determined using the Cell Titer Glo assay.  We used
nonlinear least squares to fit the data with a Gompertz curve. The
fitted curve was transformed into a GI curve using the method
described in (http://dtp.nci.nih.gov/branches/btb/ivclsp.html) and
previously described (Monks, et al, 1991).  In cases where the
underlying growth data are of high quality, but the GI50 was not
reached, the values were set to the highest concentration tested.  The
values reported are -log10 converted concentration values.

The drug response data was filtered to meet the following criteria: 1)
median standard deviation across the 9 triplicate datapoints < 0.20;
2) DT +/- 2SD of the median DT for a particular cell line; 3) slope of
the fitted curve > 0.25; 4) growth inhibition at the maximum
concentration < 50% for datasets with no clear response. For technical
or quality control reasons, please note that NOT all the cell lines
have been screened for all drugs. These cases are represented by an
"NA"

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kuo, W. L. et al. A systems analysis of the chemosensitivity of breast
cancer cells to the polyamine analogue PG-11047. BMC Med 7, 77,
(2009).

Monks, A. et al. Feasibility of a high-flux anticancer drug screen
using a diverse panel of cultured human tumor cell lines. J Natl
Cancer Inst 83, 757-766 (1991).


- Created by Laura Heiser
- Edited by Jim Costello for the DREAM7 Drug Sensitivity challenge 6/12/12
