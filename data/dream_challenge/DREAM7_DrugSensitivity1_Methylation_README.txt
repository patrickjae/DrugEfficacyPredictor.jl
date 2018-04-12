~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
DREAM7_DrugSensitivity1_Methylation.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
The methylation data has been compiled using the Illumina Human
Methylation Bead array.  The values reported are illumina beta values,
which are calculated as M/(M+U) where M and U are the intensities for
the probes measuring the abundance of methylated and unmethylated DNA
respectively, and so is interpreted as the proportion of molecules
that are methylated.

Should you wish to discretize the data as methylated/unmethylated
comparison to MSP suggests that beta=.2 is the dividing line for the
average, with any value between .1 and .2 giving similar results.

The datafile is a matrix of beta values for the listed cell lines with
the first 4 columns as follows:
-Illumina_ID
-HGNC_ID
-CGct1
-Cct1

CGct1 gives the number of CG dinucleotides in the probe sequence. The
U and M probes differ only at C's in CpG dinucleotides so there is a
good deal of cross-hyb between the two probes when the CG count is
low, and beta values shrink toward .5.  Probes are typically filtered
out with fewer than 3 CpGs as these values are unreliable.
Particularly, these probes will be less accurate for discrete M/U
calls.

Cct1 gives the number of off-CpG cytosines found in the probe. These
probes guarentee the specificity of the assay for successfully
bi-sulfite converted DNA, and probes with too little cytosine may
misinterpret unconverted DNA as methylated DNA, increasing beta
values. Probes are often filtered out that have fewer than 3
cytosines.


-- Created by Leslie Cope from the Sukumar lab  10/30/09
-- Edited by Jim Costello for the DREAM7 Drug Sensitivity challenge on 6/12/12.
