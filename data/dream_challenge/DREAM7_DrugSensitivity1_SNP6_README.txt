~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data Summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SNP6.0 data on breast cancer cell lines generated by the Joe Gray Lab
at LBNL.  Data were generated in two batches (2007 and 2009).

SNP Array and DNA copy number analysis Genome copy number was assessed
using the Affymetrix Genome-Wide Human SNP Array 6.0 analysis
platform.  Arrays were analyzed using aroma.affymetrix (Bengtsson, et
al. 2008) (http://aroma-project.org), data were normalized as
described(Bengtsson, et al. 2009) and DNA copy number ratios at each
locus were estimated relative to a set of 20 normal sample arrays.
Data were segmented using circular binary segmentation (CBS) from the
bioconductor package DNAcopy(Venkatraman & Olshen, 2007).  We used
genome build HG18 for processing and annotating.

Raw data for the 2007 batch are available in The European Genotype
Archive (EGAS00000000059), and were published in conjunction with
Heiser, et al. 2012.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_SNP6.cbs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The datafile is a tab-delimited text file that contains segmented
genome copy number calls .  The column headers can be interpreted as
follows:

CellLine: cell line name
Chromosome: chromosome number
Start: start of segment (in bp)
Stop: end of segment (in bp)
Count: number of probes in the segment
Mean: mean copy number for the segment (log2 transformed)

This datafile can be visualized with the freely available software
package called the Integrative Genomics Viewer (IGV) developed by the
Broad Institute (http://www.broadinstitute.org/igv/download).

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity_SNP6_gene_level.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data from DREAM7_DrugSensitivity1_SNP6.cbs was processed using the R
package CNTools to identify gene-level changes in copy number.

The datafile is a tab-delimited text file that contains a gene by cell
line matrix of gene level copy number quantification.  The first
column in the file is an Entrez ID and the second column is the HGNC
id.  The remaining columns are individual cell lines.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bengtsson H, Irizarry R, Carvalho B, & Speed TP (2008) Estimation and
assessment of raw copy numbers at the single locus level.
Bioinformatics (Oxford, England) 24(6):759-767.

Bengtsson H, Wirapati P, & Speed TP (2009) A single-array
preprocessing method for estimating full-resolution raw copy numbers
from all Affymetrix genotyping arrays including GenomeWideSNP 5 & 6.
Bioinformatics 25(17):2149-2156.

Venkatraman ES & Olshen AB (2007) A faster circular binary
segmentation algorithm for the analysis of array CGH data.
Bioinformatics 23(6):657-663.


-- Created by Laura Heiser on 7/25/11.
-- Edited by Jim Costello for the DREAM7 Drug Sensitivity challenge on 6/12/12.
-- Copyright 2011 Lawrence Berkeley National Laboratory. All rights reserved.

