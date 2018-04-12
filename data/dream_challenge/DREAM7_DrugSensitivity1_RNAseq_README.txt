~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data Summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Whole transcriptome shotgun sequencing (RNA-seq) was completed on
breast cancer cell lines and expression analysis was performed with
the ALEXA-seq software package as previously described (Griffith, et
al 2010). Briefly, this approach comprises (i) creation of a database
of expression and alternative expression sequence features (genes,
transcripts, exons, junctions, boundaries, introns, and intergenic
sequences) based on Ensembl gene models, (ii) mapping of short
paired-end sequence reads to these features by BLAST and BWA, (iii)
identification of features that are expressed above background noise
while taking into account locus-by-locus noise.

Please see Griffith, et al. for full details.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_RNAseq_quantification.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
As with gene expression micorarrays, an expression value can be
calculated with RNA seq data.  Accordingly, the fpkm values for
transcript expression has been calculated for each Ensembl gene model
over all the tested cell lines.  The first colum contains the Ensembl
gene id and the second column contains the HGNC gene id.  The
remaining colums represent the expression values for the listed cell
lines.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity1_RNAseq_expressed_calls.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

In addition to estimating expression values, the expression status of
the gene model was also calculated, where a binary call (1 or 0) was
made if the Ensembl gene model was detected above the background noise
level.  As with the previous file, the first colum contains the
Ensembl gene id and the second column contains the HGNC gene id.  The
remaining colums represent the expression status calls for the listed
cell lines.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Griffith, M. et al. Alternative expression analysis by RNA
sequencing. Nat Methods 7, 843-847 (2010).

-- File created by Laura Heiser on 5/1/2012
-- Edited by Jim Costello for the DREAM7 Drug Sensitivity challenge on 6/12/12.

