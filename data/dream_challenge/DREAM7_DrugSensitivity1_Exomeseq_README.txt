~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DREAM7_DrugSensitivity_Exomeseq.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Exome sequencing is a method that aims to capture sequence variation
in the form of mutations and indels (insertion/deletions) in all
coding regions of the genome.  This approach requires a hybridization
step to capture all exons and s sequencing step to the identify
sequence variation in the exons.  This procedure was done using the
Agilent SureSelect system.  Details of the technology can be found on
the Agilent website: http://www.genomics.agilent.com/. For data
processing, avery detailed list of steps to process the data is
described below, but in practice, DREAM7 participants are supplied
with a datafile that contains the following columns:

CellLine		      - Cell line name
Tissue			      - Origin of variant (Somatic or Germline). Assumed to be "Germline" for cell lines due to lack of matched normal. 
dbSNP			      - dbSNP ID, if variant overlaps
CancerGene?		      - Is the overlapping gene found in Sanger's Cancer Gene Census database.
#Cosmic			      - Number of samples in COSMIC that have mutation at this position.
Type		  	      - Variant effect (silent, missense, nonsense, frame-shifiting Indel, etc.)
HGNC_ID		 	      - Overlapping gene
Summary		  	      - Protein change
Chromosome	  	      - Chromosome
Start		    	      - Start position of variant
Stop		      	      - Stop position of variant
RefBase		      	      - Reference base
AltBase		      	      - Alternate / variant base
Confidence	      	      - Confidence in variant 
Zygosity(norm)	      	      - Zygosity (hom, het) of matched normal sample (blank)
RefCount(norm)	      	      - # Reference alleles at position in matched normal (blank)
AltCount(norm)	      	      - # Alternate alleles at position in matched normal (blank)
Zygosity(tumor)      	      - Zygosity (hom, het) of Cell Line
RefCount(tumor)      	      - # Reference alleles at position in cell line
AltCount(tumor)      	      - # Alternate alleles at position in cell line
Avg#Mismatch(ref)  	      - METRIC: Average number of other mismatching bases in reads with reference base
Avg#Mismatch(alt)  	      - METRIC: Average number of other mismatching bases in reads with variant 
MismatchQualitySum(ref)	      - METRIC: Base quality sum of mismatching bases in reads with reference base
MismatchQualitySum(alt)	      - METRIC: Base quality sum of mismatching bases in reads with variant
DistanceEffective3'End(ref)   - METRIC: Average normalized distance of reference bases from 3' end of their respective reads
DistanceEffective3'End(alt)   - METRIC: Average normalized distance of variant bases from 3' end of their respective reads
Details	 		      - Various other information collected at this position. 



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Detailed methods description
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Five3 Alignment + Mutation Calling Pipeline

(1) Fastq Files -> Initial BAM 
    a) Pairs of Fastq files (i.e. R1 & R2) sequenced from the same sample
are aligned separately using bwa aln & bwa sampe (default parameters)
to the hg19 (GRCh37) reference. Each pair of Fastq files generates a
single BAM file.
     b) The individual BAM files from (a) are merged and sorted by
chromosomal position to generate a single BAM file representing all
reads from the sequencing run. If BAM files originated from different
libraries, this information is stored and tracked in the merged BAM
file.
     c) The BAM file from (b) is then scanned for duplicate reads
using a program called superDeDuper. All read pairs with identical
start / stop coordinates are collected into a group. From this group
of duplicate reads, the read pair with the highest average base /
mapping quality is chosen as the representative read, while all others
are marked as duplicates in each of their respective BAM records. This
is performed for each library separately, and results in the Initial
BAM with duplicates marked, but not removed.

(2) Initial BAM -> Indel Realignment BAM
    a) Using the GATK routine RealignerTargetCreator, the initial BAM
from step (1) is scanned for genomic regions that may require indel
realignment using locations of known indels from the 1,000 genomes
project. A list of such regions is produced.
    b) Using the GATK routine IndelRealigner, misaligned reads in
regions identified in (a) are realigned using the knowledge that an
indel is known to be nearby. Reads that are located such that one of
their ends overlap an indel are often aligned incorrectly (i.e. indel
is not present in read's alignment). This can lead to a series of
mismatching bases where the indel is located. Since similarly
mismatched bases can occur in multiple reads due, this can quickly
lead to false positive mutation calls. This step will realign those
reads that see an improved mapping quality upon introduction of the
indel at the appropriate location in the read. A BAM file with indel
realignments is generated.


(3) Indel Realignment BAM -> Base Quality Recalibrated BAM ("Final BAM")
    a) Using the GATK routine CountCovariates, the BAM file from step
(2) is analyzed to generate the covariates necessary to perform base
quality recalibration.  Briefly, it searches for mismatching bases in
reads that do NOT overlap known heterozygous sites (1,000 genomes +
dbSNP) and collects information on the mismatching base?s quality and
a series of other covariates (e.g. base quality, read group,
neighboring bases, sequencing cycle). This is performed across the
whole BAM file, down-sampling reads to a coverage equivalent to 1X,
and generates a file containing the recalibration covariates.
    b) Using the GATK routine TableRecalibration, the recalibration
metrics from step (a) are used to recalibrate all base qualities from
the Indel Realignment BAM file generated in step (2). This step is
necessary as the base qualities generated by the sequencer often
inaccurately reflect the true frequency of mismatching bases.  For
instance, a Phred-scaled quality of 40 means that out of 10,000 bases
with Q40, only one base will be incorrectly called. If step (a) finds
that bases with Q40 are observed to have a higher miscall frequency
than 1 out of 10,000, those bases will have their base qualities
reduced accordingly. A BAM file with base quality recalibration is
generated.
    c) The BAM file generated by step (b) is the Final BAM file used
in all post- processing steps below.


(4) Final BAM -> Quality Metrics
    a) The flagstats routine of samtools is run on the BAM file,
counting total reads, duplicate reads, unmapped reads, unpaired reads,
and read pairs where the reads are aligned to different chromosomes.
    b) The number of reads & base coverage is computed for each exon,
intron, and intergenic region. The exon/intron/intergenic regions are
derived from the "Canonical" isoforms of UCSC Known Genes.
    c) Average base quality vs. base position within read is
    generated.
    d) Picard's routine EstimateLibraryComplexity is run to infer the
size/complexity of each input library.


After steps (1-4) are completed for both tumor and normal sequences
for a single patient, the steps below are performed. For all steps
below, reads overlapping the same genomic position in both tumor &
normal BAM files are collected and presented to each algorithm.


(5) Tumor vs. Normal Mutation Calling
    a) Allele counts and their associated base qualities from tumor
and normal are collected. Alleles used in the subsequent steps must
meet the following criteria:
	- Base quality (BQ) >= 10
	- Neighborhood Base Quality (NBQ) >= 10
	- Mapping Quality of associated read (MQ) >= 20
	- Read is not marked as duplicate
Any base quality exceeding the read's mapping quality is reduced to
the read?s mapping quality.
    b) If there are less than 2 reads supporting any non-reference
allele at the current position, then that position is deemed
homozygous reference and no further analysis is performed. Otherwise,
move on the steps below.
    c) The likelihoods of all possible genotypes (AA, AT, AC, etc.)
given the allelic data collected in (a) for both tumor and normal are
computed using the MAQ error model originally defined in (Li, et al.,
2008) and now available in the samtools source code.
    d) The genotype likelihoods from (c) are used in a Bayesian model
incorporating a prior probability on the reference, the heterozygous
rate of the human genome (1 in a 1,000 positions is heterozygous), the
probability of the required mutations to convert the normal genotype
to the tumor genotype, and a baseline estimate of normal contamination
in the tumor sample (10%). Each tumor+normal genotype pair is scored
using this model. The genotype pair with the highest likelihood given
the data is chosen as the most likely tumor and normal genotypes.

    e) If the tumor and normal genotypes chosen in step (d) indicates
any deviation from reference, a number of metrics are computed:
	DP: Total depth in each sample.
	AD: Depth or coverage for all alleles, including alleles not in genotype.
	BQ: Average base quality of each allele.
	MQ: Average mapping quality of reads supporting each allele. 
	MQ0: Number of mapping quality zero reads overlapping position
	MQL: Number of ?low? mapping quality reads overlapping position, low is defined as having a mapping quality between 1 and 20
	NAHP: Average number of adjacent homopolymer runs on either side of each allele in genotype
	MAHP: Longest adjacent homopolymer run on either side of each allele in genotype
	AMM: Average number of mismatches in reads supporting each allele
	MMQS: Average sum of the base qualities for all mismatching bases in reads supporting each allele
	DETP: Average effective distance to 3? end of read for each allele. If there 
	      is a run of low quality bases (BQ <= 2) at the end of the read, the 
	      "effective" end of the read is placed at the base preceding the run of low 
	      quality bases. Otherwise, the end is placed at the stop position of the 
	      read?s alignment. Distance to read?s end is normalized by the read?s 
	      length such that 0.0 is set at the beginning of the read, and 1.0 is the 
	      effective end of the read. 
	LD/MD/RD: Number of reads supporting each allele where the allele is 
		  located in the left-most third of read, middle-third of read, or right-most 
		  third of read, respectively.
	LDS/MDS/RDS: Strand-aware version of above, where left-most and 
		     right-most counts are swapped depending on the strand (forward or 
		     reverse) of each read. 
	SB: Number of reads supporting each allele aligned to the forward strand. 
	PN/NN: Previous and next nucleotides in reference, for determine 
	       genomic context of any mutation.
These metrics are used for post-processing filtering of all putative variants. 
    f) The paired genotypes from (d) also determine the variant type
at each position. If tumor and normal genotypes are identical, then
the variant is germline. Germline variants can either be homozygous
non-reference or heterozygous. If the normal genotype is heterozygous
and the tumor genotype is homozygous, and the tumor allele is in the
normal genotype, then the variant is classified as LOH (loss of
heterozygosity). Otherwise, the variant is classified as a somatic
mutation.
    g) If the variant is classified as germline or LOH in (f), then
the log-likelihood of the paired genotype is used to compute a
Phred-scaled quality / confidence of the germline variant.
    h) If the variant is classified as somatic in (f), a Somatic Score
(SS) is also calculated and used in place of the genotype quality
computed in (g). The Somatic Score was originally defined in (Larson,
et al. 2012) and determines the likelihood that the position is not
somatic by summing the likelihoods of all non- somatic genotype pairs
computed in step (c). Like the genotype quality in (g), SS is
converted to a Phred-scale.
   i) All putative variants (germline, LOH, and somatic) and
associated metrics are converted to the VCF format. When converting to
VCF, the following filters are applied to each variant:
	conf: Genotype quality / Somatic Score >= 100
	dp: Total depth (DP of normal + primary) >= 8
	mdp: Maximum depth (DP of primary + normal) < 800
	mq0: MQ0 < 5
	mql: MQL < 5 
	sb: Mutant allele strand bias p-value > 0.005, using a Binomial test on the SB metric collected
	mmqs: MMQS <= 20 
	amm: AMM <= 1.5
	detp: 0.2 <= DETP <= 0.8
	ad: AD of mutant allele in tumor >= 4
	gad: AD of mutant allele in normal <= 3
	ma: Maximum of two alleles have read support >= 2
Variants that pass all filters will be marked PASS in the FILTER
column of their VCF record. Otherwise, the names of each above filter
that the variant does not meet will be recorded in the FILTER column.


Single Sample Analysis

After the alignment pipeline has completed the alignments for a sample that has no 
matched-normal control (e.g. cell line), the steps below are performed.
(1) Mutation Calling

    a) Allele counts and their associated base qualities from the
sample are collected.  Alleles used in the subsequent steps must meet
the following criteria:
     - Base quality (BQ) >= 10
     - Neighborhood Base Quality (NBQ) >= 10
     - Mapping Quality of associated read (MQ) >= 20
     - Its associated read is not a duplicate
Any base quality exceeding the read?s mapping quality is reduced to
the read's mapping quality.
    b) If there are less than 2 reads supporting any non-reference
allele at the current position, then that position is deemed
homozygous reference and no further analysis is performed. Otherwise,
move on the steps below.
    c) The likelihoods of all possible genotypes (AA, AT, AC, etc.)
given the allelic data collected in (a) for the sample are computed
using the MAQ error model originally defined in (Li, et al., 2008) and
now available in the samtools source code.
    d) The genotype likelihoods from (c) are used in a Bayesian model
incorporating a prior probability on the reference, and the
heterozygous rate of the human genome (~1 in a 1,000 positions is
heterozygous). The genotype with the highest likelihood given the data
is chosen as most likely.

    e) If the genotype chosen in step (d) is homozygous reference, no
further analysis is performed at this position. Otherwise, a number of
metrics are computed at the variant position:
  DP: Total read depth.
  AD: Depth or coverage for all alleles, including alleles not in genotype.
  BQ: Average base quality of each allele.
  MQ: Average mapping quality of reads supporting each allele.
  MQ0: Number of mapping quality zero reads overlapping position
  MQL: Number of "low" mapping quality reads overlapping position, low is defined as having a mapping quality between 1 and 20
  NAHP: Average number of adjacent homopolymer runs on either side of each allele in genotype
  MAHP: Longest adjacent homopolymer run on either side of each allele in genotype
  AMM: Average number of mismatches in reads supporting each allele
  MMQS: Average sum of the base qualities for all mismatching bases in reads supporting each allele
  DETP: Average effective distance to 3? end of read for each allele. If there 
  	is a run of low quality bases (BQ <= 2) at the end of the read, the 
	"effective" end of the read is placed at the base preceding the run of low 
	quality bases. Otherwise, the end is placed at the stop position of the 
	read's alignment. Distance to read?s end is normalized by the read?s 
	length such that 0.0 is set at the beginning of the read, and 1.0 is the 
	effective end of the read.
  LD/MD/RD: Number of reads supporting each allele where the allele is 
  	    located in the left-most third of read, middle-third of read, or right-most 
	    third of read, respectively.
  LDS/MDS/RDS: Strand-aware version of above, where left-most and 
  	       right-most counts are swapped depending on the strand (forward or 
	       reverse) of each read.
  SB: Number of reads supporting each allele aligned to the forward strand.
  PN/NN: Previous and next nucleotides in reference, for determine genomic context of any mutation.
These metrics are used for post-processing filtering of all putative variants.
    f) Since no normal control is available, all variants are
considered germline and the genotype's log-likelihood is used to
compute a Phred-scaled quality / confidence of the germline variant.
    g) All putative variants and associated metrics are converted to
the VCF format. When converting to VCF, the following filters are
applied to each variant:
  conf: Genotype quality >= 100
  dp: Total depth >= 8
  mdp: Maximum depth < 800
  mq0: MQ0 < 5
  mql: MQL < 5
  sb: Mutant allele strand bias p-value > 0.005, using a Binomial test on the SB metric collected
  mmqs: MMQS <= 20
  amm: AMM <= 1.5
  detp: 0.2 <= DETP <= 0.8
  ad: AD of mutant allele >= 4
  ma: More than two alleles have read support >= 2
Variants that pass all filters will be marked PASS in the FILTER
column of their VCF record. Otherwise, the names of each above filter
that the variant does not meet will be recorded in the FILTER column.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
References
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Li, H., Ruan, J. and Durbin, R. (2008) Mapping short DNA sequencing
reads and calling variants using mapping quality scores, Genome Res,
18, 1851-1858.

Larson, D.E., et al. (2012) SomaticSniper: Identification of Somatic
Point Mutations in Whole Genome Sequencing Data, Bioinformatics, 28,
311-317.
	
