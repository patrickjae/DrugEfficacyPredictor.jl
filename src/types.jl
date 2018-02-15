""" An abstract type for samples, can be cell line or human tissue """
abstract type Sample end

""" A cell line comprised with its ID. """
struct CellLine <: Sample
	id::String
	cancer_type::String
end

"""
A data view, typed with a identifier and data type.
The identifier is mostly of type Gene, for RPPA data it is Protein.
In general, similarity matrices for samples are computed on the basis of
these identifiers. In other words, for each data view, a cell line is represented
as a vector of genes (or proteins) with associated values. The identifiers determine
the ordering of this vector.
"""

struct DataView{T, D}
	cell_line::CellLine
	measurements::Dict{T, D}
	DataView{T,D}(cl::CellLine) where {T,D} = new(cl, Dict{T,D}())
end

""" A gene. Mostly represented by its HGNC ID but with additional data if available."""
struct Gene
	entrez_id::String # <- primary key preferred
	hgnc_id::String
	ensembl_id::String
	cancer_gene::Bool
end

""" 
A protein. If available, paralogs of the protein that are affected by the same antibody are given.
The field antibody_validated indicates whether cross-reactivity in the assay has been checked
and the antibodies behave linearly in the dilution series.
"""
struct Protein
	hgnc_id::String
	antibody_validated::Bool
end

""" A gene expression."""
struct GeneExpression 
	expression_value::Float64
end

""" 
Illumina beta values indicating methylation.
CGct1 gives the number of CG dinucleotides in the probe sequence. The
U and M probes differ only at C's in CpG dinucleotides so there is a
good deal of cross-hyb between the two probes when the CG count is
low, and beta values shrink toward .5.  Probes are typically filtered
out with fewer than 3 CpGs as these values are unreliable.
Particularly, these probes will be less accurate for discrete M/U
calls.

Cct1 gives the number of off-CpG cytosines found in the probe. These
probes guarantee the specificity of the assay for successfully
bi-sulfite converted DNA, and probes with too little cytosine may
misinterpret unconverted DNA as methylated DNA, increasing beta
values. Probes are often filtered out that have fewer than 3
cytosines.
"""

struct Methylation 
	cgct1::Int64
	cct1::Int64
	illumina_beta_value::Float64
	is_methylated::Bool
	illumina_id::String
	function Methylation(beta_value::Float64, cgct1::Int64, cct1::Int64, methylated_threshold::Float64=.2; illumina_id::String="")
		new(cgct1, cct1, beta_value, beta_value≥threshold, illumina_id)
	end
end

"""
Whole transcriptome shotgun sequencing (RNA-seq) data.
Expression value is a fpkm value.
expression_status is 1 if the Ensembl gene model was detected above 
the background noise level.
"""
struct RNAseq 
	expression_value::Float64
	expression_status::Bool
end


"""
Whole exome sequencing data.

The different parts are:

num_cosmic - Number of samples in COSMIC that have mutation at this position.
variant_effect - Variant effect (silent, missense, nonsense, frame-shifiting Indel, etc.)
protein_change - Protein change
nucleotid_change - Change in DNAseq
variant_confidence - Confidence in variant (normalized in [0,1])
norm_zygosity - Zygosity (hom, het) of matched normal sample (blank)
norm_reference_count - # Reference alleles at position in matched normal (blank)
norm_variant_count - # Alternate alleles at position in matched normal (blank)
tumor_zygosity - Zygosity (hom, het) of Cell Line
tumor_reference_count - # Reference alleles at position in cell line
tumor_variant_count - # Alternate alleles at position in cell line
reference_mismatch_avg - METRIC: Average number of other mismatching bases in reads with reference base
variant_mismatch_avg - METRIC: Average number of other mismatching bases in reads with variant 
reference_mismatch_sum - METRIC: Base quality sum of mismatching bases in reads with reference base
variant_mismatch_sum - METRIC: Base quality sum of mismatching bases in reads with variant
reference_dist3effective_avg - METRIC: Average normalized distance of reference bases from 3' end of their respective reads
variant_dist3effective_avg - METRIC: Average normalized distance of variant bases from 3' end of their respective reads
details - Various other information collected at this position. 
"""
struct ExomeSeq 
	num_cosmic::Int64
	variant_effect::String
	protein_change::String
	nucleotid_change::String
	variant_confidence::Float64
	norm_zygosity::String
	norm_reference_count::Int64
	norm_variant_count::Int64
	tumor_zygosity::String
	tumor_reference_count::Int64
	tumor_variant_count::Int64
	reference_mismatch_avg::Float64
	variant_mismatch_avg::Float64
	reference_mismatch_sum::Float64
	variant_mismatch_sum::Float64
	reference_dist3effective_avg::Float64
	variant_dist3effective_avg::Float64
	details::String
end

""" Protein abundance value. """
struct RPPA 
	protein_abundance::Float64
end

""" The quantified value for DNA copy number analysis (aggregated across whole genes) """
struct CNV 
	gene_level_cnv::Float64
end

""" 
A drug, identified by a name.
Additional information may include knowledge about affected genes to
refine the number of features for predictions.
One useful approach could be to look at the generic pathways and select only those 
that are in a known pathway.
More specific would be to select only genes of the pathway that is affected by a drug.

The chemical structure as a string can be used to compute chemical drug similarity
to additionally aid inclusion of correlations between drugs.
"""
struct Drug
	id::String
	affected_genes::Vector{Gene}
	chemical_structure::String
end

"""
An experimental outcome.
An outcome should be associated with a drug.
The dictionary of outcome values associates the results to cell lines. 
The reason for this is, that we believe the outcome across different
cell lines is correlated according to the molecular structure of the cell line tissue.
In other words, we model the outcome on available cell lines jointly.
"""
struct Outcome
	outcome_values::Dict{CellLine, Float64}
	outcome_type::String
	Outcome() = new(Dict{CellLine, Float64}())
end

""" 
An experiment comprised of measurements (where each is associated 
with a cell line) and results.
"""
struct Experiment{T,D}
	input::Set{DataView{T,D}}
	results::Dict{Drug, Outcome}
	cell_lines::Vector{CellLine}
	genes::Vector{Gene}
	proteins::Vector{}

	function Experiment{T,D}() where {T,D}
		new(Set{DataView{T,D}}(), Dict{Drug, Outcome}(), Vector{CellLine}(), Vector{Gene}(), Vector{Protein}())
	end
end
