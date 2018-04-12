##############################################################################################################
##############################################################################################################
########################################         MODEL         ###############################################
##############################################################################################################
##############################################################################################################
mutable struct PredictionModel
	# precision parameters, T-dimensional, one per drug/task
	# ɣ::IsometricPrecisionParameter #precision to b <-- could be a matrix introducing covariance between drugs
	ɣ::Vector{GammaParameter}
	# drug specific scalars, T-dimensional
	b::Vector{NormalParameter}

	# λ::Vector{IsometricPrecisionParameter} #precision to a <-- introduce covariance between cell lines
	λ::Vector{Vector{GammaParameter}}
	# drug specific vectors, TxN_t dimensional
	a::Vector{Vector{NormalParameter}}

	# ε::IsometricPrecisionParameter #precision to y <-- between drugs
	ε::Vector{GammaParameter}

	# ν::IsometricPrecisionParameter #precision to g <-- between drugs and views?
	ν::Vector{GammaParameter}
	# intermediate results, N_t dimensional, one per drug and view
	G::Matrix{MvNormalParameter}

	# precision parameter, K-dimensional, one per view
	# ⍵::IsometricPrecisionParameter
	⍵::Vector{GammaParameter}
	# view specific scalars, K-dimensional
	e::Vector{NormalParameter}

	T::Int64
	K::Int64
	N::Vector{Int64}
	
	function PredictionModel(T::Int64, K::Int64, N::Vector{Int64})
		p = new()
		p.T = T
		p.K = K
		p.N = N
		# TODO pass in as parameters
		a_0 = 1e-3
		b_0 = 1e-3
		m_0 = 0.
		s_0 = 10.

		p.ɣ = VectorGammaParameter(T)
		p.b = Vector{NormalParameter}(T)
		for t in 1:T
			p.ɣ[t] = GammaParameter(a_0, b_0)
			p.b[t] = NormalParameter(m_0, s_0)
		end

		p.λ = Vector{VectorGammaParameter}(T)
		p.a = Vector{Vector{NormalParameter}}(T)
		for t in 1:T
			p.λ[t] = VectorGammaParameter(N[t])
			p.a[t] = MvNormalParameter(m_0, s_0, N[t])
			for n in 1:N[t]
				p.λ[t][n] = GammaParameter(a_0, b_0)
				# p.a[t][n] = NormalParameter(m_0, s_0)
			end
		end

		p.ε = VectorGammaParameter(T)
		[p.ε[t] = GammaParameter(a_0, b_0) for t in 1:T]

		p.ν = VectorGammaParameter(T)
		[p.ν[t] = GammaParameter(a_0, b_0) for t in 1:T]
		p.G = Matrix{MvNormalParameter}(T,K)
		for t in 1:T, k in 1:K
			p.G[t,k] = MvNormalParameter(m_0, s_0, N[t])
		end

		p.⍵ = VectorGammaParameter(K)
		p.e = Vector{NormalParameter}(K)
		for k in 1:K
			p.⍵[k] = GammaParameter(a_0, b_0)
			p.e[k] = NormalParameter(m_0, s_0)
		end
	end
end


##############################################################################################################
##############################################################################################################
########################################         DATA          ###############################################
##############################################################################################################
##############################################################################################################
""" An abstract type for samples, can be cell line or human tissue (latter not yet implemented) """
abstract type Sample end

abstract type ViewType end

"""
A data view, typed with a identifier and data type.
The identifier is mostly of type Gene, for RPPA data it is Protein.
In general, similarity matrices for samples are computed on the basis of
these identifiers. In other words, for each data view, a cell line is represented
as a vector of genes (or proteins) with associated values. The identifiers determine
the ordering of this vector.
"""

struct DataView{K, V <: ViewType}
	cell_line_id::String
	measurements::OrderedDict{K, V}
	DataView{K,V}(cell_line_id::String) where {K,V} = new(cell_line_id, OrderedDict{K,V}())
end

""" 
A cell line comprised with its ID, information on the cancer type and 
dictionary of views, i.e. different data collected for this cell line.
"""
struct CellLine <: Sample
	id::String
	cancer_type::String
	views::OrderedDict{Type, DataView}
	CellLine(id::String, cancer_type::String="unknown") = new(id, cancer_type, OrderedDict{Type, DataView}())
end


""" 
A gene. Mostly represented by its HGNC ID in the motivating Dreamchallenge data set.
Preferred ID is Entrez ID, Ensembl-ID can be stored as well.
As an additional data if available, can hold information on being a cancer gene.
"""
# is still mutable here, since we are continuously collecting information, maybe make immutable in the future
mutable struct Gene
	entrez_id::Int64 # <- primary key preferred
	hgnc_id::String
	ensembl_id::String
	cancer_gene::Bool
	Gene() = new(-1,"","",false)
	function Gene(id::String;id_type::String="hgnc_id", is_cancer_gene::Bool=false)
		if id_type == "ensembl_id"
	 		return new(-1, "", id, is_cancer_gene)
	 	end
	 	# hgnc_id is standard type
	 	new(-1, id, "", is_cancer_gene)
	end
	Gene(entrez_id::Int64, hgnc_id::String, ensembl_id::String; is_cancer_gene::Bool=false) = new(entrez_id, hgnc_id, ensembl_id, is_cancer_gene)
	Gene(entrez_id::Int64, is_cancer_gene::Bool=false) = new(entrez_id, "", "", is_cancer_gene)
end

""" 
A protein. If available, paralogs of the protein that are affected by the same antibody are encoded in the hgnc_id.
The field antibody_validated indicates whether cross-reactivity in the assay has been checked
and the antibodies behave linearly in the dilution series.
"""
struct Protein
	hgnc_id::String
	antibody_validated::Bool
end

""" A gene expression."""
struct GeneExpression <: ViewType
	expression_value::Float64
	normalized_value::Float64
	GeneExpression(expression_value::Float64) = new(expression_value, expression_value)
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

struct Methylation <: ViewType
	cgct1::Int64
	cct1::Int64
	illumina_beta_value::Float64
	is_methylated::Bool
	illumina_id::String
	normalized_value::Float64
	function Methylation(beta_value::Float64, cgct1::Int64, cct1::Int64, methylated_threshold::Float64=.2; illumina_id::String="")
		new(cgct1, cct1, beta_value, beta_value≥methylated_threshold, illumina_id, beta_value)
	end
end

"""
Whole transcriptome shotgun sequencing (RNA-seq) data.
Expression value is a fpkm value.
expression_status is 1 if the Ensembl gene model was detected above 
the background noise level.
"""
struct RNASeq <: ViewType
	expression_value::Float64
	normalized_value::Float64
	expression_status::Bool
	RNASeq(expression_value::Float64, expression_status::Bool) = new(expression_value, expression_value, expression_status)
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
struct ExomeSeq <: ViewType
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
	normalized_value::Float64
end

""" Protein abundance value. """
struct RPPA <: ViewType
	protein_abundance::Float64
	normalized_value::Float64
	RPPA(protein_abundance::Float64) = new(protein_abundance, protein_abundance)
end

""" The quantified value for DNA copy number analysis (aggregated across whole genes) """
struct CNV <: ViewType
	gene_level_cnv::Float64
	normalized_value::Float64
	CNV(gene_level_cnv::Float64) = new(gene_level_cnv, gene_level_cnv)
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
mutable struct Drug
	id::String
	affected_genes::Vector{Gene}
	chemical_structure::String
	Drug(id::String) = new(id, Gene[], "")
end

"""
An experimental outcome.
An outcome should be associated with a drug.
The dictionary of outcome values associates the results to cell lines. 
The reason for this is, that we believe the outcome across different
cell lines is correlated according to the molecular structure of the cell line tissue.
In other words, we model the outcome on available cell lines jointly.
"""
mutable struct Outcome
	outcome_values::OrderedDict{CellLine, Float64}
	outcome_type::String
	normalized_outcome_values::OrderedDict{CellLine, Float64}
	outcome_mean::Float64
	outcome_std::Float64
	Outcome(outcome_type::String) = new(OrderedDict{CellLine, Float64}(), outcome_type, OrderedDict{CellLine, Float64}(), .0, 1.)
end

""" 
An experiment comprised of measurements (where each is associated 
with a cell line) and results.
"""
struct Experiment
	results::OrderedDict{Drug, Outcome}
	cell_lines::OrderedDict{String, CellLine}
	genes::OrderedDict{Int64, Gene}
	genes_by_hgnc::OrderedDict{String, Gene}
	genes_by_ensembl::OrderedDict{String, Gene}
	proteins::OrderedDict{String, Protein}
	views::OrderedSet{Type{<:ViewType}}
	function Experiment()
		new(OrderedDict{Drug, Outcome}(), OrderedDict{String, CellLine}(), OrderedDict{Int64, Gene}(), OrderedDict{String, Gene}(), OrderedDict{String, Gene}(), OrderedDict{String, Protein}(), OrderedSet{Type{<:ViewType}}())
	end
end

##############################################################################################################
##############################################################################################################
########################################      PREDICTOR        ###############################################
##############################################################################################################
##############################################################################################################
struct DrugEfficacyPrediction
	experiment::Experiment
	model::PredictionModel
	kernels::OrderedDict{Drug, Vector{Matrix{Float64}}}
	targets::OrderedDict{Drug, Vector{Float64}}
	continuous_kernel::Function
	discrete_kernel::Function
	DrugEfficacyPrediction(experiment::Experiment, model::PredictionModel) = new(experiment, model, Dict{Drug, Vector{Matrix{Float64}}}(), OrderedDict{Drug, Vector{Float64}}())
end
