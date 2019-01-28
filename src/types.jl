##############################################################################################################
##############################################################################################################
########################################         DATA          ###############################################
##############################################################################################################
##############################################################################################################
export Gene, Protein, GeneExpression, Methylation,
	RNASeq, RNASeqCall, ExomeSeq, RPPA, CNV,
	Drug, Pathway, DataView, CellLine, Outcome,
	Experiment, DrugEfficacyPrediction,
	ModelConfiguration, InferenceConfiguration,
	PredictionModel

""" An abstract type for samples, can be cell line or human tissue (latter not yet implemented) """
abstract type Sample end

""" Generic view type, i.e. a way to describe a cell line """
abstract type ViewType end

abstract type NAViewType <: ViewType end

abstract type ValueViewType <: ViewType end

abstract type SingleViewType <: ValueViewType end

abstract type VectorViewType <: ValueViewType end

"""
A gene. Mostly represented by its HGNC ID in the motivating Dreamchallenge data set.
Preferred ID is Entrez ID, Ensembl-ID can be stored as well.
As an additional data if available, can hold information on being a cancer gene.
"""
mutable struct Gene
	# is still mutable here, since we are continuously collecting information, maybe make immutable in the future
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

import Base.isless
function isless(x::Gene, y::Gene)
	if x.entrez_id != -1 || y.entrez_id != -1
		return isless(x.entrez_id, y.entrez_id)
	end
	if x.ensembl_id != "" || y.ensembl_id != ""
		return isless(x.ensembl_id, y.ensembl_id)
	end
	if x.hgnc_id != "" || y.hgnc_id != ""
		return isless(x.hgnc_id, y.hgnc_id)
	end
	# if isless(x.entrez_id, y.entrez_id) return true
	# elseif x.entrez_id == y.entrez_id && (isless(x.ensembl_id, y.ensembl_id) || isless(x.hgnc_id, y.hgnc_id)) return true
	# else return false end
end
isless(x::Protein, y::Protein) = isless(x.hgnc_id, y.hgnc_id)

""" The types allowed as keys in data views """
const KeyType = Union{Gene, Protein}

""" A gene expression."""
mutable struct GeneExpression <: SingleViewType
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
mutable struct Methylation <: SingleViewType
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
mutable struct RNASeq <: SingleViewType
	expression_value::Float64
	normalized_value::Float64
	RNASeq(expression_value::Float64) = new(expression_value, expression_value)
end

"""
Whole transcriptome shotgun sequencing (RNA-seq) data.
expression_status is 1 if the Ensembl gene model was detected above
the background noise level.
"""
mutable struct RNASeqCall <: SingleViewType
	expression_status::Bool
	normalized_value::Bool
	RNASeqCall(expression_status::Bool) = new(expression_status, expression_status)
end

"""
Whole exome sequencing data. This is a vector view type because we may have multiple measurements for a single Gene key.

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
mutable struct ExomeSeq <: VectorViewType
	num_cosmic::Int64
	variant_effect::String
	protein_change::String
	nucleotid_change::String
	variant_confidence::Float64
	norm_zygosity::String
	norm_reference_count::Float64
	norm_variant_count::Float64
	tumor_zygosity::String
	tumor_reference_count::Float64
	tumor_variant_count::Float64
	reference_mismatch_avg::Float64
	variant_mismatch_avg::Float64
	reference_mismatch_sum::Float64
	variant_mismatch_sum::Float64
	reference_dist3effective_avg::Float64
	variant_dist3effective_avg::Float64
	details::String
	normalized_value::Float64
	function ExomeSeq(protein_change::AbstractString; reference_mismatch_sum::Float64=0., reference_mismatch_avg::Float64=0.,
		reference_dist3effective_avg::Float64=0., variant_mismatch_sum::Float64=0., variant_mismatch_avg::Float64=0.,
		variant_dist3effective_avg::Float64=0., num_cosmic::Int64=0, variant_effect::AbstractString="", nucleotid_change::AbstractString="",
		variant_confidence::Float64=1., norm_zygosity::AbstractString="", norm_reference_count::Float64=0., norm_variant_count::Float64=0.,
		tumor_zygosity::AbstractString="", tumor_reference_count::Float64=0., tumor_variant_count::Float64=0., details::AbstractString="")
		data_summary = reference_mismatch_avg + variant_mismatch_avg
						+ reference_mismatch_sum + variant_mismatch_sum
						+ reference_dist3effective_avg + variant_dist3effective_avg
		new(num_cosmic, variant_effect, protein_change, nucleotid_change, variant_confidence,
			norm_zygosity, norm_reference_count, norm_variant_count, tumor_zygosity, tumor_reference_count, tumor_variant_count,
			reference_mismatch_avg, variant_mismatch_avg, reference_mismatch_sum, variant_mismatch_sum,
			reference_dist3effective_avg, variant_dist3effective_avg, details,
			data_summary)
	end
end


""" Protein abundance value. """
mutable struct RPPA <: SingleViewType
	protein_abundance::Float64
	normalized_value::Float64
	RPPA(protein_abundance::Float64) = new(protein_abundance, protein_abundance)
end

""" The quantified value for DNA copy number analysis (aggregated across whole genes) """
mutable struct CNV <: SingleViewType
	gene_level_cnv::Float64
	normalized_value::Float64
	CNV(gene_level_cnv::Float64) = new(gene_level_cnv, gene_level_cnv)
end

"""
A drug, identified by a name.
Additional information may include knowledge about affected genes to
refine the number of features for predictions, this is useful for incorporating pathway information.
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
	Drug(id::String; affected_genes::Vector{Gene}=Vector{Gene}(), chemical_structure::String="") = new(id, affected_genes, chemical_structure)
end

"""
Pathway information.
Pathway information is subsumed into a set of genes that affect the signaling pathway.
"""
struct Pathway{K<:KeyType}
	id::String
	name::String
	genes::Vector{K}
	Pathway{K}(genes::Vector{K}; name::String="", id::String="") where K = new(id, name, genes)
end


"""
A data view, typed with a identifier and data type.
The identifier is mostly of type Gene, for RPPA data it is Protein.
In general, similarity matrices for samples are computed on the basis of
these identifiers. In other words, for each data view, a cell line is represented
as a vector of genes (or proteins) with associated values. The identifiers determine
the ordering of this vector.
"""
struct DataView{K <: KeyType, V <: ViewType}
	cell_line_id::String
	key_type::Type{K}
	view_type::Type{V}
	used_keys::SortedSet{K}
	common_keys::Dict{String, Vector{K}}
	measurements::Union{Dict{K, V}, MultiDict{K, V}}
	DataView{K,V}(cell_line_id::String) where {K,V <: SingleViewType} = new(cell_line_id, K, V, SortedSet{K}(), Dict{String, Vector{K}}(), Dict{K,V}())
	DataView{K,V}(cell_line_id::String) where {K,V <: VectorViewType} = new(cell_line_id, K, V, SortedSet{K}(), Dict{String, Vector{K}}(), MultiDict{K,V}())
	DataView{K,V}(cell_line_id::String) where {K, V <: NAViewType} = new(cell_line_id, K, V, SortedSet{K}(), Dict{String, Vector{K}}(), Dict{K,V}())
end


"""
A cell line comprised with its ID, information on the cancer type and
dictionary of views, i.e. different data collected for this cell line.
"""
mutable struct CellLine <: Sample
	id::String
	cancer_type::String
	views::Dict{Type{<:ViewType}, DataView{<:KeyType, <:ViewType}}
	in_test_set::Bool
	CellLine(id::String, cancer_type::String="unknown"; in_test_set::Bool=false) = new(id, cancer_type, Dict{Type{<:ViewType}, DataView{<:KeyType, <:ViewType}}(), in_test_set)
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
	drug_id::String
	outcome_values::OrderedDict{CellLine, Float64}
	outcome_type::String
	normalized_outcome_values::OrderedDict{CellLine, Float64}
	outcome_mean::Float64
	outcome_std::Float64
	Outcome(drug_id::String, outcome_type::String) = new(drug_id, OrderedDict{CellLine, Float64}(), outcome_type, OrderedDict{CellLine, Float64}(), .0, 1.)
end



# mutable struct Prediction
# 	predicted_outcome::Dict{Drug, Float64}
# 	predicted_rank::Dict{Drug, Int64}
# 	Prediction() = new(Dict{Drug, Float64}(), Dict{Drug,Int64}())
# end

"""
An experiment comprised of measurements (where each is associated
with a cell line) and results. This data structure holds all data that is provided
by the user, i.e. including pathway information, different data views, training and test responses etc.
"""
mutable struct Experiment
	results::OrderedDict{Drug, Outcome}
	# test_results::OrderedDict{Drug, Outcome}
	cell_lines::OrderedDict{String, CellLine}
	drugs::Dict{String, Drug}
	genes::Dict{Int64, Gene}
	genes_by_hgnc::Dict{String, Gene}
	genes_by_ensembl::Dict{String, Gene}
	proteins::Dict{String, Protein}
	views::OrderedSet{Type{<:ViewType}}
	statistics::Dict{Type{<:ViewType}, OrderedDict{KeyType, Tuple{Float64, Float64}}}
	pathway_information::OrderedDict{String, Pathway}
	is_normalized::Bool
	internal_id::String
	function Experiment(id::String)
		new(
			OrderedDict{Drug, Outcome}(), # results
			# OrderedDict{Drug, Outcome}(), # test results
			OrderedDict{String, CellLine}(), # cell_lines
			Dict{String, Drug}(), # drugs
			Dict{Int64, Gene}(), # genes
			Dict{String, Gene}(), # genes by hgnc
			Dict{String, Gene}(), # genes by ensembl
			Dict{String, Protein}(), # proteins
			OrderedSet{Type{<:ViewType}}(), # views
			Dict{Type{<:ViewType}, OrderedDict{KeyType, Tuple{Float64, Float64}}}(), # statistics
			OrderedDict{String, Pathway}(), #pathways
			false, # is normalized
			id
			)
	end
end

##############################################################################################################
##############################################################################################################
########################################      PREDICTOR        ###############################################
##############################################################################################################
##############################################################################################################
mutable struct DrugEfficacyPrediction
	experiment::Experiment

	base_kernels::OrderedDict{Type{<:ViewType}, Matrix{Float64}}
	pathway_specific_kernels::OrderedDict{Type{<:ViewType}, Vector{Matrix{Float64}}}

	# training data
	kernels::OrderedDict{Drug, Vector{Matrix{Float64}}}
	targets::OrderedDict{Drug, Vector{Float64}}

	# validation data
	validation_kernels::OrderedDict{Drug, Vector{Matrix{Float64}}}
	validation_targets::OrderedDict{Drug, Vector{Float64}}

	# test data
	cross_kernels::OrderedDict{Drug, Vector{Matrix{Float64}}}
	test_targets::OrderedDict{Drug, Vector{Float64}}

	T::Int64
	K::Int64
	N::Vector{Int64}
	subsume_pathways::Bool
	function DrugEfficacyPrediction(experiment::Experiment, T::Int64, K::Int64, N::Vector{Int64})
		new(
			experiment,
			OrderedDict{Type{<:ViewType}, Matrix{Float64}}(),
			OrderedDict{Type{<:ViewType}, Vector{Matrix{Float64}}}(),
			OrderedDict{Drug, Vector{Matrix{Float64}}}(),
			OrderedDict{Drug, Vector{Float64}}(),
			OrderedDict{Drug, Vector{Matrix{Float64}}}(),
			OrderedDict{Drug, Vector{Float64}}(),
			OrderedDict{Drug, Vector{Matrix{Float64}}}(),
			OrderedDict{Drug, Vector{Float64}}(),
			T, K, N, true
		)
	end
end

##############################################################################################################
##############################################################################################################
########################################         MODEL         ###############################################
##############################################################################################################
##############################################################################################################
mutable struct ModelConfiguration
	# Gamma params for precision on bias term b
	⍺_ɣ::Float64
	β_ɣ::Float64
	# Gamma params for precision on weights a
	⍺_λ::Float64
	β_λ::Float64
	# Gamma params for precision on outcome y
	⍺_ε::Float64
	β_ε::Float64
	# Gamma params for precision on intermediate results g
	⍺_ν::Float64
	β_ν::Float64
	# Gamma params for precision on kernel weights e
	⍺_⍵::Float64
	β_⍵::Float64
	# normal params for mean on bias
	μ_b::Float64
	𝜎_b::Float64
	# normal params for mean on kernel weights
	μ_e::Float64
	𝜎_e::Float64
	# mv normal params for weights
	μ_a::Float64
	Σ_a::Float64
	# mv normal weights for intermediate results
	μ_g::Float64
	Σ_g::Float64

	ModelConfiguration(⍺::Float64, β::Float64, μ::Float64, 𝜎::Float64) = new(⍺, β, ⍺, β, ⍺, β, ⍺, β, ⍺, β, 0., 𝜎, μ, 𝜎, μ, 𝜎, μ, 𝜎)
	ModelConfiguration() = new(1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 0., 2., 1., 2., 1., 2., 1., 2.)
end

mutable struct InferenceConfiguration
	convergence_criterion::Float64
	min_iter::Int64
	max_iter::Int64
	target_dir::String
	fold_num::Int64
	do_gridsearch::Bool
	do_cross_validation::Bool
	compute_wpc_index::Bool
	subsume_pathways::Bool
	do_variance_filtering::Bool

	InferenceConfiguration() = new(1e-3, 5, 200, "results/", 1, false, false, false, true, true)
end
mutable struct PredictionModel
	# precision parameters, T-dimensional, one per drug/task
	# ɣ::IsometricPrecisionParameter #precision to b <-- could be a matrix introducing covariance between drugs
	ɣ::Vector{GammaParameter}
	# drug specific scalars, T-dimensional
	b::Vector{NormalParameter}

	# λ::Vector{IsometricPrecisionParameter} #precision to a <-- introduce covariance between cell lines
	λ::Vector{Vector{GammaParameter}}
	# drug specific vectors, TxN_t dimensional
	# a::Vector{Vector{NormalParameter}}
	a::Vector{MvNormalParameter}

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

	function PredictionModel(T::Int64, K::Int64, N::Vector{Int64}; model_config::ModelConfiguration = ModelConfiguration())
		p = new()
		p.T = T
		p.K = K
		p.N = N

		p.ɣ = VectorGammaParameter(undef, T)
		p.b = Vector{NormalParameter}(undef, T)
		for t in 1:T
			p.ɣ[t] = GammaParameter(model_config.⍺_ɣ, model_config.β_ɣ)
			p.b[t] = NormalParameter(model_config.μ_b, model_config.𝜎_b)
			set_variable_name(p.ɣ[t], "ɣ[$t]")
			set_variable_name(p.b[t], "b[$t]")
		end

		p.λ = Vector{VectorGammaParameter}(undef, T)
		# p.a = Vector{Vector{NormalParameter}}(T)
		p.a = Vector{MvNormalParameter}(undef, T)
		for t in 1:T
			p.λ[t] = VectorGammaParameter(undef, N[t])
			p.a[t] = MvNormalParameter(model_config.μ_a, model_config.Σ_a, N[t])
			# p.a[t] = Vector{NormalParameter}(N[t])
			set_variable_name(p.a[t], "a[$t]")
			for n in 1:N[t]
				p.λ[t][n] = GammaParameter(model_config.⍺_λ, model_config.β_λ)
				set_variable_name(p.ɣ[t], "λ[$t][$n]")
				# p.a[t][n] = NormalParameter(m_0, s_0)
			end
		end

		p.ε = VectorGammaParameter(undef, T)
		for t in 1:T
			p.ε[t] = GammaParameter(model_config.⍺_ε, model_config.β_ε)
			set_variable_name(p.ε[t], "ε[$t]")
		end

		p.ν = VectorGammaParameter(undef, T)
		for t in 1:T
			p.ν[t] = GammaParameter(model_config.⍺_ν, model_config.β_ν)
			set_variable_name(p.ν[t], "ν[$t")
		end
		p.G = Matrix{MvNormalParameter}(undef, T, K)
		for t in 1:T, k in 1:K
			p.G[t,k] = MvNormalParameter(model_config.μ_g, model_config.Σ_g, N[t])
			set_variable_name(p.G[t,k], "G[$t,$k]")
		end

		p.⍵ = VectorGammaParameter(undef, K)
		p.e = Vector{NormalParameter}(undef, K)
		for k in 1:K
			p.⍵[k] = GammaParameter(model_config.⍺_⍵, model_config.β_⍵)
			p.e[k] = NormalParameter(model_config.μ_e, model_config.𝜎_e)
			set_variable_name(p.⍵[k], "⍵[$k]")
			set_variable_name(p.e[k], "e[$k]")
		end
		p
	end
end
