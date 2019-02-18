##############################################################################################################
##############################################################################################################
########################################         MODEL         ###############################################
##############################################################################################################
##############################################################################################################
# abstract type PredictionModel end
abstract type Model end

abstract type PredictionModelParameters end

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

mutable struct PredictionModel{M<:Model}
	# model_type::Type{<:Model}
	data::Experiment
	precomputations::Dict{String, Any}

	PredictionModel(::Type{M}, data::Experiment) where {M <: Model} = new{M}(data, Dict{String, Any}())
end

struct ModelConfiguration
	parameters::Dict{String, Any}
	function ModelConfiguration()
		new(Dict{String, Any}())
	end
end
