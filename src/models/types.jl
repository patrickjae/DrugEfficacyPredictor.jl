##############################################################################################################
##############################################################################################################
########################################         MODEL         ###############################################
##############################################################################################################
##############################################################################################################
# abstract type PredictionModel end
abstract type Model end
abstract type BMTMKLModel <: Model end

abstract type PredictionModelParameters end

mutable struct InferenceConfiguration
	model_type::Type{<:Model}
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

	InferenceConfiguration() = new(BMTMKLModel, 1e-3, 5, 200, "results/", 1, false, false, false, true, true)
end

mutable struct PredictionModel
	model_type::Type{<:Model}
	data::Experiment
	precomputations::Dict{String, Any}

	function PredictionModel(data::Experiment, ic::InferenceConfiguration)
		m = new()
		m.model_type = ic.model_type
		m.data = data

		# make sure that measurements are normalized across cell lines
		if !m.data.is_normalized
			normalize_data_views(m.data)
		end

		#do variance filtering if required
		if ic.do_variance_filtering
			filter_data_views(m.data)
		else
			# make sure we are using all keys
			unfilter_data_views(m.data)
		end

		m.precomputations = Dict{String, Any}()
		return init(m, ic)
	end
end

struct ModelConfiguration
	parameters::Dict{String, Any}
	ModelConfiguration() = new(Dict{String, Any}())
end
