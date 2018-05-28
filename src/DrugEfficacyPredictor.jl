module DrugEfficacyPredictor
using CSV, DataFrames, DataStructures, Distributions, SpecialFunctions
include("distributional_parameters.jl")
include("types.jl")
include("data_handling.jl")
include("dream_challenge_data.jl")
include("webservice.jl")

include("import_json.jl")


"""
Adjust info and error functions to include time stamp.
"""
info(message::String) = Base.info("\t[$(now())] $message")
error(message::String) = Base.error("\t[$(now())] $message")
warn(message::String) = Base.warn("[$(now())] $message")

"""
Initialization of the probability model.
"""
function initialize_probability_model(experiment::Experiment)
	# set counting variables
	T = length(experiment.results)
	K = length(experiment.views)
	N = Vector{Int64}(T)
	# for each drug, store mean and standard deviation, used for normalization
	for (t, drug) in enumerate(keys(experiment.results))
		N[t] = length(experiment.results[drug].outcome_values)
		vals = collect(values(experiment.results[drug].outcome_values))
		experiment.results[drug].outcome_mean = mean(vals)
		experiment.results[drug].outcome_std = stdm(vals, experiment.results[drug].outcome_mean)
	end
	# create predictor model
	PredictionModel(T, K, N)
end

"""
Creates a DrugEfficacyPredictor object, holding the probabilistic model, all the data and
the kernel functions.
"""
function create_drug_efficacy_predictor(experiment::Experiment)
	# initialize object
	dep = DrugEfficacyPrediction(experiment, initialize_probability_model(experiment))
	# make sure that measurements are normalized across cell lines
	if !experiment.is_normalized
		tic()
		normalize_data_views(experiment)
		info("normalizing training data took $(toq()) seconds")
	end
	all_drugs = collect(keys(dep.experiment.results)) # the tasks
	# for all drugs normalize target data across available cell lines
	info("normalizing response data and computing kernel matrices per drug")
	for drug in all_drugs
		info("processing drug $drug")
		# find cell lines for which we have outcomes for the current drug
		cell_lines = collect(keys(dep.experiment.results[drug].outcome_values))
		#normalize the outcome data data for each drug across the cell lines
		dep.targets[drug] = (collect(values(dep.experiment.results[drug].outcome_values)) .- dep.experiment.results[drug].outcome_mean)./dep.experiment.results[drug].outcome_std
		# compute the kernels for all views in the current subset of cell lines
		tic()
		dep.kernels[drug] = create_kernel_matrices(dep, cell_lines)
		info("computing kernels took $(toq()) seconds")
	end
	dep
end

"""API function"""
create_drug_efficacy_predictor(experiment::DrugEfficacyPredictor.Experiment, data::Dict{String, Any}) = create_drug_efficacy_predictor(experiment)

""" Creates the kernel matrices (one for each "view") that encode cell line similarity. """
function create_kernel_matrices(dep::DrugEfficacyPrediction, cell_lines::Vector{CellLine})
	kernels = Vector{Matrix{Float64}}()
	info("found $(length(cell_lines)) cell lines in experiment with results for current drug")
	for v in dep.experiment.views
		# get views of this type for all cell lines
		info("extracting view $v from cell lines if available")
		data_views = map(cl -> cl.views[v], 
			filter(cl -> haskey(cl.views, v), cell_lines))
		#workaround for having both expression level and status for RNASeq in one data structure
		push!(kernels, compute_kernel(dep, data_views)...)
	end
	kernels
end

function compute_kernel(dep::DrugEfficacyPrediction, dvs::Vector{DataView{K, V}}) where {K <: KeyType, V <: RNASeq}
	N = length(dvs)
	k_cont = Matrix{Float64}(N,N)
	k_dist = Matrix{Float64}(N,N)
	key_compare_time = 0
	kernel_compute_time = 0

	for i in 1:N, j in i:N
		tic()
		(x, x_prime) = prepare_kernel(dvs[i], dvs[j])
		key_compare_time += toq()
		tic()
		(sim_cont, sim_dist) = kernel_function(x, x_prime)
		kernel_compute_time = toq()

		k_cont[i,j] = sim_cont
		k_dist[i,j] = sim_dist
		if i≠j
			k_cont[j,i] = k_cont[i,j]
			k_dist[j,i] = k_dist[i,j]
		end
	end
	info("comparison of keys took $key_compare_time seconds, kernel computations took $kernel_compute_time seconds")
	[k_cont, k_dist]
end

function compute_kernel(dep::DrugEfficacyPrediction, dvs::Vector{DataView{K, V}}) where {K <: KeyType, V <: ViewType}
	N = length(dvs)
	k = Matrix{Float64}(N,N)
	key_compare_time = 0
	kernel_compute_time = 0

	# find overlapping genes for each pairing
	for i in 1:N, j in i:N
		tic()
		kernel_params = prepare_kernel(dvs[i], dvs[j])
		key_compare_time += toq()
		# create kernel matrix entries, kernel matrix will be symmetric
		tic()
		k[i,j] = kernel_function(kernel_params...)
		kernel_compute_time = toq()
		if i≠j k[j,i] = k[i,j] end
	end
	info("comparison of keys took $key_compare_time seconds, kernel computations took $kernel_compute_time seconds")
	[k]
end

function prepare_kernel(a::DataView{<:KeyType, V}, b::DataView{<:KeyType, V}) where {V <:ViewType}
	common_keys = unique(intersect(a.used_keys, b.used_keys))
	x = map(k -> a.measurements[k], common_keys)
	x_prime = map(k -> b.measurements[k], common_keys)
	(x, x_prime)
end

function prepare_kernel(a::DataView{<:KeyType, ExomeSeq}, b::DataView{<:KeyType, ExomeSeq})
	common_keys = unique(intersect(a.used_keys, b.used_keys))
	x = Vector{Float64}()
	x_prime = Vector{Float64}()
	factor_exact_matches = 1.
	for k in common_keys
		# check for match in protein change (this happens rarely and we want to emphasize this)
		# if there is match, make sure to include only the ExomeSeq measurement containing it
		measurements_a = a.measurements[k]
		measurements_b = b.measurements[k]
		protein_change_intersection = intersect(map(es -> es.protein_change, measurements_a), map(es -> es.protein_change, measurements_b))
		if length(protein_change_intersection) > 0
			for pc in protein_change_intersection
				push!(x, get_measurement_value(filter(m -> m.protein_change == pc, measurements_a))...)
				push!(x_prime, get_measurement_value(filter(m -> m.protein_change == pc, measurements_a))...)
				if pc != "None" factor_exact_matches += 1 end
			end
		else
			push!(x, mean(get_measurement_value(a, k)))
			push!(x_prime, mean(get_measurement_value(b, k)))
		end

	end
	(x, x_prime, Float64(length(x) * factor_exact_matches))
end

# TODO: implement kernel functions for different view types: GeneExpression, 
# Methylation, RNASeq, ExomeSeq, RPPA, CNV
# we employ a Gaussian kernel on all continuous values and a Jaccard kernel for binary-valued data.

function kernel_function(x::Vector{<:ViewType}, x_prime::Vector{<:ViewType})

	kernel_function(map(item -> item.normalized_value, x), map(item -> item.normalized_value, x_prime))
end

function kernel_function(x::Vector{ExomeSeq}, x_prime::Vector{ExomeSeq})
	# find overlap of genes
	genes_union = unique(union(map(item -> item.gene, x.used_keys), map(item -> item.gene, x_prime.used_keys)))
end

function kernel_function(x::Vector{RNASeq}, x_prime::Vector{RNASeq})
	sim_cont = kernel_function(map(item -> item.normalized_value, x), map(item -> item.normalized_value, x_prime))
	sim_dist = kernel_function(map(item -> item.expression_status, x), map(item -> item.expression_status, x_prime))
	(sim_cont, sim_dist)
end

#actual kernel functions
kernel_function(x::Vector{Float64}, x_prime::Vector{Float64}, l::Float64 = Float64(length(x))) = exp(-dot(x-x_prime, x-x_prime)/l)
kernel_function(x::Vector{Bool}, x_prime::Vector{Bool}) = x'*x_prime/(x'*x + x_prime'*x_prime - x'*x_prime)

include("inference.jl")

function prepare_predictor(experiment::Experiment)
	# compute the kernel matrices
end




########################
######### TODO #########
########################

# steps

# populate model

# expose API

# cross-validation during training

# split application and API


end # module
