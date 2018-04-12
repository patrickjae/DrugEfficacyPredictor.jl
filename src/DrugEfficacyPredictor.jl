module DrugEfficacyPredictor
using CSV, DataFrames, DataStructures, Distributions, SpecialFunctions
include("distributional_parameters.jl")
include("types.jl")
include("data_handling.jl")
include("dream_challenge_data.jl")


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
	all_drugs = collect(keys(dep.experiment.results)) # the tasks
	# for all drugs normalize target data across available cell lines
	for drug in all_drugs
		cell_lines = collect(keys(dep.experiment.results[drug].outcome_values))
		#normalize the data
		dep.targets[drug] = (collect(values(dep.experiment.results[drug].outcome_values)) .- dep.experiment.results[drug].outcome_mean)./dep.experiment.results[drug].outcome_std
		dep.kernels[drug] = create_kernel_matrices(dep, cell_lines)
	end
end

""" Creates the kernal matrices (one for each "view") that encode cell line similarity. """
function create_kernel_matrices(dep::DrugEfficacyPrediction, cell_lines::Vector{CellLine})
	kernels = Vector{Matrix{Float64}}()
	for v in dep.experiment.views
		data_views = map(cl -> cl.views[v], cell_lines)
		push!(kernels, compute_kernel(dep, data_views))
	end
	kernels
end

function compute_kernel(dep::DrugEfficacyPrediction, dvs::Vector{DataView{K, V}}) where {K, V <: ViewType}
	N = length(dvs)
	k = Matrix{Float64}(N,N)

	# find overlapping genes for each pairing
	for i in 1:N, j in i:N
		common_keys = intersect(collect(keys(dvs[i].measurements)), collect(keys(dvs[j].measurements)))
		x = map(k -> dvs[i].measurements[k], common_keys)
		x_prime = map(k -> dvs[j].measurements[k], common_keys)

		# create kernel matrix entries, kernel matrix will be symmetric
		k[i,j] = kernel_function(dep, x, x_prime)
		if i≠j k[j,i] = k[i,j] end
	end
	k
end

# TODO: implement kernel functions for different view types: GeneExpression, 
# Methylation, RNASeq, ExomeSeq, RPPA, CNV

function kernel_function(dep::DrugEfficacyPrediction, x::T, x_prime::T) where {T <: ViewType}

end


# actual variational inference algorithm
function parameter_inference(dep::DrugEfficacyPrediction)
	m = dep.model
	all_tasks = collect(keys(dep.experiment.results))

	kernel_products = Dict{Drug, Matrix{Float64}}()
	for (t, d) in enumerate(all_tasks)
		views = dep.kernels[d]
		kp = zeros(m.N[t], m.N[t])
		[kp += kernel'*kernel for kernel in views]
		kernel_products[d] = kp
	end
	# updates for model parameters in turn
	ll = 0
	# lambda, a
	for (t, drug) in enumerate(all_tasks)
		aaT = expected_squared_value(m.a[t])
		for n in 1:m.N[t]
			m.λ[t][n].variational_a = m.λ[t][n].prior_a + .5
			m.λ[t][n].variational_b = m.λ[t][n].prior_b + .5*aaT[n,n]
		end
		m.a[t].variational_covariance = inv(diagm(expected_value.(m.λ[t])) + expected_value(m.ν[t]).*kernel_products[drug])
		
		m.a[t].variational_mean = m.a[t].variational_covariance*
	end
	# gamma
	for t in 1:m.T
		m.ɣ.variational_a = m.ɣ.prior_a + .5
		m.ɣ.variational_b = m.ɣ.prior_b + expected_squared_value(m.b[t])
		#TODO: compute log likelihood E_q[ln p(gamma)] - E_q[ln q(gamma)]
	end


end


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
