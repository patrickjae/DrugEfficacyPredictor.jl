"""
Initialization of the probability model.
"""
function initialize_probability_model(experiment::DrugEfficacyPredictor.Experiment)
    # set counting variables
    T = length(experiment.results)
    K = length(experiment.views)+1
    N = Vector{Int64}(T)
    # for each drug, store mean and standard deviation, used for normalization
    for (t, drug) in enumerate(keys(experiment.results))
        N[t] = length(experiment.results[drug].outcome_values)
        vals = collect(values(experiment.results[drug].outcome_values))
        experiment.results[drug].outcome_mean = mean(vals)
        experiment.results[drug].outcome_std = stdm(vals, experiment.results[drug].outcome_mean)
    end
    # create predictor model
    DrugEfficacyPredictor.PredictionModel(T, K, N)
end

"""
Creates a DrugEfficacyPredictor object, holding the probabilistic model, all the data and
the kernel functions.
"""
function create_drug_efficacy_predictor(experiment::DrugEfficacyPredictor.Experiment, ::Dict{String, Any}=Dict{String, Any}())
    # make sure that measurements are normalized across cell lines
    if !experiment.is_normalized
        tic()
        normalize_data_views(experiment)
        info("normalizing training data took $(toq()) seconds")
    end

    # normalize outcomes as well
    for (t, drug) in enumerate(keys(experiment.results))
        # N[t] = length(experiment.results[drug].outcome_values)
        vals = collect(values(experiment.results[drug].outcome_values))
        experiment.results[drug].outcome_mean = mean(vals)
        experiment.results[drug].outcome_std = stdm(vals, experiment.results[drug].outcome_mean)
    end


    all_drugs = collect(keys(experiment.results)) # the tasks

    # compute kernels between cell lines once
    base_kernels = Dict{Type{<:DrugEfficacyPredictor.ViewType}, Vector{Matrix{Float64}}}()
    cell_lines = collect(values(experiment.cell_lines))
    # find cell lines that are present in all views
    cell_lines_in_views = filter(cl -> length(cl.views) == length(experiment.views), cell_lines)

    info("all cell lines: $(length(cell_lines_in_views))")
    for v in experiment.views
        # get views of this type for all cell lines
        # cell_lines_with_view = filter(cl -> haskey(cl.views, v), cell_lines)
        data_views = map(cl -> cl.views[v], cell_lines_in_views)
        #workaround for having both expression level and status for RNASeq in one data structure
        kernels = compute_kernel(data_views)
        info("view $v has $(length(cell_lines_in_views)) cell lines if available")
        base_kernels[v] = kernels
        info("computed $(length(kernels)) kernel(s) for view $v")
        # push!(base_kernels, compute_kernel(data_views)...)
    end

    T = length(experiment.results)
    N = Vector{Int64}(T)

    kernels = DataStructures.OrderedDict{DrugEfficacyPredictor.Drug, Vector{Matrix{Float64}}}()
    targets = DataStructures.OrderedDict{DrugEfficacyPredictor.Drug, Vector{Float64}}()
    # for all drugs normalize target data across available cell lines
    for (t, drug) in enumerate(all_drugs)
        info("processing drug $drug")
        kernels[drug] = Vector{Matrix{Float64}}()
        # find cell lines for which we have outcomes for the current drug
        idx_in_outcome = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines_in_views)
        idx_in_cell_lines = findin(cell_lines_in_views, collect(keys(experiment.results[drug].outcome_values)))
        N[t] = length(idx_in_cell_lines)
        tic()
        for v in experiment.views
            # cell_lines_with_view = filter(cl -> haskey(cl.views, v), cell_lines)
            # info("cell lines that have view $v: $(length(cell_lines_with_view))")
            # idx_in_outcome2 = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines)
            # info("view $v: $(length(idx_in_outcome)), overall: $(length(idx2)), N[t]: $(model.N[t])")
            for k in base_kernels[v]
                push!(kernels[drug], k[idx_in_cell_lines, idx_in_cell_lines])
            end
        end
        info("found $(length(kernels[drug])) data views")
        info("computing kernels took $(toq()) sec onds")
        #normalize the outcome data data for each drug across the cell lines
        info("normalizing response data and computing kernel matrices for drug $drug")
        targets[drug] = ((collect(values(experiment.results[drug].outcome_values))[idx_in_outcome] .- experiment.results[drug].outcome_mean)./experiment.results[drug].outcome_std)
        # compute the kernels for all views in the current subset of cell lines
        # kernels[drug] = create_kernel_matrices( cell_lines)
    end

    K = length(kernels[all_drugs[1]])
    pm = DrugEfficacyPredictor.PredictionModel(T, K, N)
    DrugEfficacyPredictor.DrugEfficacyPrediction(experiment, pm)
end

# """API function"""
# create_drug_efficacy_predictor(experiment::DrugEfficacyPredictor.Experiment, data::Dict{String, Any}) = create_drug_efficacy_predictor(experiment)

""" Creates the kernel matrices (one for each "view") that encode cell line similarity. """
function create_kernel_matrices(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, cell_lines::Vector{DrugEfficacyPredictor.CellLine})
    kernels = Vector{Matrix{Float64}}()
    info("found $(length(cell_lines)) cell lines in experiment with results for current drug")
    for v in dep.experiment.views
        # get views of this type for all cell lines
        info("extracting view $v from cell lines if available")
        data_views = map(cl -> cl.views[v], 
            filter(cl -> haskey(cl.views, v), cell_lines))
        #workaround for having both expression level and status for RNASeq in one data structure
        push!(kernels, compute_kernel(data_views)...)
    end
    # TODO
    # add kernel combinations between views and pooling views via gene lists
    kernels
end

function compute_kernel(dvs::Vector{DrugEfficacyPredictor.DataView{K, V}}) where {K <: DrugEfficacyPredictor.KeyType, V <: DrugEfficacyPredictor.RNASeq}
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

function compute_kernel(dvs::Vector{DrugEfficacyPredictor.DataView{K, V}}) where {K <: DrugEfficacyPredictor.KeyType, V <: DrugEfficacyPredictor.ViewType}
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

function prepare_kernel(a::DrugEfficacyPredictor.DataView{<:DrugEfficacyPredictor.KeyType, V}, b::DrugEfficacyPredictor.DataView{<:DrugEfficacyPredictor.KeyType, V}) where {V <:DrugEfficacyPredictor.ViewType}
    common_keys = unique(intersect(a.used_keys, b.used_keys))
    x = map(k -> a.measurements[k], common_keys)
    x_prime = map(k -> b.measurements[k], common_keys)
    (x, x_prime)
end

function prepare_kernel(a::DrugEfficacyPredictor.DataView{<:DrugEfficacyPredictor.KeyType, DrugEfficacyPredictor.ExomeSeq}, b::DrugEfficacyPredictor.DataView{<:DrugEfficacyPredictor.KeyType, DrugEfficacyPredictor.ExomeSeq})
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
                push!(x, DrugEfficacyPredictor.get_measurement_value(filter(m -> m.protein_change == pc, measurements_a))...)
                push!(x_prime, DrugEfficacyPredictor.get_measurement_value(filter(m -> m.protein_change == pc, measurements_a))...)
                if pc != "None" factor_exact_matches += 1 end
            end
        else
            push!(x, mean(DrugEfficacyPredictor.get_measurement_value(a, k)))
            push!(x_prime, mean(DrugEfficacyPredictor.get_measurement_value(b, k)))
        end

    end
    (x, x_prime, Float64(length(x) * factor_exact_matches))
end

# TODO: implement kernel functions for different view types: GeneExpression, 
# Methylation, RNASeq, ExomeSeq, RPPA, CNV
# we employ a Gaussian kernel on all continuous values and a Jaccard kernel for binary-valued data.

function kernel_function(x::Vector{<:DrugEfficacyPredictor.ViewType}, x_prime::Vector{<:DrugEfficacyPredictor.ViewType})

    kernel_function(map(item -> item.normalized_value, x), map(item -> item.normalized_value, x_prime))
end

function kernel_function(x::Vector{DrugEfficacyPredictor.ExomeSeq}, x_prime::Vector{DrugEfficacyPredictor.ExomeSeq})
    # find overlap of genes
    genes_union = unique(union(map(item -> item.gene, x.used_keys), map(item -> item.gene, x_prime.used_keys)))
end

function kernel_function(x::Vector{DrugEfficacyPredictor.RNASeq}, x_prime::Vector{DrugEfficacyPredictor.RNASeq})
    sim_cont = kernel_function(map(item -> item.normalized_value, x), map(item -> item.normalized_value, x_prime))
    sim_dist = kernel_function(map(item -> item.expression_status, x), map(item -> item.expression_status, x_prime))
    (sim_cont, sim_dist)
end

#actual kernel functions
kernel_function(x::Vector{Float64}, x_prime::Vector{Float64}, l::Float64 = Float64(length(x))) = exp(-dot(x-x_prime, x-x_prime)/l)
kernel_function(x::Vector{Bool}, x_prime::Vector{Bool}) = x'*x_prime/(x'*x + x_prime'*x_prime - x'*x_prime)
