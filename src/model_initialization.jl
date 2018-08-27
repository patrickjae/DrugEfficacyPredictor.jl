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
function create_drug_efficacy_predictor(experiment::Experiment, ::Dict{String, Any}=Dict{String, Any}())
    # make sure that measurements are normalized across cell lines
    if !experiment.is_normalized
        tic()
        normalize_data_views(experiment)
        info("normalizing training data took $(toq()) seconds")
    end

    # collect outcomes statistics
    for (t, drug) in enumerate(keys(experiment.results))
        vals = collect(values(experiment.results[drug].outcome_values))
        experiment.results[drug].outcome_mean = mean(vals)
        experiment.results[drug].outcome_std = stdm(vals, experiment.results[drug].outcome_mean)
    end


    all_drugs = collect(keys(experiment.results)) # the tasks

    # compute base kernels between cell lines once, one for each ViewType
    base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
    pathway_specific_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    grouped_data_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    cell_lines = collect(values(experiment.cell_lines))
    # find cell lines that are present in all views, i.e. we need all data views for a cell line
    cell_lines_in_views = filter(cl -> length(cl.views) == length(experiment.views), cell_lines)

    (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(experiment, cell_lines)

    # K = 0
    # # compute a base kernel for each view, containing all cell lines that are present in all views
    # for v in experiment.views
    #     info("processing view $v")
    #     # get views of this type for all cell lines
    #     data_views = map_data_views(cell_lines, v)
    #     info("found $(length(data_views)) data records")
    #     # compute the base kernel for this view
    #     base_kernels[v] = compute_kernel(data_views)
    #     K += 1
    #     pathway_specific_kernels[v] = Vector{Matrix{Float64}}()
    #     for pathway in experiment.pathway_information
    #         push!(pathway_specific_kernels[v], compute_kernel(data_views, pathway=pathway))
    #         K += 1
    #     end
    #     # push!(base_kernels, compute_kernel(data_views)...)
    # end

    T = length(experiment.results)
    N = Vector{Int64}(T)

    # compute drug specific kernels
    kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    cross_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    test_targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    for (t, drug) in enumerate(all_drugs)
        kernels[drug] = Vector{Matrix{Float64}}()
        cross_kernels[drug] = Vector{Matrix{Float64}}()

        # find cell lines for which we have outcomes for the current drug
        # indices of the cell lines with outcome that are available for all views (in outcome structure)
        # idx_in_outcome = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines_in_views)
        idx_in_outcome = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines)
        # find cell lines in test set for which we a measurement for the current drug
        # idx_in_test_outcome = findin(collect(keys(experiment.test_results[drug].outcome_values)), cell_lines_in_views)
        idx_in_test_outcome = findin(collect(keys(experiment.test_results[drug].outcome_values)), cell_lines)
        # indices of present cell lines that have outcome data (in cell lines vector)
        # idx_in_cell_lines = findin(cell_lines_in_views, collect(keys(experiment.results[drug].outcome_values)))
        idx_in_cell_lines = findin(cell_lines, collect(keys(experiment.results[drug].outcome_values)))
        # indices of present cell lines that have a test outcome
        # test_idx_in_cell_lines = findin(cell_lines_in_views, collect(keys(experiment.test_results[drug].outcome_values)))
        test_idx_in_cell_lines = findin(cell_lines, collect(keys(experiment.test_results[drug].outcome_values)))
        # only deal with cell lines that are available for all views AND have an outcome for the current drug
        N[t] = length(idx_in_cell_lines)

        for v in experiment.views
            # cell_lines_with_view = filter(cl -> haskey(cl.views, v), cell_lines)
            # info("cell lines that have view $v: $(length(cell_lines_with_view))")
            # idx_in_outcome2 = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines)
            # info("view $v: $(length(idx_in_outcome)), overall: $(length(idx2)), N[t]: $(model.N[t])")

            # add base kernel restricted to cell lines for which we have outcome
            push!(kernels[drug], base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            push!(cross_kernels[drug], base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])

            for pw_kernel in pathway_specific_kernels[v]
                push!(kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                push!(cross_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
            end
        end
        # compute additional kernels
        # grouped kernels
        for v in experiment.views

        end
        #normalize the outcome data data for each drug across the cell lines
        targets[drug] = ((collect(values(experiment.results[drug].outcome_values))[idx_in_outcome] .- experiment.results[drug].outcome_mean)./experiment.results[drug].outcome_std)
        # test outcomes, not normalized
        test_targets[drug] = collect(values(experiment.test_results[drug].outcome_values))[idx_in_test_outcome]
    end
    pm = PredictionModel(T, K, N)
    dep = DrugEfficacyPrediction(experiment, pm)
    for drug in all_drugs
        dep.kernels[drug] = kernels[drug]
        dep.targets[drug] = targets[drug]
        dep.cross_kernels[drug] = cross_kernels[drug]
        dep.test_targets[drug] = test_targets[drug]
    end
    dep
end

function compute_all_kernels(experiment::Experiment, cell_lines::Vector{CellLine}, cell_lines_test::Union{Vector{CellLine}, Void}=nothing)
    base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
    pathway_specific_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    K = 0
    # compute a base kernel for each view, containing all cell lines that are present in all views
    num_views = length(experiment.views)
    println("computing base kernels")
    Threads.@threads for v_id in 1:num_views
        # get views of this type for all cell lines
        tic()
        v = experiment.views[v_id]
        data_views = [map_data_views(cell_lines, v)]
        if cell_lines_test != nothing
            push!(data_views, map_data_views(cell_lines_test, v))
        end
        # compute the base kernel for this view
        base_kernels[v] = compute_kernel(data_views...)
        # info("computing base kernel for $v (took $(toq()) seconds)")
        K += 1
    end
    println("computing pathway specific kernels")
    for v in experiment.views
        # get views of this type for all cell lines
        data_views = [map_data_views(cell_lines, v)]
        if cell_lines_test != nothing
            push!(data_views, map_data_views(cell_lines_test, v))
        end

        num_pws = length(experiment.pathway_information)
        pathway_specific_kernels[v] = Vector{Matrix{Float64}}(num_pws)
        Threads.@threads for p in 1:num_pws
            pathway = experiment.pathway_information[p]
            tic()
            pathway_specific_kernels[v][p] = compute_kernel(data_views..., pathway=pathway)
            # info("computing pathay specific kernel for $v ($(pathway.id)) (took $(toq()) seconds)")
        end
        K += num_pws
        # push!(base_kernels, compute_kernel(data_views)...)
    end
    println("done computing kernels")
    (K, base_kernels, pathway_specific_kernels)
end

# function compute_all_kernels(experiment::Experiment, cell_lines_training::Vector{CellLine}, cell_lines_test::Vector{CellLine})
#     base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
#     pathway_specific_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
#     K = 0
#     for v in experiment.views
#         data_views_train = map_data_views(cell_lines_training, v)
#         data_views_test = map_data_views(cell_lines_test, v)
#         info("computing base kernel for $v")
#         base_kernels[v] = compute_kernel(data_views_train, data_views_test)
#         K += 1
#         pathway_specific_kernels[v] = Vector{Matrix{Float64}}(length(experiment.pathway_information))
#         Threads.@threads for (p, pathway) in enumerate(experiment.pathway_information)
#             info("computing kernel for pathway $(pathway.id)")
#             pathway_specific_kernels[v][p] = compute_kernel(data_views_train, data_views_test, pathway=pathway)
#         end
#         K += length(pathway_specific_kernels[v])
#     end
#     (K, base_kernels, pathway_specific_kernels)
# end


function map_data_views(cell_lines::Vector{CellLine}, view_type::Type{<:ViewType})
    ret = Vector{DataView{<:KeyType, <:ViewType}}()
    # ret = Vector{DataView{Gene, Union{view_type, NAViewType}}}()
    # if view_type == RPPA
    #     ret = Vector{DataView{Protein, Union{view_type, NAViewType}}}()
    # end
    for cl in cell_lines
        if haskey(cl.views, view_type)
            push!(ret, cl.views[view_type])
        else
            if view_type == RPPA
                push!(ret, DataView{Protein, NAViewType}(cl.id))
            else
                push!(ret, DataView{Gene, NAViewType}(cl.id))
            end
        end
    end
    ret
end

# function map_data_views(cell_lines::Vector{CellLine}, view_type::Type{RPPA})
#     map(cl -> if haskey(cl.views, view_type) cl.views[view_type] else DataView{Protein, NAViewType}(cl.id) end, cell_lines)
# end

get_complete_cell_lines(experiment::Experiment) = filter(cl -> length(cl.views) == length(experiment.views), collect(values(experiment.cell_lines)))

function compute_kernel(dvs::Vector{DataView}; pathway::Union{Pathway, Void} = nothing)
    N = length(dvs)
    k = Matrix{Float64}(N,N)
    key_compare_time = 0
    # kernel_compute_time = 0

    # find overlapping genes for each pairing
    tic()
    for i in 1:N, j in i:N
        k[i,j] = compute_kernel_value(dvs[i], dvs[j], pathway)

        if isnan(k[i,j])
            k[i,j] = 1e-7
        end
        if i≠j k[j,i] = k[i,j] end
    end
    # kernel_compute_time = toq()

    # if !issymmetric(k) warn("kernel is not symmetric") end
    if !isposdef(k) 
        h_k = Array(Hermitian(k))
        # warn("kernel is not positive definite, using only upper triangular form (is pos def: $(isposdef(h_k)), sum diff: $(sum(abs.(h_k - k))))") 
        k = h_k
    end
    # info("kernel computations took $kernel_compute_time seconds")
    k
end

function compute_kernel(dvs_train::Vector{DataView}, dvs_test::Vector{DataView}; pathway::Union{Pathway, Void}=nothing)
    N = length(dvs_train)
    M = length(dvs_test)
    k = Matrix{Float64}(N,M)
    for i in 1:N, j in 1:M
        compute_kernel_value(dvs_train[i], dvs_test[j], pathway)
        if isnan(k[i,j])
            k[i,j] = 1e-7
        end
    end
    k
end

function compute_kernel_value(a::DataView, b::DataView, pathway::Union{Pathway, Void})
    kernel_params = prepare_kernel(a,b,pathway)
    kernel_function(kernel_params...)
end

# general view type vs. NAViewType
function prepare_kernel(a::DataView{K, <:ValueViewType}, b::DataView{K, NAViewType}, pathway::Union{Pathway, Void}) where {K<:KeyType}
    keys_to_use = pathway == nothing ? collect(keys(a.measurements)) : extract_pathway_keys(a, b, pathway)
    (map(m -> a.measurements[m].normalized_value, keys_to_use), zeros(length(keys_to_use)))
end

# RNASeqCall vs. NAViewType
function prepare_kernel(a::DataView{K, RNASeqCall}, b::DataView{K, NAViewType}, pathway::Union{Pathway{K}, Void}) where {K<:KeyType}
    keys_to_use = pathway == nothing ? collect(keys(a.measurements)) : extract_pathway_keys(a, b, pathway)
    (map(m -> a.measurements[m].normalized_value, keys_to_use), [falses(length(keys_to_use))...])
end

# NAViewType vs. general view type
function prepare_kernel(a::DataView{K, NAViewType}, b::DataView{K, <:ValueViewType}, pathway::Union{Pathway, Void}) where {K<:KeyType}
    (x_prime, x) = prepare_kernel(b,a, pathway)
    (x, x_prime)
end

# general view type vs. general view type
function prepare_kernel(a::DataView{K, V}, b::DataView{K, V}, pathway::Union{Pathway, Void}) where {K<:KeyType, V <:ValueViewType}
    common_keys = get_common_keys(a, b, pathway)
    x = map(k -> a.measurements[k].normalized_value, common_keys)
    x_prime = map(k -> b.measurements[k].normalized_value, common_keys)
    if length(common_keys) == 0
        x = [-2.]
        x_prime = [2.]
    end
    (x, x_prime)
end

function prepare_kernel(a::DataView{K, ExomeSeq}, b::DataView{K, ExomeSeq}, pathway::Union{Pathway{K}, Void})  where {K<:KeyType}
    common_keys = get_common_keys(a, b, pathway)
    x = Vector{Float64}()
    x_prime = Vector{Float64}()
    # x = Vector{Bool}()
    # x_prime = Vector{Bool}()
    factor_exact_matches = 1.
    for k in common_keys
        # check for match in protein change (this happens rarely and we want to emphasize this)
        # if there is match, make sure to include only the ExomeSeq measurement containing it
        measurements_a = a.measurements[k]
        measurements_b = b.measurements[k]

        pcs_a = map(es -> es.protein_change, measurements_a)
        pcs_b = map(es -> es.protein_change, measurements_b)
        # for pc in pcs_a
        #     if pc != "None"
        #         push!(x, true)
        #         break
        #     end
        # end
        # for pc in pcs_b
        #     if pc != "None"
        #         push!(x_prime, true)
        #         break
        #     end
        # end

        protein_change_intersection = intersect(pcs_a, pcs_b)
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
    if length(common_keys) == 0
        x = [-2.]
        x_prime = [2.]
    end
    (x, x_prime, Float64(length(x) * factor_exact_matches))
end

function prepare_kernel(a::DataView{K, ExomeSeq}, b::DataView{K, NAViewType}, pathway::Union{Pathway{K}, Void}) where {K<:KeyType}
    x = map(k -> mean(get_measurement_value(a, k)), [a.used_keys...])
    x_prime = zeros(length(x))
    (x, x_prime)
end

prepare_kernel(a::DataView{K,V}, b::DataView{K,V}, pathway::Union{Pathway, Void}) where {K <: KeyType, V <: NAViewType} = ([randn()],[randn()])

# function get_common_keys(a::DataView{K,V}, b::DataView{K,V}) where {K<:KeyType,V<:ViewType}
#     common_keys = Vector{K}()
#     for key_a in collect(keys(a.measurements))
#         if haskey(b.measurements, key_a)
#             push!(common_keys, key_a)
#         end
#     end
#     common_keys
# end

function get_common_keys(a::DataView{K,V}, b::DataView{K,V}, pathway::Union{Pathway, Void}) where {K<:KeyType,V<:ViewType}
    if pathway == nothing 
        if haskey(a.common_keys, b.cell_line_id)
            return a.common_keys[b.cell_line_id]
        end
        common_keys = Vector{K}()
        for key_a in collect(keys(a.measurements))
            if haskey(b.measurements, key_a)
                push!(common_keys, key_a)
            end
        end
        a.common_keys[b.cell_line_id] = common_keys
        b.common_keys[a.cell_line_id] = common_keys
        return common_keys
    end
    extract_pathway_keys(a, b, pathway)
end

function extract_pathway_keys(a::DataView, b::DataView, pathway::Pathway)
    keys_to_use = Vector{KeyType}()
    keys_in_a = b.view_type != NAViewType && haskey(a.common_keys, b.cell_line_id) ? a.common_keys[b.cell_line_id] : collect(keys(a.measurements))
    for key in pathway.genes
        if key ∈ keys_in_a
            push!(keys_to_use, key)
        end
    end
    keys_to_use
end


function extract_pathway_keys(a::DataView{<:KeyType, RPPA}, b::DataView{<:KeyType, RPPA}, pathway::Pathway)
    keys_to_use = Vector{KeyType}()
    pathway_hgnc_ids = map(k -> k.hgnc_id, pathway.genes)
    keys_in_a = haskey(a.common_keys, b.cell_line_id) ? a.common_keys[b.cell_line_id] : collect(keys(a.measurements))
    hgnc_ids_in_a = map(g -> g.hgnc_id, keys_in_a)
    for hgnc in pathway_hgnc_ids
        if hgnc ∈ hgnc_ids_in_a
            found = findin(hgnc_ids_in_a, hgnc)
            if length(found) > 0
                idx = found[1]
                push!(keys_to_use, keys_in_a[idx])
            end
        end
    end
    keys_to_use
end

#actual kernel functions
# Gaussian
kernel_function(x::Vector{Float64}, x_prime::Vector{Float64}, l::Float64 = Float64(length(x))) = exp(-dot(x-x_prime, x-x_prime)/l)
# Jaccard 
kernel_function(x::Vector{Bool}, x_prime::Vector{Bool}) = x'*x_prime/(x'*x + x_prime'*x_prime - x'*x_prime)
