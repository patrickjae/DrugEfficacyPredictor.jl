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
function create_drug_efficacy_predictor(experiment::Experiment, ::Dict{String, Any}=Dict{String, Any}(); subsume_pathways::Bool=true, do_variance_filtering::Bool=true)
    # make sure that measurements are normalized across cell lines
    if !experiment.is_normalized
        tt = @elapsed normalize_data_views(experiment)
        @info "normalizing training data took $tt seconds"
    end

    #do variance filtering if required
    if do_variance_filtering
        ft = @elapsed filter_data_views(experiment)
        @info "filtering data tool $ft seconds"
    end

    # collect outcomes statistics
    for (t, drug) in enumerate(keys(experiment.results))
        vals = collect(values(experiment.results[drug].outcome_values))
        experiment.results[drug].outcome_mean = mean(vals)
        experiment.results[drug].outcome_std = stdm(vals, experiment.results[drug].outcome_mean)
    end


    all_drugs = collect(keys(experiment.results)) # the tasks

    # compute base kernels between cell lines once, one for each ViewType
    # dep_base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
    # dep_pathway_specific_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    # grouped_data_kernels = Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    # dep_base_kernels = Dict{Drug, Vector{Matrix{Float64}}}()
    # dep_pathway_specific_kernels = Dict{Drug, Vector{Matrix{Float64}}}()
    # grouped_data_kernels = Dict{Drug, Vector{Matrix{Float64}}}()

    cell_lines = collect(values(experiment.cell_lines))
    # find cell lines that are present in all views, i.e. we need all data views for a cell line
    cell_lines_in_views = filter(cl -> length(cl.views) == length(experiment.views), cell_lines)

    (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(experiment, cell_lines, subsume_pathways=subsume_pathways)
    T = length(experiment.results)
    N = Vector{Int64}(undef, T)

    # compute drug specific kernels
    kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    cross_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # cross_base_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # cross_pathway_specific_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    test_targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    for (t, drug) in enumerate(all_drugs)
        #init
        # dep_base_kernels[drug] = Vector{Matrix{Float64}}()
        # dep_pathway_specific_kernels[drug] = Vector{Matrix{Float64}}()
        # grouped_data_kernels[drug] = Vector{Matrix{Float64}}()
        # cross_base_kernels[drug] = Vector{Matrix{Float64}}()
        # cross_pathway_specific_kernels[drug] = Vector{Matrix{Float64}}()
        kernels[drug] = Vector{Matrix{Float64}}()
        cross_kernels[drug] = Vector{Matrix{Float64}}()

        result_cell_lines = collect(keys(experiment.results[drug].outcome_values))
        test_result_cell_lines = collect(keys(experiment.test_results[drug].outcome_values))
        # find cell lines for which we have outcomes for the current drug
        # indices of the cell lines with outcome that are available for all views (in outcome structure)
        # idx_in_outcome = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines_in_views)
        idx_in_cell_lines = findall((in)(result_cell_lines), cell_lines)
        # find cell lines in test set for which we a measurement for the current drug
        # idx_in_test_outcome = findin(collect(keys(experiment.test_results[drug].outcome_values)), cell_lines_in_views)
        test_idx_in_cell_lines = findall((in)(test_result_cell_lines), cell_lines)
        # indices of present cell lines that have outcome data (in cell lines vector)
        # idx_in_cell_lines = findin(cell_lines_in_views, collect(keys(experiment.results[drug].outcome_values)))
        idx_in_outcome = findall((in)(cell_lines), result_cell_lines)
        # indices of present cell lines that have a test outcome
        # test_idx_in_cell_lines = findin(cell_lines_in_views, collect(keys(experiment.test_results[drug].outcome_values)))
        idx_in_test_outcome = findall((in)(cell_lines), test_result_cell_lines)
        # only deal with cell lines that are available for all views AND have an outcome for the current drug
        N[t] = length(idx_in_cell_lines)

        for v in experiment.views
            # cell_lines_with_view = filter(cl -> haskey(cl.views, v), cell_lines)
            # @info "cell lines that have view $v: $(length(cell_lines_with_view))"
            # idx_in_outcome2 = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines)
            # @info "view $v: $(length(idx_in_outcome)), overall: $(length(idx2)), N[t]: $(model.N[t])"

            # add base kernel restricted to cell lines for which we have outcome
            # push!(dep_base_kernels[drug], base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            # push!(cross_base_kernels[drug], base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])
            push!(kernels[drug], base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            push!(cross_kernels[drug], base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])

            for pw_kernel in pathway_specific_kernels[v]
                # push!(dep_pathway_specific_kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                # push!(cross_pathway_specific_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
                push!(kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                push!(cross_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
            end
        end
        # compute additional kernels
        # grouped kernels
        gene_expression = base_kernels[GeneExpression][idx_in_cell_lines, idx_in_cell_lines]
        methylation = base_kernels[Methylation][idx_in_cell_lines, idx_in_cell_lines]
        cnv = base_kernels[CNV][idx_in_cell_lines, idx_in_cell_lines]

        push!(kernels[drug], gene_expression .* methylation)
        push!(kernels[drug], gene_expression .* cnv)
        push!(kernels[drug], cnv .* methylation)
        push!(kernels[drug], gene_expression .* methylation .* cnv)

        cross_gene_expression = base_kernels[GeneExpression][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_methylation = base_kernels[Methylation][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_cnv = base_kernels[CNV][test_idx_in_cell_lines, idx_in_cell_lines]

        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_cnv)
        push!(cross_kernels[drug], cross_cnv .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation .* cross_cnv)

        #normalize the outcome data data for each drug across the cell lines
        targets[drug] = ((collect(values(experiment.results[drug].outcome_values))[idx_in_outcome] .- experiment.results[drug].outcome_mean)./experiment.results[drug].outcome_std)
        # test outcomes, not normalized
        test_targets[drug] = collect(values(experiment.test_results[drug].outcome_values))[idx_in_test_outcome]
    end
    # accommodate for additional mixed kernels
    K += 4
    # pm = PredictionModel(T, K, N)
    dep = DrugEfficacyPrediction(experiment, T, K, N)
    dep.base_kernels = base_kernels
    dep.pathway_specific_kernels = pathway_specific_kernels
    dep.kernels = kernels
    dep.cross_kernels = cross_kernels
    dep.targets = targets
    dep.test_targets = test_targets
    # for drug in all_drugs
    #     # dep.base_kernels[drug] = dep_base_kernels[drug]
    #     # dep.pathway_specific_kernels[drug] = dep_pathway_specific_kernels[drug]
    #     dep.kernels[drug] = kernels[drug]
    #     dep.targets[drug] = targets[drug]
    #     # dep.cross_base_kernels[drug] = cross_base_kernels[drug]
    #     # dep.cross_pathway_specific_kernels[drug] = cross_pathway_specific_kernels[drug]
    #     dep.cross_kernels[drug] = cross_kernels[drug]
    #     dep.test_targets[drug] = test_targets[drug]
    # end
    dep
end

function compute_all_kernels(experiment::Experiment, cell_lines::Vector{CellLine}, cell_lines_test::Union{Vector{CellLine}, Nothing}=nothing; subsume_pathways::Bool=true)
    base_kernels = OrderedDict{Type{<:ViewType}, Matrix{Float64}}()
    pathway_specific_kernels = OrderedDict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    K = 0
    # compute a base kernel for each view, containing all cell lines that are present in all views
    num_views = length(experiment.views)
    @info "computing base kernels"
    for v in experiment.views
        # get views of this type for all cell lines
        # v = experiment.views[v_id]
        data_views = [map_data_views(cell_lines, v)]
        if cell_lines_test != nothing
            push!(data_views, map_data_views(cell_lines_test, v))
        end
        # compute the base kernel for this view
        kct = @elapsed base_kernels[v] = compute_kernel(data_views...)
        @info "computing base kernel for $v $kct seconds)"
        K += 1
        @info "computing pathway specific kernels..."
        pathway_specific_kernels[v] = Vector{Matrix{Float64}}()

        # # pathway specific kernels method 1: compute a kernel for each pathway
        if !subsume_pathways
            num_pws = length(experiment.pathway_information)
            tt = 0
            for pathway in experiment.pathway_information
                tt += @elapsed push!(pathway_specific_kernels[v], compute_kernel(data_views..., pathway=pathway))
            end
            K += num_pws

            # method 1 end
        else
            # pathway specific kernels method 2:
            # - for each cell line:
            #   - for each pathway:
            #       - compute average (expression data) or max value (otherwise) over pathway genes
            # - compute cell line similarities as usual
            tt = @elapsed begin
                cell_line_gene_set_views = Vector{Vector{Float64}}(undef, length(cell_lines))
                for (cl_idx, dv) in enumerate(data_views...)
                    cell_line_gene_set_views[cl_idx] = zeros(length(experiment.pathway_information))
                    for (pw_idx, pathway) in enumerate(experiment.pathway_information)

                        (m_values, _) = prepare_kernel(dv, dv, pathway)
                        if v in [RNASeq, GeneExpression]
                            cell_line_gene_set_views[cl_idx][pw_idx] = mean(m_values)
                        else
                            cell_line_gene_set_views[cl_idx][pw_idx] = maximum(m_values)
                        end
                    end
                end
                pw_kernel = zeros(length(cell_lines), length(cell_lines))
                for i in 1:length(cell_lines), j in 1:length(cell_lines)
                    pw_kernel[i,j] = kernel_function(cell_line_gene_set_views[i], cell_line_gene_set_views[j])
                end
                push!(pathway_specific_kernels[v], pw_kernel)
                K += 1
            end
        # # method 2 end
        end        
        @info "computing pathway specific kernels for $v took $tt seconds)"

    end
    @info "done computing kernels"
    (K, base_kernels, pathway_specific_kernels)
end


function map_data_views(cell_lines::Vector{CellLine}, view_type::Type{<:ViewType})
    ret = Vector{DataView{<:KeyType, <:ViewType}}()
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

get_complete_cell_lines(experiment::Experiment) = filter(cl -> length(cl.views) == length(experiment.views), collect(values(experiment.cell_lines)))

function compute_kernel(dvs::Vector{DataView}; pathway::Union{Pathway, Nothing} = nothing)
    N = length(dvs)
    k = zeros(N,N)
    prep_time = 0
    comp_time = 0
    kernel_compute_time = @elapsed begin
        Threads.@threads for i in 1:N
            for  j in i:N
                kpt, kct, k[i,j] = compute_kernel_value(dvs[i], dvs[j], pathway)
                prep_time += kpt
                comp_time += kct
                if isnan(k[i,j])
                    k[i,j] = 1e-7
                end
                if i≠j k[j,i] = k[i,j] end
            end
        end
        if !isposdef(k) 
            h_k = Array(Hermitian(k))
            # @warn "kernel is not positive definite, using only upper triangular form" isposdef(h_k) sum_diff=sum(abs.(h_k - k))
            k = h_k
        end
    end
    @info "computed kernel ($(N*(N-1)) computations)" preparation_time=prep_time computation_time=comp_time
    # @info "kernel computations took $kernel_compute_time seconds"
    k
end

function compute_kernel(dvs_train::Vector{DataView}, dvs_test::Vector{DataView}; pathway::Union{Pathway, Nothing}=nothing)
    N = length(dvs_train)
    M = length(dvs_test)
    k = zeros(N,M)
    Threads.@threads for i in 1:N
        for j in 1:M
            kpt, kct, k[i,j] = compute_kernel_value(dvs_train[i], dvs_test[j], pathway)
            if isnan(k[i,j])
                k[i,j] = 1e-7
            end
        end
    end
    k
end

function compute_kernel_value(a::DataView, b::DataView, pathway::Union{Pathway, Nothing})
    kpc = @elapsed kernel_params = prepare_kernel(a,b,pathway)
    kcc = @elapsed kfv = kernel_function(kernel_params...)
    kpc, kcc, kfv
end

# general view type vs. NAViewType
function prepare_kernel(a::DataView{K, <:ValueViewType}, b::DataView{K, NAViewType}, pathway::Union{Pathway, Nothing}) where {K<:KeyType}
    keys_to_use = get_common_keys(a, b, pathway)
    (map(m -> a.measurements[m].normalized_value, keys_to_use), zeros(length(keys_to_use)))
end

# NAViewType vs. general view type
function prepare_kernel(a::DataView{K, NAViewType}, b::DataView{K, <:ValueViewType}, pathway::Union{Pathway, Nothing}) where {K<:KeyType}
    (x_prime, x) = prepare_kernel(b,a, pathway)
    (x, x_prime)
end

# RNASeqCall vs. NAViewType
function prepare_kernel(a::DataView{K, RNASeqCall}, b::DataView{K, NAViewType}, pathway::Union{Pathway{K}, Nothing}) where {K<:KeyType}
    keys_to_use = get_common_keys(a, b, pathway)
    (map(m -> a.measurements[m].normalized_value, keys_to_use), [falses(length(keys_to_use))...])
end

# general view type vs. general view type
function prepare_kernel(a::DataView{K, V}, b::DataView{K, V}, pathway::Union{Pathway, Nothing}) where {K<:KeyType, V <:ValueViewType}
    common_keys = get_common_keys(a, b, pathway)
    # hack: when no keys are shared, we have small similarity
    if length(common_keys) == 0
        x = [-2.]
        x_prime = [2.]
        return (x, x_prime  )
    end
    x = map(k -> a.measurements[k].normalized_value, common_keys)
    x_prime = map(k -> b.measurements[k].normalized_value, common_keys)
    (x, x_prime)
end

function prepare_kernel(a::DataView{K, ExomeSeq}, b::DataView{K, ExomeSeq}, pathway::Union{Pathway{K}, Nothing})  where {K<:KeyType}
    common_keys = get_common_keys(a, b, pathway)
    x = Vector{Float64}()
    x_prime = Vector{Float64}()
    factor_exact_matches = 1.
    for k in common_keys
        # check for match in protein change (this happens rarely and we want to emphasize this)
        # if there is match, make sure to include only the ExomeSeq measurement containing it
        measurements_a = a.measurements[k]
        measurements_b = b.measurements[k]

        pcs_a = map(es -> es.protein_change, measurements_a)
        pcs_b = map(es -> es.protein_change, measurements_b)

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

function prepare_kernel(a::DataView{K, ExomeSeq}, b::DataView{K, NAViewType}, pathway::Union{Pathway{K}, Nothing}) where {K<:KeyType}
    x = map(k -> mean(get_measurement_value(a, k)), [a.used_keys...])
    x_prime = zeros(length(x))
    (x, x_prime)
end

prepare_kernel(a::DataView{K,V}, b::DataView{K,V}, pathway::Union{Pathway, Nothing}) where {K <: KeyType, V <: NAViewType} = ([randn()],[randn()])

function get_common_keys(a::DataView, b::DataView, pathway::Union{Pathway, Nothing})
    # @info "getting common keys" pathway
    if pathway != nothing 
        return extract_pathway_keys(a, b, pathway)
    end
    if (b.view_type == NAViewType) || (b.used_keys == a.used_keys)
        return collect(a.used_keys)
    elseif haskey(a.common_keys, b.cell_line_id)
        return a.common_keys[b.cell_line_id]
    elseif haskey(b.common_keys, a.cell_line_id)
        return b.common_keys[a.cell_line_id]
    end
    common_keys = collect(intersect(a.used_keys, b.used_keys))
    a.common_keys[b.cell_line_id] = common_keys
    b.common_keys[a.cell_line_id] = common_keys
    common_keys
end

function extract_pathway_keys(a::DataView, b::DataView, pathway::Pathway)
    # if b.view_type == NAViewType
    #     return intersect(pathway.genes, a.used_keys)
    # end
    # if haskey(a.common_keys, b.cell_line_id)
    #     # @info "found precomputed key intersection"
    #     return intersect(pathway.genes, a.common_keys[b.cell_line_id])
    # end
    # if haskey(b.common_keys, a.cell_line_id)
    #     # @info "found precomputed key intersection"
    #     return intersect(pathway.genes, b.common_keys[a.cell_line_id])
    # end
    common_keys = get_common_keys(a,b,nothing)
    intersect(pathway.genes, common_keys)
end


function extract_pathway_keys(a::DataView{<:KeyType, RPPA}, b::DataView{<:KeyType, RPPA}, pathway::Pathway)
    common_keys = get_common_keys(a,b,nothing)

    keys_to_use = Vector{KeyType}()
    pathway_hgnc_ids = map(k -> k.hgnc_id, pathway.genes)
    hgnc_ids_in_common_keys = map(g -> g.hgnc_id, common_keys)

    hgnc_intersection = intersect(pathway_hgnc_ids, hgnc_ids_in_common_keys)
    key_idx = findall((in)(hgnc_intersection), hgnc_ids_in_common_keys)
    common_keys[key_idx]
end

#actual kernel functions
# Gaussian
kernel_function(x::Vector{Float64}, x_prime::Vector{Float64}, l::Float64 = Float64(length(x))) = exp(-dot(x-x_prime, x-x_prime)/l)
# Jaccard 
kernel_function(x::Vector{Bool}, x_prime::Vector{Bool}) = x'*x_prime/(x'*x + x_prime'*x_prime - x'*x_prime)
