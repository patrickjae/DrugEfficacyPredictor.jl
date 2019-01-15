"""
Read in inference configuration.
"""
function read_inference_configuration(data::Dict{String, Any})
    ic = InferenceConfiguration()
    if haskey(data, "inference_configuration")
        ic_data = data["inference_configuration"]
        if haskey(ic_data, "convergence_criterion") ic.convergence_criterion = ic_data["convergence_criterion"] end
        if haskey(ic_data, "min_iter") ic.min_iter = ic_data["min_iter"] end
        if haskey(ic_data, "max_iter") ic.max_iter = ic_data["max_iter"] end
        if haskey(ic_data, "do_gridsearch") ic.do_gridsearch = ic_data["do_gridsearch"] end
        if haskey(ic_data, "do_cross_validation") ic.do_cross_validation = ic_data["do_cross_validation"] end
        if haskey(ic_data, "compute_wpc_index") ic.compute_wpc_index = ic_data["compute_wpc_index"] end
        if haskey(ic_data, "do_variance_filtering") ic.do_variance_filtering = ic_data["do_variance_filtering"] end
        if haskey(ic_data, "subsume_pathways") ic.subsume_pathways = ic_data["subsume_pathways"] end
    end
    ic
end

function read_model_configuration(data::Dict{String, Any})
    mc = ModelConfiguration()
    if haskey(data, "model_configuration")
        mc_data = data["model_configuration"]
        # general setting
        if haskey(mc_data, "alpha") && haskey(mc_data, "beta") && haskey(mc_data, "mu") && haskey(mc_data, "sigma2")
            return ModelConfiguration(mc_data["alpha"], mc_data["beta"], mc_data["mu"], mc_data["sigma2"])
        end
        # refined settings
        if haskey(mc_data, "bias_precision_alpha") && haskey(mc_data, "bias_precision_beta")
            mc.⍺_ɣ = mc_data["bias_precision_alpha"]
            mc.β_ɣ = mc_data["bias_precision_beta"]
        end
        if haskey(mc_data, "weights_precision_alpha") && haskey(mc_data, "weights_precision_beta")
            mc.⍺_λ = mc_data["weights_precision_alpha"]
            mc.β_λ = mc_data["weights_precision_beta"]
        end
        if haskey(mc_data, "outcome_precision_alpha") && haskey(mc_data, "outcome_precision_beta")
            mc.⍺_ε = mc_data["outcome_precision_alpha"]
            mc.β_ε = mc_data["outcome_precision_beta"]
        end
        if haskey(mc_data, "intermed_results_precision_alpha") && haskey(mc_data, "intermed_results_precision_beta")
            mc.⍺_ν = mc_data["intermed_results_precision_alpha"]
            mc.β_ν = mc_data["intermed_results_precision_beta"]
        end
        if haskey(mc_data, "kernel_weights_precision_alpha") && haskey(mc_data, "kernel_weights_precision_beta")
            mc.⍺_⍵ = mc_data["kernel_weights_precision_alpha"]
            mc.β_⍵ = mc_data["kernel_weights_precision_beta"]
        end
        if haskey(mc_data, "bias_mean") && haskey(mc_data, "bias_sigma2")
            mc.μ_b = mc_data["bias_mean"]
            mc.𝜎_b = mc_data["bias_sigma2"]
        end
        if haskey(mc_data, "kernel_weights_mean") && haskey(mc_data, "kernel_weights_sigma2")
            mc.μ_e = mc_data["kernel_weights_mean"]
            mc.𝜎_e = mc_data["kernel_weights_sigma2"]
        end
        if haskey(mc_data, "weights_mean") && haskey(mc_data, "weights_sigma2")
            mc.μ_a = mc_data["weights_mean"]
            mc.Σ_a = mc_data["weights_sigma2"]
        end
        if haskey(mc_data, "intermed_results_mean") && haskey(mc_data, "intermed_results_sigma2")
            mc.μ_g = mc_data["intermed_results_mean"]
            mc.Σ_g = mc_data["intermed_results_sigma2"]
        end
    end
    mc
end

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
function create_drug_efficacy_predictor(experiment::Experiment, ic::InferenceConfiguration)
    # make sure that measurements are normalized across cell lines
    if !experiment.is_normalized
        tt = @elapsed normalize_data_views(experiment)
        log_message("normalizing training data took $tt seconds")
    end

    #do variance filtering if required
    if ic.do_variance_filtering
        ft = @elapsed filter_data_views(experiment)
        log_message("filtering data took $ft seconds")
    end

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

    (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(experiment, cell_lines, subsume_pathways=ic.subsume_pathways)
    T = length(experiment.results)
    N = Vector{Int64}(undef, T)


    dep = DrugEfficacyPrediction(experiment, T, K, N)
    dep.base_kernels = base_kernels
    dep.pathway_specific_kernels = pathway_specific_kernels
    dep.subsume_pathways = ic.subsume_pathways

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
    log_message("computing base kernels")
    for v in experiment.views
        # get views of this type for all cell lines
        # v = experiment.views[v_id]
        data_views = [map_data_views(cell_lines, v)]
        if cell_lines_test != nothing
            push!(data_views, map_data_views(cell_lines_test, v))
        end
        # compute the base kernel for this view
        kct = @elapsed base_kernels[v] = compute_kernel(data_views...)
        log_message("computing base kernel for $v $kct seconds)")
        K += 1
        if length(experiment.pathway_information) != 0
            log_message("computing pathway specific kernels...")
            pathway_specific_kernels[v] = Vector{Matrix{Float64}}()

            # # pathway specific kernels method 1: compute a kernel for each pathway
            if !subsume_pathways
                num_pws = length(experiment.pathway_information)
                tt = 0
                for pathway in collect(values(experiment.pathway_information))
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
                    cell_line_gene_set_views = map((x) -> Vector{Vector{Float64}}(undef, length(x)), data_views)
                    for (cl_set_idx, cl_set) in enumerate(data_views)
                        for (cl_idx, dv) in enumerate(cl_set)
                            cell_line_gene_set_views[cl_set_idx][cl_idx] = zeros(length(experiment.pathway_information))
                            for (pw_idx, pathway) in enumerate(values(experiment.pathway_information))

                                (m_values, _) = prepare_kernel(dv, dv, pathway)
                                if v in [RNASeq, GeneExpression]
                                    cell_line_gene_set_views[cl_set_idx][cl_idx][pw_idx] = mean(m_values)
                                else
                                    cell_line_gene_set_views[cl_set_idx][cl_idx][pw_idx] = maximum(m_values)
                                end
                            end
                        end
                    end
                    len_views = length(cell_line_gene_set_views)
                    pw_kernel = zeros(length(cell_line_gene_set_views[1]), length(cell_line_gene_set_views[len_views]))
                    for i in 1:size(pw_kernel, 1), j in 1:size(pw_kernel, 2)
                        pw_kernel[i,j] = kernel_function(cell_line_gene_set_views[1][i], cell_line_gene_set_views[len_views][j])
                    end
                    push!(pathway_specific_kernels[v], pw_kernel)
                    K += 1
                end
            # # method 2 end
            end
            log_message("computing pathway specific kernels for $v took $tt seconds)")
        else
            log_message("Note: no pathway information provided")
        end
    end
    log_message("done computing kernels")
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
    log_message("computed kernel ($(N*(N-1)) computations), preparation time=$prep_time computation time=$comp_time")
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
