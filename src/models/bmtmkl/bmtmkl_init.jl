function read_model_configuration(pm::PredictionModel{BMTMKLModel}, data::Dict{String, Any})
	if haskey(data, "model_configuration")
		mc_data = data["model_configuration"]
		# general setting
		if haskey(mc_data, "alpha") && haskey(mc_data, "beta") && haskey(mc_data, "mu") && haskey(mc_data, "sigma2")
			return ModelConfiguration(pm, mc_data["alpha"], mc_data["beta"], mc_data["mu"], mc_data["sigma2"])
		else
			mc = ModelConfiguration(pm)
			# refined settings if present
			if haskey(mc_data, "bias_precision_alpha") && haskey(mc_data, "bias_precision_beta")
				mc.parameters["α_ɣ"] = mc_data["bias_precision_alpha"]
				mc.parameters["β_ɣ"] = mc_data["bias_precision_beta"]
			end
			if haskey(mc_data, "weights_precision_alpha") && haskey(mc_data, "weights_precision_beta")
				mc.parameters["α_λ"] = mc_data["weights_precision_alpha"]
				mc.parameters["β_λ"] = mc_data["weights_precision_beta"]
			end
			if haskey(mc_data, "outcome_precision_alpha") && haskey(mc_data, "outcome_precision_beta")
				mc.parameters["α_ε"] = mc_data["outcome_precision_alpha"]
				mc.parameters["β_ε"] = mc_data["outcome_precision_beta"]
			end
			if haskey(mc_data, "intermed_results_precision_alpha") && haskey(mc_data, "intermed_results_precision_beta")
				mc.parameters["α_ν"] = mc_data["intermed_results_precision_alpha"]
				mc.parameters["β_ν"] = mc_data["intermed_results_precision_beta"]
			end
			if haskey(mc_data, "kernel_weights_precision_alpha") && haskey(mc_data, "kernel_weights_precision_beta")
				mc.parameters["α_⍵"] = mc_data["kernel_weights_precision_alpha"]
				mc.parameters["β_⍵"] = mc_data["kernel_weights_precision_beta"]
			end
			if haskey(mc_data, "bias_mean") && haskey(mc_data, "bias_sigma2")
				mc.parameters["μ_b"] = mc_data["bias_mean"]
				mc.parameters["σ_b"] = mc_data["bias_sigma2"]
			end
			if haskey(mc_data, "kernel_weights_mean") && haskey(mc_data, "kernel_weights_sigma2")
				mc.parameters["μ_e"] = mc_data["kernel_weights_mean"]
				mc.parameters["σ_e"] = mc_data["kernel_weights_sigma2"]
			end
			if haskey(mc_data, "weights_mean") && haskey(mc_data, "weights_sigma2")
				mc.parameters["μ_a"] = mc_data["weights_mean"]
				mc.parameters["Σ_a"] = mc_data["weights_sigma2"]
			end
			if haskey(mc_data, "intermed_results_mean") && haskey(mc_data, "intermed_results_sigma2")
				mc.parameters["μ_g"] = mc_data["intermed_results_mean"]
				mc.parameters["Σ_g"] = mc_data["intermed_results_sigma2"]
			end
		end
	end
	# standard settings if no config present
	ModelConfiguration(pm)
end

function init!(m::PredictionModel{BMTMKLModel}, ic::InferenceConfiguration)
	(base_kernels, pathway_specific_kernels) = compute_all_kernels(m.data, collect(values(m.data.cell_lines)), subsume_pathways = ic.subsume_pathways)
	m.precomputations["base_kernels"] = base_kernels
	m.precomputations["pathway_specific_kernels"] = pathway_specific_kernels
end

function post_init!(m::PredictionModel{BMTMKLModel}, mc::ModelConfiguration)
    # collect outcomes statistics, only on results used for training
    for (t, drug) in enumerate(keys(m.data.results))
        vals = Float64[]
        [if !cl.in_test_set push!(vals, m.data.results[drug].outcome_values[cl]) end for cl in keys(m.data.results[drug].outcome_values)]
        m.data.results[drug].outcome_mean = mean(vals)
        m.data.results[drug].outcome_std = stdm(vals, m.data.results[drug].outcome_mean)
    end
    all_drugs = collect(keys(m.data.results)) # the tasks
    cell_lines = collect(values(m.data.cell_lines)) # all cell lines in the experiment, excluding held-out data

    # the kernels between cell lines
    kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # the target values
    targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    # cross kernels between training data and test data
    cross_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # target values for test data
    test_targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

	# get number of results per taks (drug)
	mc.parameters["N"] = DataStructures.OrderedDict{Drug,Int64}()

	# log_message("resetting N")
    for (t, drug) in enumerate(all_drugs)
        #init
        kernels[drug] = Vector{Matrix{Float64}}()
        cross_kernels[drug] = Vector{Matrix{Float64}}()

        # determine cell lines available for this drug
        result_cell_lines = filter((cl) ->  !cl.in_test_set && haskey(m.data.results[drug].outcome_values, cl), cell_lines)
        test_result_cell_lines = filter((cl) -> cl.in_test_set && haskey(m.data.results[drug].outcome_values, cl), cell_lines)

        # get their indicex in original cell line array
        idx_in_cell_lines = findall((in)(result_cell_lines), cell_lines)
        test_idx_in_cell_lines = findall((in)(test_result_cell_lines), cell_lines)

        mc.parameters["N"][drug] = length(idx_in_cell_lines)

        #construct the kernels
        for v in m.data.views
            push!(kernels[drug], m.precomputations["base_kernels"][v][idx_in_cell_lines, idx_in_cell_lines])
            push!(cross_kernels[drug], m.precomputations["base_kernels"][v][test_idx_in_cell_lines, idx_in_cell_lines])
            if length(m.data.pathway_information) != 0
                for pw_kernel in m.precomputations["pathway_specific_kernels"][v]
                    push!(kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                    push!(cross_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
                end
            end
        end
        # compute additional kernels
        gene_expression = m.precomputations["base_kernels"][GeneExpression][idx_in_cell_lines, idx_in_cell_lines]
        methylation = m.precomputations["base_kernels"][Methylation][idx_in_cell_lines, idx_in_cell_lines]
        cnv = m.precomputations["base_kernels"][CNV][idx_in_cell_lines, idx_in_cell_lines]

        push!(kernels[drug], gene_expression .* methylation)
        push!(kernels[drug], gene_expression .* cnv)
        push!(kernels[drug], cnv .* methylation)
        push!(kernels[drug], gene_expression .* methylation .* cnv)

        cross_gene_expression = m.precomputations["base_kernels"][GeneExpression][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_methylation = m.precomputations["base_kernels"][Methylation][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_cnv = m.precomputations["base_kernels"][CNV][test_idx_in_cell_lines, idx_in_cell_lines]

        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_cnv)
        push!(cross_kernels[drug], cross_cnv .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation .* cross_cnv)

        train_outcome = map(cl -> m.data.results[drug].outcome_values[cl], result_cell_lines)
        test_outcome = map(cl -> m.data.results[drug].outcome_values[cl], test_result_cell_lines)

        #normalize the outcome data for each drug across the cell lines
        targets[drug] = (train_outcome .- m.data.results[drug].outcome_mean) ./ m.data.results[drug].outcome_std

        # test outcomes, not normalized
        test_targets[drug] = test_outcome
    end

	m.precomputations["kernels"] = kernels
	m.precomputations["cross_kernels"] = cross_kernels
	m.precomputations["targets"] = targets
	m.precomputations["test_targets"] = test_targets
	nothing
 end


function compute_all_kernels(experiment::Experiment, cell_lines::Vector{CellLine}, cell_lines_test::Union{Vector{CellLine}, Nothing}=nothing; subsume_pathways::Bool=true)
	log_message("in compute kernel method")
	base_kernels = OrderedDict{Type{<:ViewType}, Matrix{Float64}}()
    pathway_specific_kernels = OrderedDict{Type{<:ViewType}, Vector{Matrix{Float64}}}()
    # compute a base kernel for each view, containing all cell lines that are present in all views
    num_views = length(experiment.views)
    for v in experiment.views
        # get views of this type for all cell lines
        # v = experiment.views[v_id]
        data_views = [map_data_views(cell_lines, v)]
        if cell_lines_test != nothing
            push!(data_views, map_data_views(cell_lines_test, v))
        end
        # compute the base kernel for this view
        kct = @elapsed base_kernels[v] = compute_kernel(data_views...)
        if cell_lines_test == nothing
            log_message("computing base kernel for $v $kct seconds)")
        end
        if length(experiment.pathway_information) != 0
            if cell_lines_test == nothing
                log_message("computing pathway specific kernels...")
            end
            pathway_specific_kernels[v] = Vector{Matrix{Float64}}()

            # # pathway specific kernels method 1: compute a kernel for each pathway
            if !subsume_pathways
                num_pws = length(experiment.pathway_information)
                tt = 0
                for pathway in collect(values(experiment.pathway_information))
                    tt += @elapsed push!(pathway_specific_kernels[v], compute_kernel(data_views..., pathway=pathway))
                end

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
                end
            # # method 2 end
            end
            if cell_lines_test == nothing
                log_message("computing pathway specific kernels for $v took $tt seconds)")
            end
        else
            if cell_lines_test == nothing
                log_message("Note: no pathway information provided")
            end
        end
    end
    if cell_lines_test == nothing
        log_message("done computing kernels")
    end
    (base_kernels, pathway_specific_kernels)
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
