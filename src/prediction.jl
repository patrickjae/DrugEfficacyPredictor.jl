function test(dep::DrugEfficacyPrediction)
    m = dep.model
    mse = 0
    test_cell_lines = unique(union(map(outcome -> collect(keys(outcome.outcome_values)), collect(values(dep.experiment.test_results)))))
    # do predictions for all drugs we saw at training time
    rankings = Dict{Drug, Vector{Int64}}()
    predictions = Dict{Drug, Vector{Float64}}()
    for (t, drug) in enumerate(keys(dep.experiment.results))
        training_drug_outcome = dep.experiment.results[drug]

        # TODO: check if we need to test for drug, i.e. if we have test_result for it
        G = Vector{Vector{Float64}}(length(dep.cross_kernels[drug]))
        exp_a = expected_value(m.a[t])
        for (k, c_kernel) in enumerate(dep.cross_kernels[drug])
            G[k] = c_kernel * exp_a
        end

        y_mean = sum(G .* expected_value.(m.e)) .+ expected_value(m.b[t])
        
        #########################
        # compute the ranking
        #########################
        # the ranking is relative to each drug, indicating which cell line is how likely to be affected by the drug, i.e we rank the cell lines, not the drugs
        # a lower rank means that the drug has higher effect on the cell line
        # WARNING: this is somewhat against the final use, since we are interested in ranking drugs conditioned on cell lines/tissues
        # lower IC50 is better, so lowest ic50 (or whatever measure) results in highest rank (actually, the lowest number)
        # !!!!!!!!  BUT we are dealing with the -log(IC50) so higher value results in higher rank
        #########################
        rankings[drug] = sortperm(y_mean, rev=true) 
        y_mean_rescaled = y_mean .* training_drug_outcome.outcome_std .+ training_drug_outcome.outcome_mean
        predictions[drug] = y_mean_rescaled

        actual_outcomes = dep.test_targets[drug]
        mse += sum((actual_outcomes .- y_mean_rescaled).^2)/m.N[t]
    end
    # @info view_weights=expected_value.(m.e))
    (mse, predictions, rankings)
end


function predict_outcomes(dep::DrugEfficacyPrediction, 
            cell_lines::Vector{CellLine}, 
            base_kernels::Union{Nothing, Dict{Type{<:ViewType}, Matrix{Float64}}}=nothing,
            pathway_specific_kernels::Union{Nothing, Dict{Type{<:ViewType}, Vector{Matrix{Float64}}}}=nothing)
    m = dep.model
    predictions = Dict{Drug, Vector{Float64}}()
    training_cell_lines = collect(values(dep.experiment.cell_lines))  

    if base_kernels == nothing || pathway_specific_kernels == nothing
        (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(dep.experiment, training_cell_lines, cell_lines)

        # base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
        # for v in dep.experiment.views
        #     dataviews_to_predict = map_data_views(cell_lines, v)
        #     dataviews_from_training = map_data_views(training_cell_lines, v)
        #     # dataviews_from_training = map(cl -> cl.views[v], collect(keys(dep.experiment.results[drug].outcome_values)))
        #     k = compute_kernel(dataviews_from_training, dataviews_to_predict)
        #     base_kernels[v] = k
        #     # @info "computed kernels for view $v, num data views to predict: $(length(dataviews_to_predict)), data views in training: $(length(dataviews_from_training)), size of kernel: $(size(k))"
        # end
    end
    # do predictions for all drugs we saw at training time
    for (t, drug) in enumerate(keys(dep.experiment.results))
        # compute similarity with all other cell lines in the training set
        kernels = Vector{Matrix{Float64}}()
        dataviews_from_training_idx = findin(training_cell_lines, collect(keys(dep.experiment.results[drug].outcome_values)))

        for v in dep.experiment.views
            # dataviews_from_training = map(cl -> cl.views[v], collect(keys(dep.experiment.results[drug].outcome_values)))
            k = base_kernels[v][dataviews_from_training_idx,:]
            push!(kernels, k)
            for pw_kernel in pathway_specific_kernels[v]
                push!(kernels, pw_kernel[dataviews_from_training_idx,:])
            end
            # @info "computed kernels for view $v, size of kernel: $(size(k))"
        end
        # TODO: compute kernel combinations

        a_expected = expected_value(m.a[t])

        G = Vector{Vector{Float64}}(length(kernels))
        exp_a = expected_value(m.a[t])
        for (k, kernel) in enumerate(kernels)
            G[k] = kernel' * exp_a
            # if k == 1
            #     println("G[$k] : $(typeof(G[k]))")
            #     display(kernel)
            # end
        end
        # G = transpose.(kernels) .* expected_value(m.a[t]) # results in a vector of scalar intermediate results
        # println("G: $(size(G)), Gs: $(size.(G))")
        e_expected = expected_value.(m.e)
        # println("e: $(size(e_expected)), e type: $(typeof(e_expected)), $e_expected")
        pred_y = sum(G .* e_expected) + expected_value(m.b[t])
        # rescale the normalized prediction
        pred_y_rescaled = pred_y * dep.experiment.results[drug].outcome_std + dep.experiment.results[drug].outcome_mean
        predictions[drug] = pred_y_rescaled
    end
    (predictions, base_kernels, pathway_specific_kernels)
end

function write_prediction_file(filename::String, dep::DrugEfficacyPrediction, predictions::Dict{Drug, Vector{Float64}})
    ranks = Dict{Drug, Vector{Int64}}()
    for d in collect(keys(predictions))
        ranks[d] = sortperm(predictions[d], rev=true)
    end
    f = open(filename, "w")
    @printf(f, "DrugAnonID,Drug1,Drug2,Drug3,Drug4,Drug5,Drug6,Drug7,Drug8,Drug9,Drug10,Drug11,Drug12,Drug13,Drug14,Drug15,Drug16,Drug17,Drug18,Drug19,Drug20,Drug21,Drug22,Drug23,Drug24,Drug25,Drug26,Drug27,Drug28,Drug29,Drug30,Drug31\n")
    for (cl_id, cl) in enumerate(collect(values(dep.experiment.cell_lines)))
       @printf(f, "%s", cl.id)
       for d_id in 1:length(dep.experiment.drugs)
           drug = dep.experiment.drugs["Drug$d_id"]
           @printf(f, ",%d", ranks[drug][cl_id])
       end
       @printf(f, "\n")
   end
   close(f)
end