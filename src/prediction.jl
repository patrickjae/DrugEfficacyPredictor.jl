function test(dep::DrugEfficacyPrediction, m::PredictionModel)
    mse = 0

    rankings = Dict{Drug, Vector{Int64}}()
    predictions = Dict{Drug, Vector{Float64}}()
    sum_data_points = 0
    for (t, drug) in enumerate(keys(dep.experiment.results))

        G = Vector{Vector{Float64}}(undef, length(dep.cross_kernels[drug]))
        exp_a = expected_value(m.a[t])
        for (k, c_kernel) in enumerate(dep.cross_kernels[drug])
            G[k] = c_kernel * exp_a
        end

        y_mean = sum([G[k] .* expected_value(m.e[k]) for k in 1:length(G)]) .+ expected_value(m.b[t])

        #########################
        # compute the ranking
        #########################
        # the ranking is relative to each drug, indicating which cell line is how likely to be affected by the drug, i.e we rank the cell lines, not the drugs
        # a lower rank means that the drug has higher effect on the cell line
        # WARNING: this is somewhat against the final use, since we are interested in ranking drugs conditioned on cell lines/tissues
        # lower IC50 is better, so lowest ic50 (or whatever measure) results in highest rank (actually, the lowest number)
        # !!!!!!!!  BUT we are dealing with the -log(IC50) so higher value results in higher rank
        #########################
        rankings[drug] = zeros(Int64, length(y_mean))
        rankings[drug][sortperm(y_mean, rev=true)] = collect(1:length(y_mean))

        y_mean_rescaled = y_mean .* dep.experiment.results[drug].outcome_std .+ dep.experiment.results[drug].outcome_mean
        predictions[drug] = y_mean_rescaled

        # @info "######### Drug $(drug.id) #########" targets=dep.targets[drug] predictions=y_mean_rescaled target_ranking prediction_ranking=rankings[drug] prediction_variance=(1 ./ expected_value(m.ε[t])) gamma=expected_value(m.ɣ[t]) nu=expected_value(m.ν[t])

        #get unnormalized test targets
        actual_outcomes = dep.test_targets[drug]
        mse += sum((actual_outcomes .- y_mean_rescaled).^2)
        sum_data_points += length(actual_outcomes)
        # @info "current test error setting" actual_outcomes y_mean_rescaled mse=sum((actual_outcomes .- y_mean_rescaled).^2)/length(actual_outcomes)
    end
    mse /= sum_data_points
    # @info view_weights=expected_value.(m.e))
    (mse, predictions, rankings)
end


function predict_outcomes(dep::DrugEfficacyPrediction, m::PredictionModel,
            cell_lines::Vector{CellLine})
    predictions = Dict{Drug, Vector{Float64}}()
    ranks = Dict{Drug, Vector{Int64}}()
    # all cell lines
    all_cell_lines = collect(values(dep.experiment.cell_lines))



    # if base_kernels == nothing || pathway_specific_kernels == nothing
    #     (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(dep.experiment, training_cell_lines, cell_lines)

    #     # base_kernels = Dict{Type{<:ViewType}, Matrix{Float64}}()
    #     # for v in dep.experiment.views
    #     #     dataviews_to_predict = map_data_views(cell_lines, v)
    #     #     dataviews_from_training = map_data_views(training_cell_lines, v)
    #     #     # dataviews_from_training = map(cl -> cl.views[v], collect(keys(dep.experiment.results[drug].outcome_values)))
    #     #     k = compute_kernel(dataviews_from_training, dataviews_to_predict)
    #     #     base_kernels[v] = k
    #     #     # @info "computed kernels for view $v, num data views to predict: $(length(dataviews_to_predict)), data views in training: $(length(dataviews_from_training)), size of kernel: $(size(k))"
    #     # end
    # end
    # do predictions for all drugs we saw at training time
    for (t, drug) in enumerate(keys(dep.experiment.results))
        # find cell lines that were used for training
        training_cell_lines = filter(cl -> !cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), all_cell_lines)

        training_set_cell_line_idx = findall((in)(training_cell_lines), all_cell_lines)
        predict_cell_line_idx = findall((in)(cell_lines), all_cell_lines)
        # compute similarity with all other cell lines in the training set
        kernels = Vector{Matrix{Float64}}()
        # dataviews_from_training_idx = findin(training_cell_lines, collect(keys(dep.experiment.results[drug].outcome_values)))

        for v in dep.experiment.views
            # dataviews_from_training = map(cl -> cl.views[v], collect(keys(dep.experiment.results[drug].outcome_values)))
            k = dep.base_kernels[v][training_set_cell_line_idx,predict_cell_line_idx]
            push!(kernels, k)
            if length(dep.experiment.pathway_information) != 0
                for pw_kernel in dep.pathway_specific_kernels[v]
                    push!(kernels, pw_kernel[training_set_cell_line_idx,predict_cell_line_idx])
                end
            end
            # @info "computed kernels for view $v, size of kernel: $(size(k))"
        end

        gene_expression = dep.base_kernels[GeneExpression][training_set_cell_line_idx, predict_cell_line_idx]
        methylation = dep.base_kernels[Methylation][training_set_cell_line_idx, predict_cell_line_idx]
        cnv = dep.base_kernels[CNV][training_set_cell_line_idx, predict_cell_line_idx]

        push!(kernels, gene_expression .* methylation)
        push!(kernels, gene_expression .* cnv)
        push!(kernels, cnv .* methylation)
        push!(kernels, gene_expression .* methylation .* cnv)

        a_expected = expected_value(m.a[t])

        G = Vector{Vector{Float64}}(undef, length(kernels))
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
        pred_y = sum(G .* e_expected) .+ expected_value(m.b[t])

        #compute ranking of predicted cell line responses
        ranks[drug] = zeros(Int64, length(pred_y))
        ranks[drug][sortperm(pred_y, rev=true)] = collect(1:length(pred_y))

        # rescale the normalized prediction
        pred_y_rescaled = pred_y * dep.experiment.results[drug].outcome_std .+ dep.experiment.results[drug].outcome_mean
        predictions[drug] = pred_y_rescaled
    end
    (predictions, ranks)
end
