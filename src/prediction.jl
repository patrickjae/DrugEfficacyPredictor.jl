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


function predict_outcomes(dep::DrugEfficacyPrediction, m::PredictionModel, cell_lines::Vector{CellLine}; held_out::Bool = false)
    predictions = Dict{Drug, Vector{Float64}}()
    ranks = Dict{Drug, Vector{Int64}}()
    # all cell lines present for training
    all_cell_lines = collect(values(dep.experiment.cell_lines))
    # only compute the kernels when dealing with held-out data
    num_cross_kernels = dep.K
    base_cross_kernels = dep.base_kernels
    pw_specific_cross_kernels = dep.pathway_specific_kernels
    if held_out
        (num_cross_kernels, base_cross_kernels, pw_specific_cross_kernels) = compute_all_kernels(dep.experiment, all_cell_lines, cell_lines, subsume_pathways = dep.subsume_pathways)
    end
    # do predictions for all drugs we saw at training time
    for (t, drug) in enumerate(keys(dep.experiment.results))
        # find cell lines that have outcome for this drug and are in the training set
        training_cell_lines = filter(cl -> !cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), all_cell_lines)
        training_set_cell_line_idx = findall((in)(training_cell_lines), all_cell_lines)
        # find the ids of these cell_lines
        predict_cell_line_idx = collect(1:length(cell_lines))
        if !held_out
            # these cell lines should be included in the experiment, determine their position
            predict_cell_line_idx = findall((in)(cell_lines), all_cell_lines)
        end

        # compute similarity with all other cell lines in the training set
        kernels = Vector{Matrix{Float64}}()

        for v in dep.experiment.views
            push!(kernels, base_cross_kernels[v][training_set_cell_line_idx, predict_cell_line_idx])
            if length(dep.experiment.pathway_information) != 0
                for pw_kernel in pw_specific_cross_kernels[v]
                    push!(kernels, pw_kernel[training_set_cell_line_idx, predict_cell_line_idx])
                end
            end
        end

        gene_expression = base_cross_kernels[GeneExpression][training_set_cell_line_idx, :]
        methylation = base_cross_kernels[Methylation][training_set_cell_line_idx, :]
        cnv = base_cross_kernels[CNV][training_set_cell_line_idx, :]

        push!(kernels, gene_expression .* methylation)
        push!(kernels, gene_expression .* cnv)
        push!(kernels, cnv .* methylation)
        push!(kernels, gene_expression .* methylation .* cnv)

        G = Vector{Vector{Float64}}(undef, length(kernels))
        exp_a = expected_value(m.a[t])
        [G[k] = kernel' * exp_a for (k, kernel) in enumerate(kernels)]

        pred_y = sum(G .* expected_value.(m.e)) .+ expected_value(m.b[t])

        #compute ranking of predicted cell line responses
        ranks[drug] = zeros(Int64, length(pred_y))
        ranks[drug][sortperm(pred_y, rev=true)] = collect(1:length(pred_y))

        # rescale the normalized prediction
        pred_y_rescaled = pred_y * dep.experiment.results[drug].outcome_std .+ dep.experiment.results[drug].outcome_mean
        predictions[drug] = pred_y_rescaled
    end
    (predictions, ranks)
end
