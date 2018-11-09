function test(dep::DrugEfficacyPrediction, m::PredictionModel)
    mse = 0
    test_cell_lines = unique(union(map(outcome -> collect(keys(outcome.outcome_values)), collect(values(dep.experiment.test_results)))))
    # do predictions for all drugs we saw at training time
    rankings = Dict{Drug, Vector{Int64}}()
    predictions = Dict{Drug, Vector{Float64}}()
    sum_data_points = 0
    for (t, drug) in enumerate(keys(dep.experiment.results))
        training_drug_outcome = dep.experiment.results[drug]

        # TODO: check if we need to test for drug, i.e. if we have test_result for it
        G = Vector{Vector{Float64}}(undef, length(dep.cross_kernels[drug]))
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
        rankings[drug] = zeros(Int64, length(y_mean))
        rankings[drug][sortperm(y_mean, rev=true)] = collect(1:length(y_mean))
        y_mean_rescaled = y_mean .* training_drug_outcome.outcome_std .+ training_drug_outcome.outcome_mean
        predictions[drug] = y_mean_rescaled

        #get unnormalized test targets
        actual_outcomes = dep.test_targets[drug]
        mse += sum((actual_outcomes .- y_mean_rescaled).^2)
        sum_data_points += length(actual_outcomes)
    end
    mse /= sum_data_points
    # @info view_weights=expected_value.(m.e))
    (mse, predictions, rankings)
end


function predict_outcomes(dep::DrugEfficacyPrediction, m::PredictionModel,
            cell_lines::Vector{CellLine})
    predictions = Dict{Drug, Vector{Float64}}()
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
        training_cell_lines = collect(keys(dep.experiment.results[drug].outcome_values))

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
        # rescale the normalized prediction
        pred_y_rescaled = pred_y * dep.experiment.results[drug].outcome_std + dep.experiment.results[drug].outcome_mean
        predictions[drug] = pred_y_rescaled
    end
    predictions
end

function write_results(parent_dir::String, filename::String, dep::DrugEfficacyPrediction, predictions::Dict{Drug, Vector{Float64}}, m::PredictionModel)
    # write predictions and true values
    mkpath(joinpath(parent_dir,"prediction"))
    prediction_file = joinpath(parent_dir, "prediction", filename)
    f = open(prediction_file, "w")
    @printf(f, "DrugAnonID")
    for i in 1:length(dep.experiment.drugs)
        @printf(f, "Drug%d\t", i)
    end
    @printf(f, "\n")
    for (cl_id, cl) in enumerate(collect(values(dep.experiment.cell_lines)))
        @printf(f, "%s", cl.id)
        for d_id in 1:length(dep.experiment.drugs)
            drug = dep.experiment.drugs["Drug$d_id"]
            target = 0.
            if haskey(dep.experiment.results, drug) && haskey(dep.experiment.results[drug].outcome_values, cl)
                target = dep.experiment.results[drug].outcome_values[cl]
            elseif haskey(dep.experiment.test_results, drug) && haskey(dep.experiment.test_results[drug].outcome_values, cl)
                target = dep.experiment.test_results[drug].outcome_values[cl]
            end
            @printf(f, "\t%.5f(%.5f)", predictions[drug][cl_id], target)
        end
        @printf(f, "\n")
    end
    close(f)
    # write rankings
    mkpath(joinpath(parent_dir,"ranking"))
    ranking_file = joinpath(parent_dir, "ranking", filename)
    ranks = Dict{Drug, Vector{Int64}}()
    # we predict neg log values, i.e. higher value means lower concentration 
    # which means higher susceptibility of the cell line to the drug and thus higher rank (lower number)
    # hence, we sort in reverse order
    for d in collect(keys(predictions))
        ranks[d] = zeros(Int64, length(predictions[d]))
        ranks[d][sortperm(predictions[d], rev=true)] = collect(1:length(predictions[d]))
    end
    f = open(ranking_file, "w")
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

    # write model
    mkpath(joinpath(parent_dir, "model"))
    model_file = joinpath(parent_dir, "model", filename)
    f = open(model_file, "w")
    # T, K, N
    @printf(f, "T\tK\n")
    @printf(f, "%d\t%d\n", m.T, m.K)
    @printf(f, "N\n")
    [@printf(f, "%d\t", n) for n in m.N[1:end-1]]
    @printf(f, "%d\n", m.N[end])
    # gamma
    @printf(f, "ɣ\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.ɣ[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.ɣ[end]))
    # b
    @printf(f, "b\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.b[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.b[end]))
    # a
    @printf(f, "a\n")
    for a_t in m.a
        @printf(f, "%s\n", a_t.var_name)
        a_t_exp = expected_value(a_t)
        [@printf(f, "%.5f\t", val) for val in a_t_exp[1:end-1]]
        @printf(f, "%.5f\n", a_t_exp[end])
    end
    # lambda
    @printf(f, "λ\n")
    for (t, λ_t) in enumerate(m.λ)
        λ_t_exp = expected_value.(λ_t)
        @printf(f, "λ[%d]:\t", t)
        for val in λ_t_exp[1:end-1]
            @printf(f, "%.5f\t", val)
        end
        @printf(f, "%.5f\n", λ_t_exp[end])
    end
    # epsilon
    @printf(f, "ε\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.ε[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.ε[end]))
    # nu
    @printf(f, "ν\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.ν[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.ν[end]))
    # G
    @printf(f, "G\n")
    for t in 1:m.T, k in 1:m.K
        exp_g = expected_value(m.G[t,k])
        @printf(f, "%s:\t", m.G[t,k].var_name)
        [@printf(f, "%.5f\t", val) for val in exp_g[1:end-1]]
        @printf(f, "%.5f\n", exp_g[end])
    end
    # omega
    @printf(f, "⍵\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.⍵[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.⍵[end]))
    # e
    @printf(f, "e\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in m.e[1:end-1]]
    @printf(f, "%.5f\n", expected_value(m.e[end]))

    close(f)

end