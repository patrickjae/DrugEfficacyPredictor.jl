function run_cross_validation(dep::DrugEfficacyPrediction; num_folds = 10)
    num_cl = length(dep.experiment.cell_lines)
    cls_per_fold = Int64(floor(num_cl/num_folds))

    for i in 1:num_folds
        start_idx = (i-1)*cls_per_fold + 1
        end_idx = i==num_folds ? num_cl : i*cls_per_fold

        # mark cell lines as test and training set
        for (cl_idx, cl) in enumerate(collect(values(dep.experiment.cell_lines)))
            if start_idx ≤ cl_idx ≤ end_idx
                cl.in_test_set = true
            else
                cl.in_test_set = false
            end
        end

        set_training_test_kernels(dep)
    end
    # running cross validation means, we are discarding the original training set/test set designations

end


function set_training_test_kernels(dep::DrugEfficacyPrediction)
    # collect outcomes statistics
    for (t, drug) in enumerate(keys(dep.experiment.results))
        vals = Float64[]
        [if !cl.in_test_set push!(vals, dep.experiment.results[drug].outcome_values[cl]) end for cl in keys(dep.experiment.results[drug].outcome_values)]
        # vals = collect(values(experiment.results[drug].outcome_values))
        dep.experiment.results[drug].outcome_mean = mean(vals)
        dep.experiment.results[drug].outcome_std = stdm(vals, dep.experiment.results[drug].outcome_mean)
    end

    all_drugs = collect(keys(dep.experiment.results)) # the tasks
    cell_lines = collect(values(dep.experiment.cell_lines))
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
        result_cell_lines = filter((cl) ->  !cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), cell_lines)
        test_result_cell_lines = filter((cl) -> cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), cell_lines)
        # result_cell_lines = collect(keys(experiment.results[drug].outcome_values))
        # test_result_cell_lines = collect(keys(experiment.test_results[drug].outcome_values))
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
        dep.N[t] = length(idx_in_cell_lines)

        for v in dep.experiment.views
            # cell_lines_with_view = filter(cl -> haskey(cl.views, v), cell_lines)
            # @info "cell lines that have view $v: $(length(cell_lines_with_view))"
            # idx_in_outcome2 = findin(collect(keys(experiment.results[drug].outcome_values)), cell_lines)
            # @info "view $v: $(length(idx_in_outcome)), overall: $(length(idx2)), N[t]: $(model.N[t])"

            # add base kernel restricted to cell lines for which we have outcome
            # push!(dep_base_kernels[drug], base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            # push!(cross_base_kernels[drug], base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])
            push!(kernels[drug], dep.base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            push!(cross_kernels[drug], dep.base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])
            if length(dep.experiment.pathway_information) != 0
                for pw_kernel in dep.pathway_specific_kernels[v]
                    # push!(dep_pathway_specific_kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                    # push!(cross_pathway_specific_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
                    push!(kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                    push!(cross_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
                end
            end
        end
        # compute additional kernels
        # grouped kernels
        gene_expression = dep.base_kernels[GeneExpression][idx_in_cell_lines, idx_in_cell_lines]
        methylation = dep.base_kernels[Methylation][idx_in_cell_lines, idx_in_cell_lines]
        cnv = dep.base_kernels[CNV][idx_in_cell_lines, idx_in_cell_lines]

        push!(kernels[drug], gene_expression .* methylation)
        push!(kernels[drug], gene_expression .* cnv)
        push!(kernels[drug], cnv .* methylation)
        push!(kernels[drug], gene_expression .* methylation .* cnv)

        cross_gene_expression = dep.base_kernels[GeneExpression][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_methylation = dep.base_kernels[Methylation][test_idx_in_cell_lines, idx_in_cell_lines]
        cross_cnv = dep.base_kernels[CNV][test_idx_in_cell_lines, idx_in_cell_lines]

        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_cnv)
        push!(cross_kernels[drug], cross_cnv .* cross_methylation)
        push!(cross_kernels[drug], cross_gene_expression .* cross_methylation .* cross_cnv)

        #normalize the outcome data for each drug across the cell lines
        # train_outcome = Float64[]
        # test_outcome = Float64[]
        # [if haskey(dep.experiment.results[drug].outcome_values, cl) push!(train_outcome, dep.experiment.results[drug].outcome_values[cl]) end for cl in result_cell_lines]
        # [if haskey(dep.experiment.results[drug].outcome_values, cl) push!(testn_outcome, dep.experiment.results[drug].outcome_values[cl]) end for cl in test_result_cell_lines]
        train_outcome = map(cl -> dep.experiment.results[drug].outcome_values[cl], result_cell_lines)
        test_outcome = map(cl -> dep.experiment.results[drug].outcome_values[cl], test_result_cell_lines)

        # targets[drug] = ((collect(values(experiment.results[drug].outcome_values))[idx_in_outcome] .- experiment.results[drug].outcome_mean)./experiment.results[drug].outcome_std)
        targets[drug] = (train_outcome .- dep.experiment.results[drug].outcome_mean)./dep.experiment.results[drug].outcome_std
        # test outcomes, not normalized
        # test_targets[drug] = collect(values(experiment.test_results[drug].outcome_values))[idx_in_test_outcome]
        test_targets[drug] = test_outcome
    end
    # accommodate for additional mixed kernels
    dep.K = length(dep.base_kernels) + sum(map(length, collect(values(dep.pathway_specific_kernels)))) + 4

    dep.kernels = kernels
    dep.cross_kernels = cross_kernels
    dep.targets = targets
    dep.test_targets = test_targets

    dep
 end
