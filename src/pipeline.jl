function log_progress(experiment::Experiment, s::String)
    @spawnat server_process push!(training_progress[experiment.internal_id], s)
end


function train(experiment_id::String, data::Dict{String, Any})
    experiment = get_experiment(experiment_id)
    log_message("in train method")
    # try
        # extract inference config
        inference_config = read_inference_configuration(data)
        log_progress(experiment, "read inference configuration")
        model_config = nothing
        # if do_gridsearch is false, extract the model config
        if !inference_config.do_gridsearch
            model_config = read_model_configuration(data)
        end
        log_progress(experiment, "read model configuration")
        # run the model
        dep = create_drug_efficacy_predictor(experiment, inference_config)
        inference_config.target_dir = mktempdir()

        log_progress(experiment, "running model computations")
        if inference_config.do_cross_validation
            run_cross_validation(dep, inference_config = inference_config, model_config = model_config)
        else
            run_model(dep, inference_config = inference_config, model_config = model_config)
        end
        log_progress(experiment, "preparing results for download")
        # compile the results
        tar_dir = mktempdir()
        tar_file = joinpath(tar_dir, "results.tgz")
        log_message("archiving results at $tar_file")
        tar_results_cmd = `tar czf $(tar_file) $(inference_config.target_dir)`
        run(tar_results_cmd)
        # remove the output files
        rm(inference_config.target_dir, force = true, recursive = true)
        result_file_dictionary[experiment.internal_id] = tar_file
        log_progress(experiment, "finished, results ready for download")
    # catch ex
    #     st = map(string, stacktrace(catch_backtrace()))
    #     # log_progress(experiment, "exception has occurred: $(typeof(ex))")
    #     log_progress(experiment, "stacktrace: $(join(st, "\n"))")
    #     @info "caught exception..." st
    # end
    tar_file
end

function predict(experiment::Experiment, data::Dict{String, Any})
    inference_config = read_inference_configuration(data)
    log_progress(experiment, "read inference configuration")
    if inference_config.do_gridsearch
        throw(ArgumentError("Prediction not implemented for gridsearch"))
    end
    model_config = read_model_configuration(data)
    log_progress(experiment, "read model configuration")
    dep = create_drug_efficacy_predictor(experiment, inference_config)

    log_progress(experiment, "training model")
    (lls, errs, test_errs, model, convergence) = parameter_inference(dep, inference_config = inference_config, model_config = model_config)
    log_progress(experiment, "predicting efficacy on provided cell lines")
    prediction_cell_lines = create_cell_lines(experiment, data, for_prediction = true)
    (predictions, ranks) = predict_outcomes(dep, model, prediction_cell_lines)
    result_dictionary = Dict{String, Any}()
    for (t, drug) in enumerate(keys(dep.experiment.results))
        result_dictionary[drug.id] = Dict{String, Any}()
        for (cl_idx, cl) in enumerate(prediction_cell_lines)
            result_dictionary[drug.id][cl.id] = Dict{String, Any}()
            result_dictionary[drug.id][cl.id]["prediction"] = predictions[drug][cl_idx]
            result_dictionary[drug.id][cl.id]["rank"] = ranks[drug][cl_idx]
        end
    end
    log_progress(experiment, "finished")
    result_dictionary
end


function run_cross_validation(dep::DrugEfficacyPrediction; num_folds::Int64 = 10, inference_config::InferenceConfiguration = InferenceConfiguration(), model_config::Union{Nothing, ModelConfiguration} = nothing)
    if !inference_config.do_cross_validation
        throw(ArgumentError("run_cross_validation called without proper inference configuration, make sure the do_cross_validation flag is set to true"))
    end
    # running cross validation means, we are discarding the original training set/test set designations
    num_cl = length(dep.experiment.cell_lines)
    cls_per_fold = Int64(floor(num_cl/num_folds))

    for i in 1:num_folds
        log_progress(dep.experiment, "starting computation of fold $i in CV")
        log_message("starting fold $i")
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
        inference_config.fold_num = i
        run_model(dep, inference_config = inference_config, model_config = model_config)
    end
end

function create_response_file(dep::DrugEfficacyPrediction, test::Bool, target_dir::String, filename::String)
    response_file = joinpath(target_dir, filename)
    f = open(response_file, "w")
    @printf(f, "CellLine")
    drugs = collect(keys(dep.experiment.results))
    for drug in drugs
        @printf(f, "\t%s", drug.id)
    end
    @printf(f, "\n")

    for cl in collect(values(dep.experiment.cell_lines))
        # if these two are equal
        if test == cl.in_test_set
            @printf(f, "%s", cl.id)
            for drug in drugs
                if haskey(dep.experiment.results[drug].outcome_values, cl)
                    @printf(f, "\t%.9f", dep.experiment.results[drug].outcome_values[cl])
                else
                    @printf(f, "\tNA")
                end
            end
            @printf(f, "\n")
        end
    end
    close(f)
    response_file
end

function run_model(dep; inference_config::InferenceConfiguration = InferenceConfiguration(), model_config::Union{Nothing, ModelConfiguration} = nothing)
    # things to do only once (or once per fold)
    # compute kernels and cross kernels depending on which cell lines are in the training and test set
    set_training_test_kernels(dep)
    # result dir is the target directory amended by the fold
    result_dir = inference_config.do_cross_validation ? joinpath(inference_config.target_dir, string(inference_config.fold_num)) : inference_config.target_dir
    mkpath(result_dir)

    if inference_config.compute_wpc_index
        # preparation for wpc index computations

        # 1) create files with test responses and training responses respectively (for this fold)
        # training
        training_response_file = create_response_file(dep, false, result_dir, "training_response.txt")
        # test
        test_response_file = create_response_file(dep, true, result_dir, "test_response.txt")

        # 2) create random predictions from these
        script_name = joinpath(PROJECT_ROOT, "analytics", "create_random_predictions.pl")
        rand_pred_dir = joinpath(result_dir, "random_predictions")
        mkpath(rand_pred_dir)
        create_rand_pred_cmd = `perl $script_name $training_response_file $test_response_file $rand_pred_dir 100`
        run(create_rand_pred_cmd)
        # 3) calculate z-scores
        script_name = joinpath(PROJECT_ROOT, "analytics", "calculate_zscore_weights.pl")
        sd_file = joinpath(PROJECT_ROOT, "analytics", "pooled_sd.txt")
        gold_standard_file = joinpath(PROJECT_ROOT, "analytics", "gold_standard.csv")
        zscore_file = joinpath(result_dir, "zscores.txt")
        calc_zscore_cmd = `perl $script_name $rand_pred_dir $sd_file $test_response_file $gold_standard_file $zscore_file`
        run(calc_zscore_cmd)
        # 4) with random prediction and zscores, compute the wpc index
        script_name = joinpath(PROJECT_ROOT, "analytics", "weighted_average_concordance_index.pl")
        mkpath(joinpath(result_dir, "wpc"))
        wpc_target_file = joinpath(result_dir, "wpc", "gold_standard_wpc.txt")

        wpc_summary_file = joinpath(inference_config.target_dir, "wpc_summary.txt")

        file = open(wpc_summary_file, "w")
        @printf(file, "ranking file\taverage probabilistic c-index\tweighted average probabilistic c-index\twpc p-value\n")
        close(file)

        wpc_script_gs_cmd = `perl $script_name $sd_file $test_response_file $zscore_file $gold_standard_file $wpc_target_file $wpc_summary_file`
    end
    all_model_configs = ModelConfiguration[]
    # if doing gridsearch, create model configs to use
    if inference_config.do_gridsearch
        gamma_dist_alphas = [250., 500., 750., 1000.]
        alpha_beta_ratios = [.1, .5, 1., 2.5, 5.]
        normal_vars = [.1, .5, 1., 2.]
        # normal_means = [1.]
        for alpha in gamma_dist_alphas, ratio in alpha_beta_ratios, v in normal_vars
            mc = ModelConfiguration(alpha, alpha/ratio, 1., v)
            # assume no bias, hence zero-mean priors on b
            mc.μ_b = 0.
            push!(all_model_configs, mc)
        end
        log_message("generated gridsearch configs ($(length(all_model_configs)))")
    #  else just use the model config as provided in the call
    else
        model_config == nothing ?
            throw(ArgumentError("not doing a grid search but no model configuration specified")) :
            push!(all_model_configs, model_config)
    end
    all_errors = Vector{String}(undef, length(all_model_configs))
    # predictions = Vector{Dict{Drug, Vector{Float64}}}(undef, length(all_model_configs))
    # ranks = Vector{Dict{Drug, Vector{Int64}}}(undef, length(all_model_configs))
    # models = Vector{PredictionModel}(undef, length(all_model_configs))
    # do the actual inference sweep
    cnt = Base.Threads.Atomic{Int64}(0)
    Threads.@threads for i in 1:length(all_model_configs)
        mc = all_model_configs[i]
        # TODO: this assumes we are using the same values for all params, maybe change later
        # same for the generic filenam below
        alpha = mc.⍺_ɣ
        beta = mc.β_ɣ
        mu = mc.μ_e
        v = mc.𝜎_e
        Base.Threads.atomic_add!(cnt, 1)
        log_message("parameter setting $(cnt[]) (of $(length(all_model_configs))): alpha=$alpha, beta=$beta mu=$mu var=$v")
        # TODO: run each model multiple times (e.g. 10) and collect prediction results
        # report mean and variance of those predictions
        # log_message("calling parameter_inference method")
        (lls, errs, test_errs, model, convergence) = parameter_inference(dep, inference_config = inference_config, model_config = mc)
        log_message("inference stats: likelihood=$(lls[end]) train_error=$(errs[end]) test_error=$(test_errs[end]) convergence=$convergence")
        # models[i] = model
        error_string = inference_config.do_cross_validation ?
            "$alpha\t$beta\t$mu\t$v\t$(inference_config.fold_num)\t$(errs[end])\t$(test_errs[end])\t$(lls[end])\n" :
            "$alpha\t$beta\t$mu\t$v\t$(errs[end])\t$(test_errs[end])\t$(lls[end])\n"
        all_errors[i] = error_string

        (p, r) = predict_outcomes(dep, model, collect(values(dep.experiment.cell_lines)))
        # predictions[i] = p
        # ranks[i] = r


        generic_filename =  "alpha_$(alpha)_beta_$(beta)_mean_$(mu)_var_$(v).txt"
        # writing result files
        write_results(dep, result_dir, generic_filename, p, r, model)

        resulting_ranking = joinpath(result_dir, "ranking", generic_filename)

        if inference_config.compute_wpc_index
            wpc_script_cmd = `perl $script_name $sd_file $test_response_file $zscore_file $resulting_ranking $wpc_target_file $wpc_summary_file`
            run(wpc_script_cmd)
        end

        p = nothing
        r = nothing
        model = nothing
        lls = nothing
        errs = nothing
        test_errs = nothing

        Base.GC.gc(true)
    end

    # for i in 1:length(all_model_configs)
    #     mc = all_model_configs[i]
    #     alpha = mc.⍺_ɣ
    #     beta = mc.β_ɣ
    #     mu = mc.μ_e
    #     v = mc.𝜎_e
    #
    # end

    # write the errors file
    f = open(joinpath(inference_config.target_dir, "errors.txt"), "a")
    for s in all_errors
        @printf(f, "%s", s)
    end
    close(f)

end

function set_training_test_kernels(dep::DrugEfficacyPrediction)
    # collect outcomes statistics, only on results used for training
    for (t, drug) in enumerate(keys(dep.experiment.results))
        vals = Float64[]
        [if !cl.in_test_set push!(vals, dep.experiment.results[drug].outcome_values[cl]) end for cl in keys(dep.experiment.results[drug].outcome_values)]
        dep.experiment.results[drug].outcome_mean = mean(vals)
        dep.experiment.results[drug].outcome_std = stdm(vals, dep.experiment.results[drug].outcome_mean)
    end

    all_drugs = collect(keys(dep.experiment.results)) # the tasks
    cell_lines = collect(values(dep.experiment.cell_lines)) # all cell lines in the experiment, excluding held-out data

    # the kernels between cell lines
    kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # the target values
    targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()

    # cross kernels between training data and test data
    cross_kernels = DataStructures.OrderedDict{Drug, Vector{Matrix{Float64}}}()
    # target values for test data
    test_targets = DataStructures.OrderedDict{Drug, Vector{Float64}}()


    for (t, drug) in enumerate(all_drugs)
        #init
        kernels[drug] = Vector{Matrix{Float64}}()
        cross_kernels[drug] = Vector{Matrix{Float64}}()

        # determine cell lines available for this drug
        result_cell_lines = filter((cl) ->  !cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), cell_lines)
        test_result_cell_lines = filter((cl) -> cl.in_test_set && haskey(dep.experiment.results[drug].outcome_values, cl), cell_lines)

        # get their indicex in original cell line array
        idx_in_cell_lines = findall((in)(result_cell_lines), cell_lines)
        test_idx_in_cell_lines = findall((in)(test_result_cell_lines), cell_lines)

        dep.N[t] = length(idx_in_cell_lines)

        #construct the kernels
        for v in dep.experiment.views
            push!(kernels[drug], dep.base_kernels[v][idx_in_cell_lines, idx_in_cell_lines])
            push!(cross_kernels[drug], dep.base_kernels[v][test_idx_in_cell_lines, idx_in_cell_lines])
            if length(dep.experiment.pathway_information) != 0
                for pw_kernel in dep.pathway_specific_kernels[v]
                    push!(kernels[drug], pw_kernel[idx_in_cell_lines, idx_in_cell_lines])
                    push!(cross_kernels[drug], pw_kernel[test_idx_in_cell_lines, idx_in_cell_lines])
                end
            end
        end
        # compute additional kernels
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

        train_outcome = map(cl -> dep.experiment.results[drug].outcome_values[cl], result_cell_lines)
        test_outcome = map(cl -> dep.experiment.results[drug].outcome_values[cl], test_result_cell_lines)

        #normalize the outcome data for each drug across the cell lines
        targets[drug] = (train_outcome .- dep.experiment.results[drug].outcome_mean)./dep.experiment.results[drug].outcome_std

        # test outcomes, not normalized
        test_targets[drug] = test_outcome
    end
    # recompute overall number of kernels (+4 accomodates for combined kernels)
    dep.K = length(dep.base_kernels) + sum(map(length, collect(values(dep.pathway_specific_kernels)))) + 4

    dep.kernels = kernels
    dep.cross_kernels = cross_kernels
    dep.targets = targets
    dep.test_targets = test_targets

    dep
 end



function write_results(dep::DrugEfficacyPrediction, parent_dir::String, filename::String, predictions::Dict{Drug, Vector{Float64}}, ranks::Dict{Drug, Vector{Int64}}, m::PredictionModel)
    # write predictions and true values
    mkpath(joinpath(parent_dir,"prediction"))
    prediction_file = joinpath(parent_dir, "prediction", filename)
    f = open(prediction_file, "w")
    @printf(f, "DrugAnonID\tType")
    for d in collect(keys(dep.experiment.results))
        @printf(f, "\t%s (prediction)\t%s (measured)", d.id, d.id)
    end
    @printf(f, "\n")
    for (cl_id, cl) in enumerate(collect(values(dep.experiment.cell_lines)))
        @printf(f, "%s", cl.id)
        if cl.in_test_set @printf(f, "\tTest") else @printf(f, "\tTraining") end
        for drug in collect(keys(dep.experiment.results))
            # drug = dep.experiment.drugs["Drug$d_id"]
            target = 0.
            if haskey(dep.experiment.results, drug) && haskey(dep.experiment.results[drug].outcome_values, cl)
                target = dep.experiment.results[drug].outcome_values[cl]
            # elseif haskey(dep.experiment.test_results, drug) && haskey(dep.experiment.test_results[drug].outcome_values, cl)
            #     target = dep.experiment.test_results[drug].outcome_values[cl]
            end
            @printf(f, "\t%.5f\t%.5f", predictions[drug][cl_id], target)
        end
        @printf(f, "\n")
    end
    close(f)
    # write rankings
    mkpath(joinpath(parent_dir,"ranking"))
    ranking_file = joinpath(parent_dir, "ranking", filename)
    # ranks = Dict{Drug, Vector{Int64}}()
    # we predict neg log values, i.e. higher value means lower concentration
    # which means higher susceptibility of the cell line to the drug and thus higher rank (lower number)
    # hence, we sort in reverse order
    for d in collect(keys(predictions))
        ranks[d] = zeros(Int64, length(predictions[d]))
        ranks[d][sortperm(predictions[d], rev=true)] = collect(1:length(predictions[d]))
    end
    f = open(ranking_file, "w")
    @printf(f, "DrugAnonID")
    for drug in collect(keys(dep.experiment.results))
        @printf(f, ",%s", drug.id)
    end
    @printf(f, "\n")
    for (cl_id, cl) in enumerate(collect(values(dep.experiment.cell_lines)))
        @printf(f, "%s", cl.id)
        for drug in collect(keys(dep.experiment.results))
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
