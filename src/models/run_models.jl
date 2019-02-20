function run_cross_validation(pm::PredictionModel, inference_config::InferenceConfiguration, model_config::Union{Nothing, ModelConfiguration}; num_folds::Int64 = 10)
    if !inference_config.do_cross_validation
        throw(ArgumentError("run_cross_validation called without proper inference configuration, make sure the do_cross_validation flag is set to true"))
    end
    # running cross validation means, we are discarding the original training set/test set designations
    num_cl = length(pm.data.cell_lines)
    cls_per_fold = Int64(floor(num_cl/num_folds))

    for i in 1:num_folds
        log_progress(pm.data.internal_id, "starting computation of fold $i in CV")
        log_message("starting fold $i")
        start_idx = (i-1)*cls_per_fold + 1
        end_idx = i == num_folds ? num_cl : i * cls_per_fold

        # mark cell lines as test and training set
        for (cl_idx, cl) in enumerate(collect(values(pm.data.cell_lines)))
            if start_idx ≤ cl_idx ≤ end_idx
                cl.in_test_set = true
            else
                cl.in_test_set = false
            end
        end
        inference_config.fold_num = i
        log_message("reshuffled training/test set and set fold to $i")
        run_model(pm, inference_config, model_config)
    end
end

function run_model(pm::PredictionModel, inference_config::InferenceConfiguration, model_config::Union{Nothing, ModelConfiguration})
    # things to do only once (or once per fold)
    # for BMTMKL, this computes kernels and cross kernels depending on which cell lines
    # are in the training and test set
    log_message("running post_init!")
    post_init!(pm, model_config)
    # result dir is the target directory amended by the fold
    result_dir = inference_config.do_cross_validation ? joinpath(inference_config.target_dir, string(inference_config.fold_num)) : inference_config.target_dir
    mkpath(result_dir)

    if inference_config.compute_wpc_index
        # preparation for wpc index computations

        # 1) create files with test responses and training responses respectively (for this fold)
        # training
        training_response_file = create_response_file(pm, false, result_dir, "training_response.txt")
        # test
        test_response_file = create_response_file(pm, true, result_dir, "test_response.txt")

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
        run(wpc_script_gs_cmd)
    end
    all_model_configs = ModelConfiguration[]
    # if doing gridsearch, create model configs to use
    if inference_config.do_gridsearch
        gamma_dist_alphas = [250., 500., 750., 1000.]
        alpha_beta_ratios = [.1, .5, 1., 2.5, 5.]
        normal_vars = [.1, .5, 1., 2.]
        # normal_means = [1.]
        for alpha in gamma_dist_alphas, ratio in alpha_beta_ratios, v in normal_vars
            mc = ModelConfiguration(alpha, alpha/ratio, 1., v, pm)
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

    # do the actual inference sweep
    cnt = Base.Threads.Atomic{Int64}(0)
    Threads.@threads for i in 1:length(all_model_configs)
        mc = all_model_configs[i]
        # TODO: this assumes we are using the same values for all params, maybe change later
        # same for the generic filenam below
        alpha = mc.parameters["α_ɣ"]
        beta = mc.parameters["β_ɣ"]
        mu = mc.parameters["μ_e"]
        v = mc.parameters["σ_e"]
        Base.Threads.atomic_add!(cnt, 1)
        log_message("parameter setting $(cnt[]) (of $(length(all_model_configs))): alpha=$alpha, beta=$beta mu=$mu var=$v")
        # TODO: run each model multiple times (e.g. 10) and collect prediction results
        # report mean and variance of those predictions
        # log_message("calling parameter_inference method")
        (lls, errs, test_errs, params, convergence) = inference(pm, mc, inference_config = inference_config)
        log_message("inference stats: likelihood=$(lls[end]) train_error=$(errs[end]) test_error=$(test_errs[end]) convergence=$convergence")
        # models[i] = model
        error_string = inference_config.do_cross_validation ?
            "$alpha\t$beta\t$mu\t$v\t$(inference_config.fold_num)\t$(errs[end])\t$(test_errs[end])\t$(lls[end])\n" :
            "$alpha\t$beta\t$mu\t$v\t$(errs[end])\t$(test_errs[end])\t$(lls[end])\n"
        all_errors[i] = error_string

        (p, r) = predict_outcomes(pm, params, collect(values(pm.data.cell_lines)))
        log_message("done predicting")
        # predictions[i] = p
        # ranks[i] = r


        generic_filename =  "alpha_$(alpha)_beta_$(beta)_mean_$(mu)_var_$(v).txt"
        # writing result files
        write_results(pm, params, result_dir, generic_filename, p, r)

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

    # write the errors file
    f = open(joinpath(inference_config.target_dir, "errors.txt"), "a")
    for s in all_errors
        @printf(f, "%s", s)
    end
    close(f)

end
