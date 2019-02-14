function train(experiment_id::String, data::Dict{String, Any})
    log_message("in train method, retrieving experiment object")
    experiment = get_experiment(experiment_id)
    log_message("got experiment object")
    # try
        # extract inference config
        inference_config = read_inference_configuration(data)
        log_progress(experiment_id, "read inference configuration")
        # create and initialize the model
        pm = PredictionModel(experiment, inference_config)

        model_config = nothing
        # if do_gridsearch is false, extract the model config
        if !inference_config.do_gridsearch
            model_config = read_model_configuration(pm.model_type, pm, data)
        end
        log_progress(experiment_id, "read model configuration")

        # dep = create_drug_efficacy_predictor(experiment, inference_config)
        inference_config.target_dir = mktempdir()

        log_progress(experiment_id, "running model computations")
        if inference_config.do_cross_validation
            run_cross_validation(pm, inference_config, model_config)
        else
            run_model(pm, inference_config, model_config)
        end
        log_progress(experiment_id, "preparing results for download")
        # compile the results
        tar_dir = mktempdir()
        tar_file = joinpath(tar_dir, "results.tgz")
        log_message("archiving results at $tar_file")
        tar_results_cmd = `tar czf $(tar_file) $(inference_config.target_dir)`
        run(tar_results_cmd)
        # remove the output files
        rm(inference_config.target_dir, force = true, recursive = true)
        result_file_dictionary[experiment.internal_id] = tar_file
        log_progress(experiment_id, "finished, results ready for download")
    # catch ex
    #     st = map(string, stacktrace(catch_backtrace()))
    #     # log_progress(experiment, "exception has occurred: $(typeof(ex))")
    #     log_progress(experiment, "stacktrace: $(join(st, "\n"))")
    #     @info "caught exception..." st
    # end
    tar_file
end

function predict(experiment_id::String, data::Dict{String, Any})
    inference_config = read_inference_configuration(data)
    log_progress(experiment_id, "read inference configuration")
    if inference_config.do_gridsearch
        throw(ArgumentError("Prediction not implemented for gridsearch"))
    end
    pm = PredictionModel(experiment, inference_config)

    model_config = read_model_configuration(pm, data)
    log_progress(experiment_id, "read model configuration")

    log_progress(experiment_id, "training model")
    (lls, errs, test_errs, params, convergence) = inference(pm, model_config, inference_config = inference_config)
    log_progress(experiment_id, "predicting efficacy on provided cell lines")
    prediction_cell_lines = create_cell_lines(experiment, data, for_prediction = true)
    (predictions, ranks) = predict_outcomes(pm, params, prediction_cell_lines, held_out = true)

    result_dictionary = Dict{String, Any}()
    for (t, drug) in enumerate(keys(pm.data.results))
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
