"""
Read in inference configuration.
"""
function read_inference_configuration(data::Dict{String, Any})
    ic = InferenceConfiguration()
	log_message("created inference config object")
    if haskey(data, "inference_configuration")
        ic_data = data["inference_configuration"]
        if haskey(ic_data, "model") ic.model_config = Meta.eval(Symbol(ic_data["model"])) end
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

function read_model_type(data::Dict{String, Any})
	log_message("in read_model_type")
	if haskey(data, "model")
		try
			model_type = eval(Meta.parse(data["model"]))
			log_message("model type is $model_type")
			return model_type
		catch
			log_message("unknown model type: $(mc_data["model"])")
			throw(ArgumentError("Unknown model: $(mc_data["model"])"))
		end
    end
	log_message("no model specified in $data")
    throw(ArgumentError("You must specify a model."))
end

function read_model_configuration(pm::PredictionModel, data::Dict{String, Any})
	log_message("read_model_configuration not implemented for $(typeof(pm))")
end

function init!(pm::PredictionModel, data::Dict{String, Any})
	log_message("init! method not implemented, skipping...")
end

function post_init!(pm::PredictionModel, mc::ModelConfiguration)
	log_message("post_init! method not implemented, skipping...")
end

function train(experiment_id::String, data::Dict{String, Any})
	try
	    (pm, inference_config, model_config) = create_prediction_model(experiment_id, data)
	    log_progress(experiment_id, "created prediction model")

	    # actual training procedure

	    # create a temp dir for storing the results
        temp_dir = mktempdir()
	    inference_config.target_dir = joinpath(temp_dir, "results")

	    log_progress(experiment_id, "running model computations")
	    # do cross validation or single model training
	    if inference_config.do_cross_validation
	        run_cross_validation(pm, inference_config, model_config)
	    else
	        run_model(pm, inference_config, model_config)
	    end

	    log_progress(experiment_id, "preparing results for download")
	    # compile the results into a tar ball
	    tar_dir = mktempdir()
	    tar_file = joinpath(tar_dir, "results.tgz")
	    log_message("archiving results at $tar_file")
	    tar_results_cmd = `tar czf $(tar_file) -C $temp_dir results`
	    run(tar_results_cmd)
	    # remove the output files
	    rm(temp_dir, force = true, recursive = true)
	    Utils.store_result(experiment_id, tar_file)
	    log_progress(experiment_id, "finished, results ready for download")
	catch ex
		exception_name = string(ex)
		Core.println(exception_name)
        st = map(string, stacktrace(catch_backtrace()))
		log_progress(experiment_id, exception_name)
		for ste in st
			Core.println(ste)
			log_progress(experiment_id, ste)
		end
	end
end

function predict(experiment_id::String, data::Dict{String, Any})
    (pm, inference_config, model_config) = create_prediction_model(experiment_id, data)

    if model_config == nothing
        throw(ArgumentError("Prediction not implemented for gridsearch"))
    end

    log_progress(experiment_id, "training model")
    # do one inference run on the model with provided parameters
    # TODO: do the training via additional api call, then store the most recent paramters object
    (lls, errs, test_errs, params, convergence) = inference(pm, model_config, inference_config = inference_config)

    # read cell lines from data for which we want a prediction
    log_progress(experiment_id, "predicting efficacy on provided cell lines")
    prediction_cell_lines = create_cell_lines(experiment, data, for_prediction = true)
    (predictions, ranks) = predict_outcomes(pm, params, prediction_cell_lines, held_out = true)

    # create a result object
    # TODO: do this asynchronously and download results?
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
