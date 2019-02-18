function create_prediction_model(experiment_id::String, data::Dict{String, Any})
    # get experiment object
    experiment = get_experiment(experiment_id)
    # extract inference config    log_message("reading inference config")
    inference_config = read_inference_configuration(data)
    prepare_data!(experiment, inference_config)
    # read the model type
    model_type = read_model_type(data)
    # create and initialize the model
    pm = PredictionModel(model_type, experiment)
    init!(pm, inference_config)
    # when not doing gridsearch, look for model config in request data
    model_config = nothing
    if !inference_config.do_gridsearch
        model_config = read_model_configuration(pm, data)
    end
    (pm, inference_config, model_config)
end

function prepare_data!(data::Experiment, ic::InferenceConfiguration)
    # make sure that measurements are normalized across cell lines
    if !data.is_normalized
        Data.normalize_data_views(data)
    end
    #do variance filtering if required
    if ic.do_variance_filtering
        Data.filter_data_views(data)
    else
        # make sure we are using all keys
        Data.unfilter_data_views(data)
    end
    # for each drug, store mean and standard deviation, used for normalization
    for (t, drug) in enumerate(keys(data.results))
        vals = collect(values(data.results[drug].outcome_values))
        data.results[drug].outcome_mean = mean(vals)
        data.results[drug].outcome_std = stdm(vals, data.results[drug].outcome_mean)
    end
end

function copy(o::ModelConfiguration)
    c = typeof(o)()
    for f in fieldnames(typeof(o))
        setfield!(c, f, getfield(o, f))
    end
    c
end
