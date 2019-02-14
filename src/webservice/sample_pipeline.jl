function load_iorio_data(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    host = "localhost"
    base_dir = joinpath(PROJECT_ROOT, "data", "iorio")
    force_reload = haskey(request_dictionary, "force_reload") ? request_dictionary["force_reload"] : false
    if experiment_exists("iorio")
        if force_reload
            delete_experiment("iorio")
        else
            req.response.status = 200
            req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Iorio data set already exists, returning experiment ID. If you need to reload, call with parameter force_reload = true", "experiment_id" => "iorio")))
            return req.response
        end
    end

    # create the experiment
    run(`curl -X POST http://$host:$port/experiments -d '{"experiment_id":"iorio"}'`)

    #import the cell lines
    for f in readdir(joinpath(base_dir, "cell_lines"))
        if !startswith(basename(f), "json") continue end
        filename = joinpath(base_dir, "cell_lines", f)
        cl_cmd = `curl -X POST http://$host:$port/experiments/iorio/data/cell_line -d @$filename`
        run(cl_cmd)
    end
    log_message("imported cell lines")

    #import the outcome, use AUC for the time being
    outcome_file = joinpath(PROJECT_ROOT, "data", "iorio", "json_AUCvalues.txt")
    outcome_add_cmd = `curl -X POST http://$host:$port/experiments/iorio/data/outcomes -d @$outcome_file`
    run(outcome_add_cmd)
    log_message("imported outcomes")

    # add pathway info
    pw_file = joinpath(PROJECT_ROOT, "data", "json", "pathways_current.json")
    run(`curl -X POST http://$host:$port/experiments/iorio/data/pathways -d @$pw_file`)
    log_message("imported pathway info")
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Created Iorio data set", "experiment_id" => "iorio")))
    return req.response
end


# function run_dreamchallenge_data(pathways_file::AbstractString, dest_dir::AbstractString="results/")
#     @async start_server(8888)
#     srv_shutdown_cmd = `curl http://localhost:8888/stop`
#     experiment = load_dream_challenge_data()
#
#     # stop_cmd = `curl http://localhost:8888/stop`
#     # run(stop_cmd)
#
#     # do the different data settings
#     # 1) full exome sequencing data
#     #   a) full data, no pathways
#     #   b) variance filtered, no pathways
#     #   c) full data, full pathways
#     #   d) variance filtered, full pathways
#     #   e) full data, subsumed pathways
#     #   f) variance filtered, subsumed pathways
#     # 2) all the above but with exome seq data constant 1
#
#     # 1)
#     do_experiment_stack(joinpath(dest_dir, "full_exome"), experiment, pathways_file)
#
#     # normalize again with constant exome data
#     experiment.is_normalized = false
#     # overwrite function
#     get_measurement_value(d::ExomeSeq) = 1
#
#     # 2)
#     do_experiment_stack(joinpath(dest_dir, "const_exome"), experiment, pathways_file)
#
#     @async run(srv_shutdown_cmd)
#     nothing
# end
#
# function do_experiment_stack(dest_dir::String, experiment::Experiment, pathways_file::String)
#     add_pathways_cmd = `curl -X POST http://localhost:8888/experiments/dream_challenge/pathways -d @$pathways_file`
#     ic = InferenceConfiguration()
#     ic.do_gridsearch = true
#     ic.compute_wpc_index = true
#     # a
#     ic.do_variance_filtering = false
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "full_data_no_pathways")
#     run_model(dep, inference_config = ic)
#
#     run(add_pathways_cmd)
#     # c
#     ic.subsume_pathways = false
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "full_data_full_pathways")
#     run_model(dep, inference_config = ic)
#
#     # e
#     ic.subsume_pathways = true
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "full_data_sub_pathways")
#     run_model(dep, inference_config = ic)
#
#     # remove pathways for experiment
#     empty!(experiment.pathway_information)
#
#     # b
#     ic.do_variance_filtering = true
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "var_filter_data_no_pathways")
#     run_model(dep, inference_config = ic)
#
#     # add pathways again
#     run(add_pathways_cmd)
#     # d
#     ic.subsume_pathways = false
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "var_filter_data_full_pathways")
#     run_model(dep, inference_config = ic)
#
#     # f
#     ic.subsume_pathways = true
#     dep = create_drug_efficacy_predictor(experiment, ic)
#     ic.target_dir = joinpath(dest_dir, "var_filter_data_sub_pathways")
#     run_model(dep, inference_config = ic)
#
# end
