function load_iorio_data(host::String = "localhost", port::Int64 = 8888; force_reload::Bool = false)
    base_dir = joinpath(PROJECT_ROOT, "data", "iorio")

    if haskey(experiments_dictionary, "iorio")
        if force_reload
            empty!(experiments_dictionary)
            empty!(training_progress)
        else
            log_message("Iorio data set already exists, returning experiment object. If you need to reload, call with parameter force_reload = true")
            return experiments_dictionary["iorio"]
        end
    end

    run(`curl -X POST http://$host:$port/experiments -d '{"experiment_id":"iorio"}'`)
    sleep(2)
    #import the cell lines
    for f in readdir(joinpath(base_dir, "cell_lines"))
        if !startswith(basename(f), "json") continue end
        filename = joinpath(base_dir, "cell_lines", f)
        cl_cmd = `curl -X POST http://$host:$port/experiments/iorio/data/cell_line -d @$filename`
        run(cl_cmd)
    end

    log_message("imported cell lines")
    sleep(2)
    #import the outcome, use AUC for the time being
    outcome_file = joinpath(PROJECT_ROOT, "data", "iorio", "json_AUCvalues.txt")
    outcome_add_cmd = `curl -X POST http://$host:$port/experiments/iorio/data/outcomes -d @$outcome_file`
    run(outcome_add_cmd)
    log_message("imported outcomes")
    sleep(2)

    # add pathway info
    pw_file = joinpath(PROJECT_ROOT, "data", "json", "pathways_current.json")
    run(`curl -X POST http://$host:$port/experiments/iorio/data/pathways -d @$pw_file`)
    log_message("imported pathway info")
    sleep(2)

    # stop the server
    # run(`curl http://localhost:8888/stop`)

    experiments_dictionary["iorio"]
end
