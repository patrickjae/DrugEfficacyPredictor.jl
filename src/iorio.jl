function load_iorio_data()
    @async start_server(8888)
    sleep(4)

    base_dir = joinpath(PROJECT_ROOT, "data", "iorio")

    run(`curl -X POST http://localhost:8888/experiments -d '{"experiment_id":"iorio"}'`)
    @info "check for key" check = haskey(experiments_dictionary, "iorio")
    sleep(2)
    #import the cell lines
    for f in readdir(joinpath(base_dir, "cell_lines"))
        if !startswith(basename(f), "json") continue end
        filename = joinpath(base_dir, "cell_lines", f)
        cl_cmd = `curl -X POST http://localhost:8888/experiments/iorio/cell_lines -d @$filename`
        run(cl_cmd)
    end

    @info "imported cell lines"
    sleep(2)
    #import the outcome, use AUC for the time being
    outcome_file = joinpath(PROJECT_ROOT, "data", "iorio", "json_AUCvalues.txt")
    outcome_add_cmd = `curl -X POST http://localhost:8888/experiments/iorio/outcomes -d @$outcome_file`
    run(outcome_add_cmd)
    @info "imported outcomes"
    sleep(2)

    # add pathway info
    pw_file = joinpath(PROJECT_ROOT, "data", "json", "pathway_current.json")
    run(`curl -X POST http://localhost:8888/experiments/iorio/pathways -d @$pw_file`)
    @info "imported pathway info"
    sleep(2)

    # stop the server
    run(`curl http://localhost:8888/stop`)

    experiments_dictionary["iorio"]
end