function run_dreamchallenge_data(dc_dir::AbstractString, pathways_file::AbstractString, dest_dir::AbstractString="results/")
    experiment = import_dream_challenge_data(dc_dir)
    @async start_server(8888)

    sleep(2)
    curl_cmd = `curl -X POST http://localhost:8888/experiments/dream_challenge/pathways -d @$pathways_file`
    run(curl_cmd)

    # stop_cmd = `curl http://localhost:8888/stop`
    # run(stop_cmd)

    dep = create_drug_efficacy_predictor(experiment)

    gridsearch(dep, dest_dir)

    dep
end