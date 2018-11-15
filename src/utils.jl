function run_dreamchallenge_data(dc_dir::AbstractString, pathways_file::AbstractString, dest_dir::AbstractString="results/")
    experiment = import_dream_challenge_data(dc_dir)
    @async start_server(8888)

    # stop_cmd = `curl http://localhost:8888/stop`
    # run(stop_cmd)

    # do the different data settings
    # 1) full exome sequencing data
    #   a) full data, no pathways
    #   b) variance filtered, no pathways
    #   c) full data, full pathways
    #   d) variance filtered, full pathways
    #   e) full data, subsumed pathways
    #   f) variance filtered, subsumed pathways
    # 2) all the above but with exome seq data constant 1

    # 1)
    do_experiment_stack(joinpath(dest_dir, "full_exome"), experiment)

    # normalize again with constant exome data
    experiment.is_normalized = false
    # overwrite function
    get_measurement_value(d::ExomeSeq) = 1
    
    do_experiment_stack(joinpath(dest_dir, "const_exome"), experiment)


    nothing
end

function do_experiment_stack(dest_dir::String, experiment::Experiment)
    add_pathways_cmd = `curl -X POST http://localhost:8888/experiments/dream_challenge/pathways -d @$pathways_file`
    # a
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = false)
    gridsearch(dep, joinpath(dest_dir, "full_data_no_pathways"))

    run(add_pathways_cmd)
    # c
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = false, subsume_pathways = false)
    gridsearch(dep, joinpath(dest_dir, "full_data_full_pathways"))

    # e
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = false, subsume_pathways = true)
    gridsearch(dep, joinpath(dest_dir, "full_data_sub_pathways"))

    # remove pathways for experiment
    empty!(experiment.pathway_information)

    # b
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = true)
    gridsearch(dep, joinpath(dest_dir, "var_filter_data_no_pathways"))

    # add pathways again
    run(add_pathways_cmd)
    # d
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = true, subsume_pathways = false)
    gridsearch(dep, joinpath(dest_dir, "var_filter_data_full_pathways"))

    # f
    dep = create_drug_efficacy_predictor(experiment, do_variance_filtering = true, subsume_pathways = true)
    gridsearch(dep, joinpath(dest_dir, "var_filter_data_sub_pathways"))

end