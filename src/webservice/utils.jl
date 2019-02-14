# create a global variable to keep track of experiments and their processes
const experiments_dictionary = Dict{String, Int64}()
# keep track of the number of experiments for each process
const experiment_counts = Vector{Int64}(undef, nworkers())
# create a training progress message channel
const training_progress = Dict{String, Vector{String}}()

const result_file_dictionary = Dict{String, String}()


function get_experiment_id(uri::String)
    experiment_id = ""
    try
        experiment_id = HTTP.URIs.splitpath(uri)[2]
    catch
        throw(ArgumentError("Could not determine experiment ID in URI '$uri'"))
    end
    log_message("extracted experiment id $experiment_id")
    experiment_id
end

function add_experiment(experiment_id::String)
    # find a suitable worker (that with the minimal number of experiments)
    proc_idx = argmin(experiment_counts)
    proc_id = workers()[proc_idx]
    # remotecall_wait()
    @spawnat proc_id Data.create_experiment(experiment_id)
    Core.println("Webservice.add_experiment($experiment_id)")
    experiments_dictionary[experiment_id] = proc_id
    # update counts
    experiment_counts[proc_idx] += 1
    # init message queue
    training_progress[experiment_id] = Vector{String}()
end

function delete_experiment(experiment_id::String)
    proc_id = experiments_dictionary[experiment_id]
    @spawnat proc_id delete!(experiments, experiment_id)
    delete!(experiments_dictionary, experiment_id)
    delete!(training_progress, experiment_id)
end

get_function_target(uri::String) = HTTP.URIs.splitpath(uri)[4]
get_object_id(uri::String) = HTTP.URIs.splitpath(uri)[5]

create_response(s::String) = collect(codeunits(s))

experiment_exists(experiment_id::String) = haskey(experiments_dictionary, experiment_id)
progress_exists(experiment_id::String) = haskey(training_progress, experiment_id)
result_exists(experiment_id::String) = haskey(result_file_dictionary, experiment_id)
