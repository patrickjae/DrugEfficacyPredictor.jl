module Utils
using Distributed

# project root
const PROJECT_ROOT = joinpath(@__DIR__, "..", "..")
# create a global variable to keep track of experiments and their processes
const experiments_dictionary = Dict{String, Int64}()
# keep track of the number of experiments for each process
experiment_counts = zeros(Int64, nworkers())
worker_ids = workers()
# create a training progress message channel
const training_progress = Dict{String, Vector{String}}()
# storage for result file names
const result_file_dictionary = Dict{String, String}()

function log_message(msg)
    prefix = ""
    if myid() == 1
        prefix = "      From master 1:    "
    end
    Core.println(prefix*msg)
end

function init_workers(workers::Vector{Int64})
    global worker_ids = workers
    global experiment_counts = zeros(Int64, length(workers))
end

experiment_registered(experiment_id::String) = haskey(experiments_dictionary, experiment_id)

get_experiment_process(experiment_id::String) = experiments_dictionary[experiment_id]

log_progress(experiment_id::String, s::String) = @spawnat 1 push!(training_progress[experiment_id], s)

get_progress(experiment_id::String) = training_progress[experiment_id]

store_result(experiment_id::String, result_file::String) = @spawnat 1 push!(result_file_dictionary, experiment_id => result_file)

get_result(experiment_id::String) = result_file_dictionary[experiment_id]

function register_experiment(experiment_id::String)
    # find a suitable worker (that with the minimal number of experiments)
    proc_idx = argmin(experiment_counts)
    proc_id = worker_ids[proc_idx]
    experiments_dictionary[experiment_id] = proc_id
    # update counts
    experiment_counts[proc_idx] += 1
    # init message queue
    training_progress[experiment_id] = Vector{String}()
    result_file_dictionary[experiment_id] = ""
    proc_id
end

function unregister_experiment(experiment_id::String)
    proc_id = get_experiment_process(experiment_id)
    delete!(experiments_dictionary, experiment_id)
    delete!(training_progress, experiment_id)
    result_file = get_result(experiment_id)
    if result_file != ""
        rm(result_file, force = true)
    end
    delete!(result_file_dictionary, experiment_id)
    proc_id
end

export experiment_registered, register_experiment, unregister_experiment
export log_message, log_progress, get_experiment_process, get_progress, get_result
export PROJECT_ROOT

end
