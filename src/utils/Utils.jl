module Utils
using Distributed

function log_message(msg)
    prefix = ""
    if myid() == 1
        prefix = "      From master 1:    "
    end
    Core.println(prefix*msg)
end

function log_progress(experiment_id::String, s::String)
    @spawnat 1 push!(training_progress[experiment_id], s)
end
const PROJECT_ROOT = joinpath(@__DIR__, "..", "..")

export log_message, log_progress, PROJECT_ROOT

end
