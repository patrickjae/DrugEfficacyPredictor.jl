module DrugEfficacyPredictor
using Distributed

server_process = workers()[1]
data_process = procs()[1]

log_message(msg) = Core.println("Thread $(Threads.threadid()): $msg")
function set_process_ids(sp::Int64, dp::Int64)
    global server_process = sp
    global data_process = dp
    nothing
end
# helpers
using Printf

using DataStructures

#data
using CSV, DataFrames

#computations
using Statistics, LinearAlgebra, SpecialFunctions


const PROJECT_ROOT = joinpath(@__DIR__, "..")

include(joinpath(@__DIR__, "distributional_parameters.jl"))
include(joinpath(@__DIR__, "types.jl"))
include(joinpath(@__DIR__, "data_handling.jl"))
include(joinpath(@__DIR__, "dream_challenge_data.jl"))
include(joinpath(@__DIR__, "iorio.jl"))

include(joinpath(@__DIR__, "model_initialization.jl"))
include(joinpath(@__DIR__, "inference.jl"))
include(joinpath(@__DIR__, "prediction.jl"))
include(joinpath(@__DIR__, "pipeline.jl"))
include(joinpath(@__DIR__, "utils.jl"))


# web service
using HTTP, JSON, UUIDs

include(joinpath(@__DIR__, "webservice.jl"))
include(joinpath(@__DIR__, "create.jl"))
include(joinpath(@__DIR__, "get.jl"))
include(joinpath(@__DIR__, "delete.jl"))

########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training
end # module
