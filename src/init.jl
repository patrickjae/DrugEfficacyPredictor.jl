println("loading module")
module DrugEfficacyPredictor

println("using Distributed")
using Distributed
println("using Distributed end")

const PROJECT_ROOT = joinpath(@__DIR__, "..")
println("project root dir is $PROJECT_ROOT")
server_process = workers()[1]
data_process = procs()[1]
println("setting server process id to $server_process, data process id to $data_process")

log_message(msg) = Core.println("Thread $(Threads.threadid()): $msg")
println("defined log_message")
# helpers
using Printf
println("imported Printf")
using DataStructures
println("imported DataStructures")
#data
using CSV, DataFrames
println("imported CSV, DataFrames")


#computations
using Statistics, LinearAlgebra, SpecialFunctions
println("imported Statistics, LinearAlgebra, SpecialFunctions")

# web service
using HTTP, JSON, UUIDs
println("after using statements")


include("distributional_parameters.jl")
println("loaded distributional_parameters")
include("types.jl")
println("loaded types")
include("data_handling.jl")
println("loaded data_handling")
include("dream_challenge_data.jl")
println("loaded dream_challenge_data")
include("iorio.jl")
println("loaded iorio")
include("webservice.jl")
println("loaded webservice")

include("create.jl")
println("loaded create")
include("get.jl")
println("loaded get")
include("delete.jl")
println("loaded delete")
include("model_initialization.jl")
println("loaded model_initialization")
include("inference.jl")
println("loaded inference")
include("prediction.jl")
println("loaded prediction")
include("pipeline.jl")
println("loaded pipeline")
include("utils.jl")
println("loaded utils")

println("end of module def")
########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training



end # module
