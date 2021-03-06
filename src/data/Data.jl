module Data
using Statistics, DataStructures, CSV, DataFrames, Distributed
# Core.println("loading submodule Data on process $(myid()): $(@__MODULE__)")
import ..Utils: log_message, log_progress, PROJECT_ROOT
using ..Utils
include("types.jl")
# this data structure holds all experiments on this process
const experiments = Dict{String, Experiment}()
include("data_handling.jl")
include("create.jl")
include("get.jl")
include("delete.jl")
include("dream_challenge_data.jl")

end  # module Data
