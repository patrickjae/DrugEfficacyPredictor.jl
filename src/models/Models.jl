module Models
#computations
using Statistics, LinearAlgebra, SpecialFunctions, Distributions, Printf, Distributed
using ..Utils, ..Data
# Core.println("loading submodule Models on process $(myid()): $(@__MODULE__)")

include("distributional_parameters.jl")
include("types.jl")
include("bmtmkl.jl")
include("run_models.jl")
include("io.jl")
include("utils.jl")
include("pipeline.jl")


end
