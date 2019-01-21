module DrugEfficacyPredictor

using CSV, DataFrames, SpecialFunctions, Printf,
	Statistics, LinearAlgebra, Distributed

include("distributional_parameters.jl")
include("types.jl")
include("../init.jl")
include("data_handling.jl")
include("dream_challenge_data.jl")
include("iorio.jl")
include("webservice.jl")

include("create.jl")
include("get.jl")
include("delete.jl")
include("model_initialization.jl")
include("inference.jl")
include("prediction.jl")
include("pipeline.jl")
include("utils.jl")

########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training



end # module
