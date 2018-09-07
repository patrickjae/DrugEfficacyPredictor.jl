module DrugEfficacyPredictor
using CSV, DataFrames, DataStructures, Distributions, SpecialFunctions, Printf,
	Statistics, LinearAlgebra

include("distributional_parameters.jl")
include("types.jl")
include("data_handling.jl")
include("dream_challenge_data.jl")
include("webservice.jl")

include("import_json.jl")
include("model_initialization.jl")
include("inference.jl")
include("prediction.jl")


########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training



end # module
