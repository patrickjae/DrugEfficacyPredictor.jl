module DrugEfficacyPredictor
using CSV, DataFrames, DataStructures, Distributions, SpecialFunctions
include("distributional_parameters.jl")
include("types.jl")
include("data_handling.jl")
include("dream_challenge_data.jl")
include("webservice.jl")

include("import_json.jl")
include("model_initialization.jl")
include("inference.jl")

"""
Adjust info and error functions to include time stamp.
"""
info(message::String) = Base.info("\t[$(now())] $message")
error(message::String) = Base.error("\t[$(now())] $message")
warn(message::String) = Base.warn("[$(now())] $message")





########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training



end # module
