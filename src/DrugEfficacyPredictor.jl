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
include("prediction.jl")

"""
Adjust info and error functions to include time stamp.
"""

const info_messages = Channel{Tuple}(32)
const warn_messages = Channel{Tuple}(32)
const error_messages = Channel{Tuple}(32)

"""
Start the logging processes
"""
function start_info_logging()
	for msg in info_messages
		Base.info(msg..., prefix = "[$(now())]\t")
	end
end

function start_warn_logging()
	for msg in warn_messages
		Base.warn(msg..., prefix = "[$(now())]\t")
	end
end

function start_error_logging()
	for msg in error_messages
		try 
			Base.error("[$(now())]\t", msg...)
		catch exc
			println(exc)
		end
	end
end

function start_logging()
	@async start_info_logging()
	@async start_warn_logging()
	@async start_error_logging()
end

start_logging()

import Base.info
import Base.warn
import Base.error

info(message...) = put!(info_messages, message)
warn(message...) = put!(warn_messages, message)
error(message...) = put!(error_messages, message)






########################
######### TODO #########
########################

# steps

# populate model - OK

# expose API _ OK (partly)
# split application and API - OK (partly)

# cross-validation during training



end # module
