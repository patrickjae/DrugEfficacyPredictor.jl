module Webservice
using HTTP, JSON, UUIDs, Distributed, Sockets
using ..Data, ..Models
import ..Utils: log_message, log_progress, PROJECT_ROOT

include("utils.jl")
include("request_handlers.jl")
include("sample_pipeline.jl")
include("server.jl")
include("routes.jl")

export start_server, stop_server
end
