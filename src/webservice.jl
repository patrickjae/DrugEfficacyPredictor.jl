# server.jl
using HTTP, JSON, UUIDs

# create a global variable to keep track of experiments in the server
experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()
predictor_dictionary = Dict{String, DrugEfficacyPredictor. DrugEfficacyPrediction}()

# create an experiment for testing purposes
experiments_dictionary["test_id"] = create_experiment()

# debug_mode = false

server = nothing

# base requests
function handle_base_request(req::HTTP.Request)
    target = req.target
    if endswith(target, "stop")
        @info "stopping server"
        req.response.status = 200
        req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Server has been shut down.")))
        stop_server(server)
        return req.response
    end

    # if endswith(target,"unset_debug_mode")
    #     @info "unsetting debug mode"
    #     req.response.status = 200
    #     req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Experiment ID will be set to 'test_id' permanently.")))
    #     global debug_mode = false
    #     return req.response
    # end

    # if endswith(target,"set_debug_mode")
    #     @info "setting debug mode"
    #     req.response.status = 200
    #     req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Experiment ID will be set to 'test_id' permanently.")))
    #     global debug_mode = true
    #     return req.response
    # end
    req.response.status = 404
    req.response.body = create_response(JSON.json(Dict("status" => "error", "message" => "Unknown URI, use a registered function to proceed.")))
    req.response
end

# create experiment
function handle_create_experiment_request(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    exp_id = string(UUIDs.uuid4())
    exp_id_wanted = ""
    create_message = "Created new experiment."
    if haskey(request_dictionary, "experiment_id")
        exp_id_wanted = request_dictionary["experiment_id"]
        if !haskey(experiments_dictionary, exp_id_wanted)
            exp_id = exp_id_wanted
        else
            create_message *= " Your requested experiment_id is already in use, we have created a random one."
        end
    end
    experiment = DrugEfficacyPredictor.create_experiment()
    @info "creating experiment" experiment_id = exp_id
    experiments_dictionary[exp_id] = experiment
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => create_message, "experiment_id" => exp_id)))
    req.response
end

function get_experiment_id(uri::String) 
    experiment_id = ""
    try 
        experiment_id = HTTP.URIs.splitpath(uri)[2]
    catch 
        throw(ArgumentError("Could not determine experiment ID in URI '$uri'"))
    end
    if !haskey(experiments_dictionary, experiment_id)
        throw(ArgumentError("Experiment ID '$experiment_id' not found, have you created it?"))
    end
    @info "extracted experiment object" experiment_id
    experiment_id
end

get_function_target(uri::String) = HTTP.URIs.splitpath(uri)[3]
get_object_id(uri::String) = HTTP.URIs.splitpath(uri)[4]

function handle_experiment_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    experiment = experiments_dictionary[experiment_id]
    req.response.status = 200
    req.response.body = create_response(JSON.json(to_json(experiment)))
    req.response
end

function handle_post_request(req::HTTP.Request)
    target = req.target
    experiment_id = get_experiment_id(target)
    experiment = experiments_dictionary[experiment_id]

    post_subject = get_function_target(target)
    # TODO check for whitelist of post_subject

    target_function = eval(Meta.parse("import_$post_subject"))

    @info "parsing request dictionary"
    request_dictionary = JSON.parse(transcode(String, req.body))

    @info "calling function" target_function
    try
        target_function(experiment, request_dictionary)
        req.response.status = 200
        req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id)))
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end

    req.response
end

function handle_get_collection_request(req::HTTP.Request)
    target = req.target
    experiment_id = get_experiment_id(target)
    get_target = get_function_target(target)

    #get the object
    getter_function = eval(Meta.parse("$(get_target)_to_json"))
    try 
        obj = getter_function(experiments_dictionary[experiment_id])
        req.response.status = 200
        req.response.body = create_response(JSON.json(Dict("status" => "success", "data" => to_json(obj))))
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end
    req.response
end

function handle_get_object_request(req::HTTP.Request)
    target = req.target
    experiment_id = get_experiment_id(target)
    get_target = chop(get_function_target(target))
    get_target_id = get_object_id(target)

    #use integer id if possible
    try get_target_id = parse(Int64, get_target_id) catch end

    #get the object
    getter_function = eval(Meta.parse("get_$get_target"))
    try 
        obj = getter_function(get_target_id)
        req.response.status = 200
        req.response.body = create_response(JSON.json(Dict("status" => "success", "data" => to_json(obj))))
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end
    req.response
end

# function handle_request(req::HTTP.Request)
#     target = req.target
# 	func = eval(Meta.parse(replace(target, "/" => "DrugEfficacyPredictor.", count=1)))

#     # for all requests, require the experiment id.
#     request_dictionary = JSON.parse(transcode(String, req.body))
#     if !haskey(request_dictionary, "experiment_id") && !debug_mode
#         req.response.status = 500
#     	req.response.body = create_response(JSON.json(Dict("status" => "error", "message" => "Please provide an experiment_id.")))
#         return req.response
#     end
#     exp_id = debug_mode ? "test_id" : request_dictionary["experiment_id"]
#     if !haskey(experiments_dictionary, exp_id)
#         req.response.status = 500
#     	req.response.body = create_response(JSON.json(Dict("status" => "error", "message" => "Unknown experiment ID '$exp_id'.")))
#         return req.response
#     end
#     experiment = experiments_dictionary[exp_id]

#     try
#     	func(experiment, request_dictionary)
#         req.response.status = 200
# 		req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => exp_id)))
#     catch ex
#     	st = map(string, stacktrace(catch_backtrace()))
#     	# error("Exception occurred: $(typeof(ex))")
#     	# @error "Stacktrace" st
#         rethrow(ex)
#         req.response.status = 500
# 		req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
#         return req.response
# 	end
#     req.response
# end

create_response(s::String) = Vector{UInt8}(s)

"""
Create the base routes for this webservice.
The basic interface is ```/experiments```.
Issuing a POST request to ```/experiments```, will create a new experiment and report back the experiment id.

For accessing previously created experiments, append the experiment id to the URL.
I.e. to access experiment1, aim your requests at ```/experiments/experiment1```.
A simple GET request to this URL with give back a full JSON representation of this experiment (carefull, this can be large).

Access to the different parts of data in an experiment is given by again appending to the URL.
Currently, the different data types that can be accessed or worked with are:
 - genes
 - proteins
 - cell_lines
 - drugs
 - pathways
 - outcomes

POST request to the according URLs will create the data inside of an experiment.
I.e. to import gene representations into experiment1, issue a POST request to ```/experiments/experiment1/genes``` and provide the
JSON representation of a list of genes.

"""
function create_base_router()
    r = HTTP.Handlers.Router()
    base_uri = HTTP.URI(path="/")
    HTTP.register!(r, base_uri, HTTP.Handlers.HandlerFunction(handle_base_request))
    # create experiments
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments"), HTTP.Handlers.HandlerFunction(handle_create_experiment_request))

    # create collections, can be genes, proteins, cell_lines, drugs, pathways and outcomes
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments/*/*"), HTTP.Handlers.HandlerFunction(handle_post_request))

    # query for whole experiment
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*"), HTTP.Handlers.HandlerFunction(handle_experiment_request))

    # query for collections
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/*"), HTTP.Handlers.HandlerFunction(handle_get_collection_request))

    # query single objects
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/*/*"), HTTP.Handlers.HandlerFunction(handle_get_object_request))

    # delete single objects
    # HTTP.register!(r, "DELETE", HTTP.URI(path="/experiments/*/*/*"), HTTP.Handlers.HandlerFunction(handle_delete_object_request))

    r
end 

stop_server(server::HTTP.Server) = put!(server.in, HTTP.Servers.KILL)


function start_server(port::Int64)
	# global server = HTTP.Server(create_http_handler())
    global server = HTTP.Server(create_base_router())
	try 
		HTTP.serve(server, "127.0.0.1", port)
	finally
		stop_server(server)
	end	
end