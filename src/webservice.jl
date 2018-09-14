# server.jl
using HTTP, JSON, Random

# create a global variable to keep track of experiments in the server
experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()

# create an experiment for testing purposes
experiments_dictionary["test_id"] = create_experiment()


function handle_request(req::HTTP.Request)
    target = req.target
    if endswith(target, "stop")
        req.response.status = 200
        req.response.body = JSON.json(Dict("status" => "success", "message" => "Server has been shut down."))
        close(server)
        return req.response
    end

	func = eval(Meta.parse(replace(target, "/" => "DrugEfficacyPredictor.", count=1)))
    # if request wants to create an experiment return the id of the newly created object
    if endswith(target, "create_experiment")
    	experiment = func()
    	exp_id = string(Random.uuid4())
    	experiments_dictionary[exp_id] = experiment
        req.response.status = 200
		req.response.body = JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
        return req.response
    end

    # for all other requests, require the experiment id.
    request_dictionary = JSON.parse(transcode(String, req.body))
    if !haskey(request_dictionary, "experiment_id")
        req.response.status = 500
    	req.response.body = JSON.json(Dict("status" => "error", "message" => "Please provide an experiment_id."))
        return req.response
    end
    exp_id = request_dictionary["experiment_id"]
    if !haskey(experiments_dictionary, exp_id)
        req.response.status = 500
    	req.response.body = JSON.json(Dict("status" => "error", "message" => "Unknown experiment ID '$exp_id'."))
        return req.response
    end
    experiment = experiments_dictionary[exp_id]

    try
    	func(experiment, request_dictionary)
        req.response.status = 200
		req.response.body = JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
    catch ex
    	st = map(string, stacktrace(catch_backtrace()))
    	# error("Exception occurred: $(typeof(ex))")
    	# @error "Stacktrace" st
        req.response.status = 500
		req.response.body = JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st))
        return req.response
	end
    req.response
end

server = nothing

create_http_handler() = HTTP.Handlers.HandlerFunction(handle_request)

import Base.close

close(server::HTTP.Server) = put!(server.in, HTTP.Servers.KILL)


function start_server(port::Int64)
	global server = HTTP.Server(create_http_handler())
	try 
		HTTP.serve(server, "127.0.0.1", port)
	finally
		close(server)
	end	
end