# server.jl
using HttpServer, JSON

# create a global variable to keep track of experiments in the server
experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()

# create an experiment for testing purposes
experiments_dictionary["test_id"] = create_experiment()

server = nothing
function create_http_handler()
	http = HttpHandler() do req::Request, res::Response
    	response = handle_request(req)
    	return Response(response)
	end
	http
end

function handle_request(req::Request)
	func = eval(parse(replace(req.resource, "/", "DrugEfficacyPredictor.", 1)))
    request_dictionary = JSON.parse(transcode(String, req.data))
    # if request wants to create an experiment return the id of the newly created object
    if endswith(req.resource, "create_experiment")
    	experiment = func()
    	exp_id = string(Base.Random.uuid4())
    	experiments_dictionary[exp_id] = experiment
		return JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
    end

    if endswith(req.resource, "stop")
    	close(server)
    end
    # for all other requests, require the experiment id.
    if !haskey(request_dictionary, "experiment_id")
    	return JSON.json(Dict("status" => "error", "message" => "Please provide an experiment_id."))
    end
    exp_id = request_dictionary["experiment_id"]
    if !haskey(experiments_dictionary, exp_id)
    	return JSON.json(Dict("status" => "error", "message" => "Unknown experiment ID '$exp_id'."))
    end
    experiment = experiments_dictionary[exp_id]
    try
    	func(experiment, request_dictionary)
    	info("Request handled sucessfully, returning experiment $exp_id")
		return JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
    catch e
    	st = map(string, stacktrace(catch_backtrace()))
    	# error("Exception occurred: $(typeof(e))")
    	# error("Stacktrace:")
    	# for st_entry in st
    	# 	error(st)
    	# end
    	# error(string(e))
		return JSON.json(Dict("status" => "exception", "type" => string(e), "stacktrace" => st))
	finally
		rethrow(e)	
	end
end

function start_server(port::Int64)
	global server = Server(create_http_handler())
	try 
		run(server, port)
	finally
		close(server)
	end	
end