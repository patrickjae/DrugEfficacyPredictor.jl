# server.jl
using HttpServer, JSON


# create a global variable to keep track of experiments in the server
experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()

# create an experiment for testing purposes
experiments_dictionary["test_id"] = create_experiment()
function create_http_handler()
http = HttpHandler() do req::Request, res::Response
    response = handle_request(req)
    return Response(response)
	end
	http
end

function handle_request(req::Request)
	display(req)
	func = eval(parse(replace(req.resource, "/", "DrugEfficacyPredictor.", 1)))
	display(req.data)

	println("test before")
	transcode(String, req.data)
	println("test after")
	JSON.parse(transcode(String, req.data))
	println("test after after")
    request_dictionary = JSON.parse(transcode(String, req.data))
    display(transcode(String, req.data))
    if endswith(req.resource, "create_experiment")	    
    	experiment = func()
    	exp_id = string(Base.Random.uuid4())
    	experiments_dictionary[exp_id] = experiment
		return JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
    elseif !haskey(request_dictionary, "experiment_id")
    	return JSON.json(Dict("status" => "error", "message" => "Please provide an experiment_id."))
    end
    exp_id = request_dictionary["experiment_id"]
    experiment = experiments_dictionary[exp_id]
    try
    	func(experiment, request_dictionary)
		return JSON.json(Dict("status" => "success", "experiment_id" => exp_id))
    catch e
		return JSON.json(Dict("status" => "exception", "type" => typeof(e), "stacktrace" => map(string, stacktrace(catch_backtrace()))))
	finally
    	rethrow()
	end
end

function start_server(port::Int64)
	server = Server(create_http_handler())
	run(server, port)
end