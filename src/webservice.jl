# server.jl
using HTTP, JSON, UUIDs


# create an experiment for testing purposes
experiments_dictionary["test_id"] = DrugEfficacyPredictor.create_experiment()

# debug_mode = false

server = nothing


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
    log_message("extracted experiment object $experiment_id")
    experiment_id
end

get_function_target(uri::String) = HTTP.URIs.splitpath(uri)[4]
get_object_id(uri::String) = HTTP.URIs.splitpath(uri)[5]

create_response(s::String) = collect(codeunits(s))


# base requests
function handle_base_request(req::HTTP.Request)
    target = req.target
    if endswith(target, "stop")
        log_message("stopping server")
        req.response.status = 200
        req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Server has been shut down.")))
        stop_server(server)
        return req.response
    end

    req.response.status = 404
    req.response.body = create_response(JSON.json(Dict("status" => "error", "message" => "Unknown URI, use a registered function to proceed.")))
    req.response
end

# computation requests

function train_request(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    experiment_id = get_experiment_id(req.target)
    if !haskey(experiments_dictionary, experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "Experiment with ID '$experiment_id' does not exist.")))
        return req.response
    end
    experiment = experiments_dictionary[experiment_id]
    log_message("starting train method")
    @async begin
        try
            train(experiment, request_dictionary)
        catch ex
            st = map(string, stacktrace(catch_backtrace()))
            log_message("caught exception $ex")
            log_message("stacktrace")
            [log_message(ste) for ste in st]
        end
    end
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Training in progress, check via /experiments/$experiment_id/progress. Once results are ready for download, access them via /experiments/$experiment_id/results.")))
    req.response
end

function progress_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    if !haskey(experiments_dictionary, experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "Experiment with ID '$experiment_id' does not exist.")))
        return req.response
    end
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id, "progress" => training_progress[experiment_id])))
    req.response
end

function results_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    if !haskey(result_file_dictionary, experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "No results found for experiment with ID '$experiment_id'")))
        return req.response
    end
    req.response.status = 200
    req.response.headers = [Pair("Content-type", "application/zip")]
    req.response.body = read(result_file_dictionary[experiment_id])
    req.response
end

function predict_request(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    experiment_id = get_experiment_id(req.target)
    if !haskey(experiments_dictionary, experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "No results found for experiment with ID '$experiment_id'")))
        return req.response
    end
    prediction_result = nothing
    try
        prediction_result = predict(experiments_dictionary[experiment_id], request_dictionary)
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "results" => prediction_result)))
    return req.response
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
    log_message("creating experiment $exp_id")
    experiment.internal_id = exp_id
    experiments_dictionary[exp_id] = experiment
    training_progress[exp_id] = Vector{String}()
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => create_message, "experiment_id" => exp_id)))
    req.response
end

function handle_get_experiment_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    experiment = experiments_dictionary[experiment_id]
    req.response.status = 200
    req.response.body = create_response(JSON.json(to_json(experiment)))
    req.response
end

function handle_delete_experiment_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    message = "Deleted experiment with id '$experiment_id'"
    status = "success"
    if haskey(experiments_dictionary, experiment_id)
        delete!(experiments_dictionary, experiment_id)
    else
        message = "Experiment with id '$experiment_id' does not exist."
        status = "failure"
    end
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => status, "message" => message)))
    req.response
end


handle_create_collection_request(req::HTTP.Request) = handle_collection_request(req, "create")
handle_get_collection_request(req::HTTP.Request) = handle_collection_request(req, "get")
handle_delete_collection_request(req::HTTP.Request) = handle_collection_request(req, "delete")

function handle_collection_request(req::HTTP.Request, action::String)
    target = req.target
    experiment_id = get_experiment_id(target)
    data_target = get_function_target(target)

    target_function = eval(Meta.parse("$(action)_$data_target"))

    try  catch end

    try
        obj = nothing
        if action == "create"
            request_dictionary = JSON.parse(transcode(String, req.body))
            target_function(experiments_dictionary[experiment_id], request_dictionary)
        else
            obj = target_function(experiments_dictionary[experiment_id])
        end
        req.response.status = 200
        if obj == nothing
            # either create or delete
            req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id)))
        else
            # get
            req.response.body = create_response(JSON.json(Dict("status" => "success", "data" => to_json(obj))))
        end
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end

    req.response
end

handle_create_object_request(req::HTTP.Request) = handle_object_request(req, "create")
handle_get_object_request(req::HTTP.Request) = handle_object_request(req, "get")
handle_delete_object_request(req::HTTP.Request) = handle_object_request(req, "delete")

function handle_object_request(req::HTTP.Request, action::String)
    target = req.target
    experiment_id = get_experiment_id(target)
    get_target = chop(get_function_target(target))
    get_target_id = get_object_id(target)

    #use integer id if possible
    try get_target_id = parse(Int64, get_target_id) catch end

    #get the object
    getter_function = eval(Meta.parse("$(action)_$get_target"))
    try
        obj = getter_function(experiments_dictionary[experiment_id], get_target_id)
        req.response.status = 200
        if obj == nothing
            if action == "get"
                req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "$(uppercasefirst(get_target)) with id '$get_target_id' not found in experiment with id '$experiment_id'")))
            else
                req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id)))
            end
        else
            if action == "get"
                req.response.body = create_response(JSON.json(Dict("status" => "success", "data" => to_json(obj))))
            else
                req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id)))
            end
        end
    catch ex
        st = map(string, stacktrace(catch_backtrace()))
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "exception", "type" => string(ex), "stacktrace" => st)))
        return req.response
    end
    req.response
end


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

    # experiment methods§
    # create experiments
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments"), HTTP.Handlers.HandlerFunction(handle_create_experiment_request))

    # query for whole experiment
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*"), HTTP.Handlers.HandlerFunction(handle_get_experiment_request))

    # delete experiment
    HTTP.register!(r, "DELETE", HTTP.URI(path="/experiments/*"), HTTP.Handlers.HandlerFunction(handle_delete_experiment_request))

    # data methods
    # create collections, can be genes, proteins, cell_lines, drugs, pathways and outcomes
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments/*/data/*"), HTTP.Handlers.HandlerFunction(handle_create_collection_request))

    # query for collections
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/data/*"), HTTP.Handlers.HandlerFunction(handle_get_collection_request))

    # query single objects
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/data/*/*"), HTTP.Handlers.HandlerFunction(handle_get_object_request))

    # delete collections
    HTTP.register!(r, "DELETE", HTTP.URI(path="/experiments/*/data/*"), HTTP.Handlers.HandlerFunction(handle_delete_collection_request))

    # delete single objects
    HTTP.register!(r, "DELETE", HTTP.URI(path="/experiments/*/data/*/*"), HTTP.Handlers.HandlerFunction(handle_delete_object_request))

    # computation methods
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments/*/train"), HTTP.Handlers.HandlerFunction(train_request))
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/progress"), HTTP.Handlers.HandlerFunction(progress_request))
    HTTP.register!(r, "GET", HTTP.URI(path="/experiments/*/results"), HTTP.Handlers.HandlerFunction(results_request))
    HTTP.register!(r, "POST", HTTP.URI(path="/experiments/*/predict"), HTTP.Handlers.HandlerFunction(predict_request))

    r
end

stop_server(server::HTTP.Server) = put!(server.in, HTTP.Servers.KILL)


function start_server(port::Int64)
	# global server = HTTP.Server(create_http_handler())
    global server = HTTP.Server(create_base_router())
	try
		HTTP.serve(server, "127.0.0.1", port, verbose = true)
	finally
		stop_server(server)
	end
end
