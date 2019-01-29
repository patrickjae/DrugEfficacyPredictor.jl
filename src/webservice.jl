# server.jl

# create a global variable to keep track of experiments in the server
const experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()

# create a training progress message channel
const training_progress = Dict{String, Vector{String}}()

const result_file_dictionary = Dict{String, String}()

# create an experiment for testing purposes
experiments_dictionary["test_id"] = DrugEfficacyPredictor.create_experiment("test_id")

# debug_mode = false

server = nothing

function get_experiment_id(uri::String)
    experiment_id = ""
    try
        experiment_id = HTTP.URIs.splitpath(uri)[2]
    catch
        throw(ArgumentError("Could not determine experiment ID in URI '$uri'"))
    end
    # if !experiment_exists(experiment_id)
    #     throw(ArgumentError("Experiment ID '$experiment_id' not found, have you created it?"))
    # end
    log_message("extracted experiment object $experiment_id")
    experiment_id
end

function set_experiment(id::String, experiment::Experiment)
    global experiments_dictionary[id] = experiment
    global training_progress[id] = Vector{String}()
end

get_function_target(uri::String) = HTTP.URIs.splitpath(uri)[4]
get_object_id(uri::String) = HTTP.URIs.splitpath(uri)[5]

create_response(s::String) = collect(codeunits(s))

init_training_progress(experiment_id::String) = remote_do((eid) -> training_progress[eid] = Vector{String}(), server_process, experiment_id)

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
experiment_exists(experiment_id::String) = remotecall_fetch((eid) -> haskey(experiments_dictionary, eid), data_process, experiment_id)
progress_exists(experiment_id::String) = remotecall_fetch((eid) -> haskey(training_progress, eid), server_process, experiment_id)
result_exists(experiment_id::String) = remotecall_fetch((eid) -> haskey(result_file_dictionary, eid), data_process, experiment_id)

function train_request(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    experiment_id = get_experiment_id(req.target)

    if !experiment_exists(experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "Experiment with ID '$experiment_id' does not exist.")))
        return req.response
    end
    # experiment = experiments_dictionary[experiment_id]
    log_message("starting train method")
    remote_do(train, data_process, experiment_id, request_dictionary)

    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Training in progress, check via /experiments/$experiment_id/progress. Once results are ready for download, access them via /experiments/$experiment_id/results.")))
    req.response
end

function progress_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    if !progress_exists(experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "Experiment with ID '$experiment_id' does not exist.")))
        return req.response
    end
    log_message("found experiment id in progress tracker")
    req.response.status = 200
    # progress = remotecall_fetch((eid) -> training_progress[eid], data_process, experiment_id)
    progress = @fetchfrom server_process training_progress[experiment_id]
    req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id, "progress" => progress)))
    req.response
end

function results_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    if !result_exists(experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "No results found for experiment with ID '$experiment_id'")))
        return req.response
    end
    result_file = remotecall_fetch((eid) -> result_file_dictionary[eid], data_process, experiment_id)
    req.response.status = 200
    req.response.headers = [Pair("Content-type", "application/zip")]
    req.response.body = read(result_file)
    req.response
end

function predict_request(req::HTTP.Request)
    request_dictionary = JSON.parse(transcode(String, req.body))
    experiment_id = get_experiment_id(req.target)
    if !experiment_exists(experiment_id)
        req.response.status = 500
        req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "No results found for experiment with ID '$experiment_id'")))
        return req.response
    end
    prediction_result = nothing
    try
        prediction_result = remotecall_fetch((eid, req_params) -> predict(experiments_dictionary[eid], req_params), data_process, experiment_id, request_dictionary)
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
        if !experiment_exists(exp_id_wanted)
            exp_id = exp_id_wanted
        else
            create_message *= " Your requested experiment_id is already in use, we have created a random one."
        end
    end
    experiment = remotecall_fetch(DrugEfficacyPredictor.create_experiment, data_process, exp_id)
    log_message("created experiment $(experiment.internal_id)")
    # experiment.internal_id = exp_id
    # remote_do((eid, exprmt) -> experiments_dictionary[eid] = exprmt, data_process, exp_id, experiment)
    # experiments_dictionary[exp_id] = experiment
    # remote_do((eid) -> training_progress[eid] = Vector{String}(), data_process, exp_id)
    init_training_progress(experiment_id)
    # training_progress[experiment.internal_id] = Vector{String}()
    req.response.status = 200
    req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => create_message, "experiment_id" => exp_id)))
    req.response
end

function handle_get_experiment_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    response_string = remotecall_fetch((eid) -> to_json(experiments_dictionary[eid]), data_process, experiment_id)
    # experiment = experiments_dictionary[experiment_id]
    req.response.status = 200
    # req.response.body = create_response(JSON.json(to_json(experiment)))
    req.response.body = create_response(JSON.json(response_string))
    req.response
end

function handle_delete_experiment_request(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    message = "Deleted experiment with id '$experiment_id'"
    status = "success"
    if exists(experiment_id)
        remote_do((eid) -> delete!(experiments_dictionary, eid), data_process, experiment_id)
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
    try
        obj = nothing
        if action == "create"
            request_dictionary = JSON.parse(transcode(String, req.body))
            remote_do(target_function, data_process, experiment_id, request_dictionary)
            # target_function(experiments_dictionary[experiment_id], request_dictionary)
        else
            obj = remotecall_fetch(target_function, data_process, experiment_id)
            # obj = target_function(experiments_dictionary[experiment_id])
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
        obj = remotecall_fetch((eid, tid) -> getter_function(experiments_dictionary[eid], tid), data_process, experiment_id, get_target_id)
        # obj = getter_function(experiments_dictionary[experiment_id], get_target_id)
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
    println("start create router")
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
    println("finish create router")
    r
end

stop_server(server::HTTP.Server) = @spawnat server_process put!(DrugEfficacyPredictor.server.in, HTTP.Servers.KILL)

function start_server(port::Int64)
    global server = HTTP.Server(DrugEfficacyPredictor.create_base_router())
    try
        HTTP.serve(DrugEfficacyPredictor.server, "127.0.0.1", port, verbose = true)
    finally
        stop_server(DrugEfficacyPredictor.server)
	end
end
