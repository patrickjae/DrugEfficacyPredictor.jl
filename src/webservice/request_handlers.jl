function debug_request(req::HTTP.Request)
    Core.println(req)
    req.response.status = 200
    req.response
end

# base requests
function handle_base_request(req::HTTP.Request)
   target = req.target
   if endswith(target, "stop")
       log_message("stopping server")
       req.response.status = 200
       req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Server has been shut down.")))
       stop_server()
       return req.response
   end

   req.response.status = 404
   req.response.body = create_response(JSON.json(Dict("status" => "error", "message" => "Unknown URI, use a registered function to proceed.")))
   req.response
end

function no_experiment_response(req::HTTP.Request)
    experiment_id = get_experiment_id(req.target)
    req.response.status = 400
    req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "Experiment with ID '$experiment_id' does not exist.")))
    req.response
end

# computation requests

function train_request(req::HTTP.Request)
   request_dictionary = JSON.parse(transcode(String, req.body))
   experiment_id = get_experiment_id(req.target)

   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   proc_id = get_experiment_process(experiment_id)
   @spawnat proc_id Models.train(experiment_id, request_dictionary)

   req.response.status = 200
   req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Training in progress, check via /experiments/$experiment_id/progress. Once results are ready for download, access them via /experiments/$experiment_id/results.")))
   req.response
end

function progress_request(req::HTTP.Request)
   experiment_id = get_experiment_id(req.target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   req.response.status = 200
   req.response.body = create_response(JSON.json(Dict("status" => "success", "experiment_id" => experiment_id, "progress" => Utils.get_progress(experiment_id))))
   req.response
end

function results_request(req::HTTP.Request)
   experiment_id = get_experiment_id(req.target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   result_file = Utils.get_result(experiment_id)
   if get_result(experiment_id) == ""
       req.response.status = 500
       req.response.body = create_response(JSON.json(Dict("status" => "failure", "message" => "No results found for experiment with ID '$experiment_id'")))
       return req.response
   end
   req.response.status = 200
   req.response.headers = [Pair("Content-type", "application/zip")]
   req.response.body = read(result_file)
   req.response
end

function predict_request(req::HTTP.Request)
   request_dictionary = JSON.parse(transcode(String, req.body))
   experiment_id = get_experiment_id(req.target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   proc_id = Utils.get_experiment_process(experiment_id)

   prediction_result = nothing
   try
       prediction_result = remotecall_fetch(Models.predict, proc_id, experiment_id, request_dictionary)
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
       if !Utils.experiment_registered(exp_id_wanted)
           exp_id = exp_id_wanted
       else
           create_message *= " Your requested experiment_id is already in use, we have created a random one."
       end
   end
   proc_id = Utils.register_experiment(exp_id)
   remotecall_wait(Data.create_experiment, proc_id, exp_id)
   req.response.status = 200
   req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => create_message, "experiment_id" => exp_id)))
   req.response
end

function handle_get_experiment_request(req::HTTP.Request)
   experiment_id = get_experiment_id(req.target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   proc_id = get_experiment_process(experiment_id)
   response_string = remotecall_fetch(Data.to_json, proc_id, experiment_id)
   req.response.status = 200
   req.response.body = create_response(JSON.json(response_string))
   return req.response
end

function handle_delete_experiment_request(req::HTTP.Request)
   experiment_id = get_experiment_id(req.target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   proc_id = unregister_experiment(experiment_id)
   remotecall_wait(Data.delete_experiment, proc_id, experiment_id)
   req.response.status = 200
   req.response.body = create_response(JSON.json(Dict("status" => "success", "message" => "Deleted experiment with id '$experiment_id'")))
   return req.response
end


handle_create_collection_request(req::HTTP.Request) = handle_collection_request(req, "create")
handle_get_collection_request(req::HTTP.Request) = handle_collection_request(req, "get")
handle_delete_collection_request(req::HTTP.Request) = handle_collection_request(req, "delete")

function handle_collection_request(req::HTTP.Request, action::String)
   target = req.target
   experiment_id = get_experiment_id(target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end
   proc_id = get_experiment_process(experiment_id)
   data_target = get_function_target(target)

   target_function = eval(Meta.parse("Data.$(action)_$data_target"))
   try
       obj = nothing
       if action == "create"
           request_dictionary = JSON.parse(transcode(String, req.body))
           remotecall_wait(target_function, proc_id, experiment_id, request_dictionary)
       else
           obj = remotecall_fetch(target_function, proc_id, experiment_id)
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

   return req.response
end

handle_create_object_request(req::HTTP.Request) = handle_object_request(req, "create")
handle_get_object_request(req::HTTP.Request) = handle_object_request(req, "get")
handle_delete_object_request(req::HTTP.Request) = handle_object_request(req, "delete")

function handle_object_request(req::HTTP.Request, action::String)
   target = req.target
   experiment_id = get_experiment_id(target)
   if !Utils.experiment_registered(experiment_id)
       return no_experiment_response(req)
   end

   proc_id = get_experiment_process(experiment_id)
   get_target = chop(get_function_target(target))
   get_target_id = get_object_id(target)

   #use integer id if possible
   try get_target_id = parse(Int64, get_target_id) catch end

   #get the object
   getter_function = eval(Meta.parse("Data.$(action)_$get_target"))
   try
       obj = remotecall_fetch(getter_function, proc_id, experiment_id, get_target_id)
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
