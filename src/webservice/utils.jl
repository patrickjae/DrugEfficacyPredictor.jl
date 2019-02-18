function get_experiment_id(uri::String)
    experiment_id = ""
    try
        experiment_id = HTTP.URIs.splitpath(uri)[2]
    catch
        throw(ArgumentError("Could not determine experiment ID in URI '$uri'"))
    end
    log_message("extracted experiment id $experiment_id")
    experiment_id
end

get_function_target(uri::String) = HTTP.URIs.splitpath(uri)[4]
get_object_id(uri::String) = HTTP.URIs.splitpath(uri)[5]

create_response(s::String) = collect(codeunits(s))
