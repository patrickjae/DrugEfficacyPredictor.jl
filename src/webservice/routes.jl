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
    HTTP.@register(r, "/", HTTP.Handlers.HandlerFunction(handle_base_request))
	# load sample data
	HTTP.@register(r, "/load_iorio_data", load_iorio_data)
    # experiment methods
    # create experiments
    HTTP.@register(r, "POST", "/experiments", handle_create_experiment_request)

    # query for whole experiment
    HTTP.@register(r, "GET", "/experiments/*", handle_get_experiment_request)

    # delete experiment
    HTTP.@register(r, "DELETE", "/experiments/*", handle_delete_experiment_request)

    # data methods
    # create collections, can be genes, proteins, cell_lines, drugs, pathways and outcomes
    HTTP.@register(r, "POST", "/experiments/*/data/*", handle_create_collection_request)

    # query for collections
    HTTP.@register(r, "GET", "/experiments/*/data/*", handle_get_collection_request)

    # query single objects
    HTTP.@register(r, "GET", "/experiments/*/data/*/*", handle_get_object_request)

    # delete collections
    HTTP.@register(r, "DELETE", "/experiments/*/data/*", handle_delete_collection_request)

    # delete single objects
    HTTP.@register(r, "DELETE", "/experiments/*/data/*/*", handle_delete_object_request)

    # computation methods
    HTTP.@register(r, "POST", "/experiments/*/train", train_request)
    HTTP.@register(r, "GET", "/experiments/*/progress", progress_request)
    HTTP.@register(r, "GET", "/experiments/*/results", results_request)
    HTTP.@register(r, "POST", "/experiments/*/predict", predict_request)
    r
end
