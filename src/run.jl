using Distributed
addprocs(topology=:master_worker)

#make sure we are able to load everything
@everywhere begin
    import Pkg; Pkg.activate(".")
    using DrugEfficacyPredictor
end
DrugEfficacyPredictor.Webservice.start_server(8888)

# TODO make port configurable
