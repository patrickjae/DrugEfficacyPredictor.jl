using Distributed
addprocs(topology=:master_worker)

#make sure we are able to load everything
@everywhere begin
    import Pkg; Pkg.activate(".")
    using DrugEfficacyPredictor
end
DrugEfficacyPredictor.Webservice.start_server(8888)

# sp = workers()[1]
# dp = myid()
#
# if nprocs() == nworkers()
#     # no workers have ben started yet, start a process for the server
#     println("adding server process")
#     global sp = addprocs(1)[1]
# end
#
# using DrugEfficacyPredictor
# # # define module on each process
# @everywhere using DrugEfficacyPredictor
#
# #broadcast server process id and data management process id
# @everywhere DrugEfficacyPredictor.set_process_ids($sp, $dp)
#
# # TODO make port configurable
# @spawnat sp DrugEfficacyPredictor.start_server(8888)
