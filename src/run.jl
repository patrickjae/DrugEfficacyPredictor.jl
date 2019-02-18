# set the number of worker threads for managing experiments
num_worker_procs = length(ARGS) == 0 ? 2 : parse(Int64, ARGS[1])
# setting the number of threads to number of CPUs
# this is default behavior, setting it anyway to override if set otherwise
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS
# start the workers on which experiments are loaded
using Distributed
addprocs(num_worker_procs, topology=:master_worker)
w = workers()

#make sure we are able to load everything
@everywhere begin
    Core.println("loading code")
    import Pkg; Pkg.activate(".")
    using DrugEfficacyPredictor
    Core.println("done loading code")
end

# make sure workers now their ids
@everywhere DrugEfficacyPredictor.Utils.init_workers($w)

# start the webservice
DrugEfficacyPredictor.Webservice.start_server(8888)

# TODO make port configurable
