using Distributed

sp = workers()[1]
dp = myid()

if nprocs() == nworkers()
    # no workers have ben started yet, start a process for the server
    println("adding server process")
    global sp = addprocs(1)[1]
end

using DrugEfficacyPredictor
# # define module on each process
@everywhere using DrugEfficacyPredictor

#broadcast server process id and data management process id
@everywhere DrugEfficacyPredictor.set_process_ids($sp, $dp)

# TODO make port configurable
@spawnat sp start_server(8888)
