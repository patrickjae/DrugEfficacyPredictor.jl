using DrugEfficacyPredictor

const PROJECT_ROOT = @__DIR__

# create a global variable to keep track of experiments in the server
const experiments_dictionary = Dict{String, DrugEfficacyPredictor.Experiment}()

const predictor_dictionary = Dict{String, DrugEfficacyPredictor.DrugEfficacyPrediction}()

# create a job channel
# const training_jobs = Channel{}

# create a training progress message channel
const training_progress = Dict{String, Vector{String}}()

const result_file_dictionary = Dict{String, String}()

log_message(msg) = Core.println("Thread $(Threads.threadid()): $msg")
