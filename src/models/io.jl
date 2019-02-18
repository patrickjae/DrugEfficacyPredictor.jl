function create_response_file(pm::PredictionModel, test::Bool, target_dir::String, filename::String)
    response_file = joinpath(target_dir, filename)
    f = open(response_file, "w")
    @printf(f, "CellLine")
    drugs = collect(keys(pm.data.results))
    for drug in drugs
        @printf(f, "\t%s", drug.id)
    end
    @printf(f, "\n")

    for cl in collect(values(pm.data.cell_lines))
        # if these two are equal
        if test == cl.in_test_set
            @printf(f, "%s", cl.id)
            for drug in drugs
                if haskey(pm.data.results[drug].outcome_values, cl)
                    @printf(f, "\t%.9f", pm.data.results[drug].outcome_values[cl])
                else
                    @printf(f, "\tNA")
                end
            end
            @printf(f, "\n")
        end
    end
    close(f)
    response_file
end


function write_results(pm::PredictionModel, p::PredictionModelParameters,
        parent_dir::String, filename::String,
        predictions::Dict{Drug, Vector{Float64}}, ranks::Dict{Drug, Vector{Int64}})
    # write predictions and true values
    mkpath(joinpath(parent_dir, "prediction"))
    prediction_file = joinpath(parent_dir, "prediction", filename)
    f = open(prediction_file, "w")
    @printf(f, "DrugAnonID\tType")
    for d in collect(keys(pm.data.results))
        @printf(f, "\t%s (prediction)\t%s (measured)", d.id, d.id)
    end
    @printf(f, "\n")
    for (cl_id, cl) in enumerate(collect(values(pm.data.cell_lines)))
        @printf(f, "%s", cl.id)
        if cl.in_test_set @printf(f, "\tTest") else @printf(f, "\tTraining") end
        for drug in collect(keys(pm.data.results))
            # drug = pm.data.drugs["Drug$d_id"]
            target = 0.
            if haskey(pm.data.results, drug) && haskey(pm.data.results[drug].outcome_values, cl)
                target = pm.data.results[drug].outcome_values[cl]
            # elseif haskey(pm.data.test_results, drug) && haskey(pm.data.test_results[drug].outcome_values, cl)
            #     target = pm.data.test_results[drug].outcome_values[cl]
            end
            @printf(f, "\t%.5f\t%.5f", predictions[drug][cl_id], target)
        end
        @printf(f, "\n")
    end
    close(f)
    # write rankings
    mkpath(joinpath(parent_dir, "ranking"))
    ranking_file = joinpath(parent_dir, "ranking", filename)
    # ranks = Dict{Drug, Vector{Int64}}()
    # we predict neg log values, i.e. higher value means lower concentration
    # which means higher susceptibility of the cell line to the drug and thus higher rank (lower number)
    # hence, we sort in reverse order
    for d in collect(keys(predictions))
        ranks[d] = zeros(Int64, length(predictions[d]))
        ranks[d][sortperm(predictions[d], rev=true)] = collect(1:length(predictions[d]))
    end
    f = open(ranking_file, "w")
    @printf(f, "DrugAnonID")
    for drug in collect(keys(pm.data.results))
        @printf(f, ",%s", drug.id)
    end
    @printf(f, "\n")
    for (cl_id, cl) in enumerate(collect(values(pm.data.cell_lines)))
        @printf(f, "%s", cl.id)
        for drug in collect(keys(pm.data.results))
           @printf(f, ",%d", ranks[drug][cl_id])
       end
       @printf(f, "\n")
    end
    close(f)

    # write model
    mkpath(joinpath(parent_dir, "model"))
    model_file = joinpath(parent_dir, "model", filename)
    f = open(model_file, "w")
    # T, K, N
    @printf(f, "T\tK\n")
    @printf(f, "%d\t%d\n", p.T, p.K)
    @printf(f, "N\n")
    N = collect(values(p.N))
    [@printf(f, "%d\t", n) for n in N[1:end-1]]
    @printf(f, "%d\n", N[end])
    # gamma
    @printf(f, "ɣ\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.ɣ[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.ɣ[end]))
    # b
    @printf(f, "b\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.b[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.b[end]))
    # a
    @printf(f, "a\n")
    for a_t in p.a
        @printf(f, "%s\n", a_t.var_name)
        a_t_exp = expected_value(a_t)
        [@printf(f, "%.5f\t", val) for val in a_t_exp[1:end-1]]
        @printf(f, "%.5f\n", a_t_exp[end])
    end
    # lambda
    @printf(f, "λ\n")
    for (t, λ_t) in enumerate(p.λ)
        λ_t_exp = expected_value.(λ_t)
        @printf(f, "λ[%d]:\t", t)
        for val in λ_t_exp[1:end-1]
            @printf(f, "%.5f\t", val)
        end
        @printf(f, "%.5f\n", λ_t_exp[end])
    end
    # epsilon
    @printf(f, "ε\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.ε[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.ε[end]))
    # nu
    @printf(f, "ν\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.ν[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.ν[end]))
    # G
    @printf(f, "G\n")
    for t in 1:p.T, k in 1:p.K
        exp_g = expected_value(p.G[t,k])
        @printf(f, "%s:\t", p.G[t,k].var_name)
        [@printf(f, "%.5f\t", val) for val in exp_g[1:end-1]]
        @printf(f, "%.5f\n", exp_g[end])
    end
    # omega
    @printf(f, "⍵\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.⍵[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.⍵[end]))
    # e
    @printf(f, "e\n")
    [@printf(f, "%.5f\t", expected_value(val)) for val in p.e[1:end-1]]
    @printf(f, "%.5f\n", expected_value(p.e[end]))

    close(f)

end
