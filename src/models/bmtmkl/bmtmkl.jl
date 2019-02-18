include("types.jl")
include("bmtmkl_init.jl")

function inference(m::PredictionModel{BMTMKLModel}, mc::ModelConfiguration; inference_config::InferenceConfiguration = InferenceConfiguration())
	params = BMTMKLModelParameters(mc)
	kernel_products = Dict{Drug, Matrix{Float64}}()
	for (t, d) in enumerate(keys(m.data.results))
		kernels = m.precomputations["kernels"][d]
		kp = zeros(params.N[d], params.N[d])
		for kernel in kernels
			kp += kernel'*kernel
		end
		kernel_products[d] = kp
	end
	old_ll = -1e100
	old_err = 1e100
	ll = 0.
	err = 0.
	convergence = 1
	iter = 1
	lls = Float64[]
	errs = Float64[]
	test_errs = Float64[]
	# while err_convergence > convergence_criterion || iter < min_iter
	while convergence > inference_config.convergence_criterion || iter < inference_config.min_iter
		ll, err = parameter_inference_step(m, params, kernel_products)
		convergence = (old_ll - ll)/old_ll
		# convergence = (old_err - err)/old_err
		# break
		old_ll = ll
		old_err = err
		iter += 1
		push!(lls,ll)
		push!(errs,err)
		(test_err, _, _) = test(m, params)
		push!(test_errs, test_err)
		if iter > inference_config.max_iter
			log_message("reached maxium number of iterations ($(inference_config.max_iter)), terminating")
			break
		end
	end
	(lls, errs, test_errs, params, convergence)
end

function test(m::PredictionModel{BMTMKLModel}, params::BMTMKLModelParameters)
    mse = 0

    rankings = Dict{Drug, Vector{Int64}}()
    predictions = Dict{Drug, Vector{Float64}}()
    sum_data_points = 0
    for (t, drug) in enumerate(keys(m.data.results))
		kernels = m.precomputations["cross_kernels"][drug]
		G = Vector{Vector{Float64}}(undef, length(kernels))
		exp_a = expected_value(params.a[t])
		[G[k] = kernel * exp_a for (k, kernel) in enumerate(kernels)]

        y_mean = sum(G .* expected_value.(params.e)) .+ expected_value(params.b[t])

        # y_mean = sum([G[k] .* expected_value(params.e[k]) for k in 1:length(G)]) .+ expected_value(params.b[t])

        #########################
        # compute the ranking
        #########################
        # the ranking is relative to each drug, indicating which cell line is how likely to be affected by the drug, i.e we rank the cell lines, not the drugs
        # a lower rank means that the drug has higher effect on the cell line
        # WARNING: this is somewhat against the final use, since we are interested in ranking drugs conditioned on cell lines/tissues
        # lower IC50 is better, so lowest ic50 (or whatever measure) results in highest rank (actually, the lowest number)
        # !!!!!!!!  BUT we are dealing with the -log(IC50) so higher value results in higher rank
        #########################
        rankings[drug] = zeros(Int64, length(y_mean))
        rankings[drug][sortperm(y_mean, rev=true)] = collect(1:length(y_mean))

        y_mean_rescaled = y_mean .* m.data.results[drug].outcome_std .+ m.data.results[drug].outcome_mean
        predictions[drug] = y_mean_rescaled

        # @info "######### Drug $(drug.id) #########" targets=dep.targets[drug] predictions=y_mean_rescaled target_ranking prediction_ranking=rankings[drug] prediction_variance=(1 ./ expected_value(m.ε[t])) gamma=expected_value(m.ɣ[t]) nu=expected_value(m.ν[t])

        #get unnormalized test targets
        actual_outcomes = m.precomputations["test_targets"][drug]
        mse += sum((actual_outcomes .- y_mean_rescaled).^2)
        sum_data_points += length(actual_outcomes)
        # @info "current test error setting" actual_outcomes y_mean_rescaled mse=sum((actual_outcomes .- y_mean_rescaled).^2)/length(actual_outcomes)
    end
    mse /= sum_data_points
    # @info view_weights=expected_value.(m.e))
    (mse, predictions, rankings)
end


function predict_outcomes(m::PredictionModel{BMTMKLModel}, params::BMTMKLModelParameters, cell_lines::Vector{CellLine}; held_out::Bool = false, inference_config::InferenceConfiguration = InferenceConfiguration())
    predictions = Dict{Drug, Vector{Float64}}()
    ranks = Dict{Drug, Vector{Int64}}()
    # all cell lines present for training
    all_cell_lines = collect(values(m.data.cell_lines))
    # only compute the kernels when dealing with held-out data
    base_cross_kernels = m.precomputations["base_kernels"]
    pw_specific_cross_kernels = m.precomputations["pathway_specific_kernels"]

	#if dealing with held out data, ignore validation data and recompute all kernels
    if held_out
        (base_cross_kernels, pw_specific_cross_kernels) = compute_all_kernels(m.data, all_cell_lines, cell_lines, subsume_pathways = inference_config.subsume_pathways)
    end

	num_cross_kernels = length(base_cross_kernels) + sum(map(length, collect(values(pw_specific_cross_kernels)))) + 4

    # do predictions for all drugs we saw at training time
    for (t, drug) in enumerate(keys(m.data.results))
        # find cell lines that have outcome for this drug and are in the training set
        training_cell_lines = filter(cl -> !cl.in_test_set && haskey(m.data.results[drug].outcome_values, cl), all_cell_lines)
        training_set_cell_line_idx = findall((in)(training_cell_lines), all_cell_lines)
        # find the ids of these cell_lines
        predict_cell_line_idx = collect(1:length(cell_lines))
        if !held_out
            # these cell lines should be included in the experiment, determine their position
            predict_cell_line_idx = findall((in)(cell_lines), all_cell_lines)
        end

        # compute similarity with all other cell lines in the training set
        kernels = Vector{Matrix{Float64}}()

        for v in m.data.views
            push!(kernels, base_cross_kernels[v][training_set_cell_line_idx, predict_cell_line_idx])
            if length(m.data.pathway_information) != 0
                for pw_kernel in pw_specific_cross_kernels[v]
                    push!(kernels, pw_kernel[training_set_cell_line_idx, predict_cell_line_idx])
                end
            end
        end

        gene_expression = base_cross_kernels[GeneExpression][training_set_cell_line_idx, :]
        methylation = base_cross_kernels[Methylation][training_set_cell_line_idx, :]
        cnv = base_cross_kernels[CNV][training_set_cell_line_idx, :]

        push!(kernels, gene_expression .* methylation)
        push!(kernels, gene_expression .* cnv)
        push!(kernels, cnv .* methylation)
        push!(kernels, gene_expression .* methylation .* cnv)

        G = Vector{Vector{Float64}}(undef, length(kernels))
        exp_a = expected_value(params.a[t])
        [G[k] = kernel' * exp_a for (k, kernel) in enumerate(kernels)]

        pred_y = sum(G .* expected_value.(params.e)) .+ expected_value(params.b[t])

        #compute ranking of predicted cell line responses
        ranks[drug] = zeros(Int64, length(pred_y))
        ranks[drug][sortperm(pred_y, rev=true)] = collect(1:length(pred_y))

        # rescale the normalized prediction
        pred_y_rescaled = pred_y * m.data.results[drug].outcome_std .+ m.data.results[drug].outcome_mean
        predictions[drug] = pred_y_rescaled
    end
    (predictions, ranks)
end


# actual variational inference algorithm
function parameter_inference_step(m::PredictionModel{BMTMKLModel}, params::BMTMKLModelParameters, kernel_products::Dict{Drug, Matrix{Float64}})
	g_times_kernel = update_intermed_kernel_sum(m, params)
	# updates for model parameters in turn
	ll = 0
	T = params.T
	N = params.N
	K = params.K
	for (t, drug) in enumerate(keys(m.data.results))
		exp_gte = sum(expected_value.(params.e) .* expected_value.(params.G[t,:]))
		# lambda
		a_square = expected_squared_value(params.a[t])
		for n in 1:N[drug]
			params.λ[t][n].variational_a = params.λ[t][n].prior_a + .5
			params.λ[t][n].variational_b = params.λ[t][n].prior_b + .5*a_square[n,n]
		end
		ll += sum(elbo.(params.λ[t]))

		# a
		ν_expected = expected_value(params.ν[t])
		params.a[t].variational_covariance = inv(diagm(0 => expected_value.(params.λ[t])) + ν_expected .* kernel_products[drug])
		params.a[t].variational_mean = params.a[t].variational_covariance*(ν_expected .* g_times_kernel[drug])
		ll += elbo(params.a[t])

		# gamma
		params.ɣ[t].variational_a = params.ɣ[t].prior_a + .5
		params.ɣ[t].variational_b = params.ɣ[t].prior_b + .5*expected_squared_value(params.b[t])
		ll += elbo(params.ɣ[t])

		# b
		ε_expected = expected_value(params.ε[t])
		params.b[t].variational_variance = 1. / (expected_value(params.ɣ[t]) + ε_expected * N[drug])
		params.b[t].variational_mean = params.b[t].variational_variance * ε_expected * sum(m.precomputations["targets"][drug] .- exp_gte)
		ll += elbo(params.b[t])

		# nu
		params.ν[t].variational_a = params.ν[t].prior_a + .5*(K * N[drug])
		params.ν[t].variational_b = params.ν[t].prior_b +
						.5*tr(sum(expected_squared_value.(params.G[t,:]))) -
						sum(transpose.(expected_value.(params.G[t,:])) .* m.precomputations["kernels"][drug]) * expected_value(params.a[t]) +
						.5*tr(kernel_products[drug] * expected_squared_value(params.a[t]))
		ll += elbo(params.ν[t])

		# epsilon
		params.ε[t].variational_a = params.ε[t].prior_a + .5 * N[drug]
		eps_beta_update = dot(m.precomputations["targets"][drug], m.precomputations["targets"][drug])
		eps_beta_update -= 2 * dot(m.precomputations["targets"][drug], exp_gte .+ expected_value(params.b[t]))
		eps_beta_update += sum(expected_squared_value.(params.e) .* tr.(expected_squared_value.(params.G[t,:])))
		eps_beta_update -= 2 * sum(exp_gte .* expected_value(params.b[t]))
		eps_beta_update += N[drug] * expected_squared_value(params.b[t])
		params.ε[t].variational_b = params.ε[t].prior_b + .5 * eps_beta_update
		ll += elbo(params.ε[t])
	end
	for k in 1:K
		# omega
		params.⍵[k].variational_a = params.⍵[k].prior_a + .5
		params.⍵[k].variational_b = params.⍵[k].prior_b + .5*expected_squared_value(params.e[k])
		ll += elbo(params.⍵[k])

		# e
		eta1 = sum(
			expected_value.(params.ε) .*
			dot.(
				[
					el .- expected_value(params.b[i])
					for (i, el) in enumerate(collect(values(m.precomputations["targets"])))
				], expected_value.(params.G[:,k])
				)
			)
		exp_g_sq_sum = sum(expected_value.(params.ε) .* tr.(expected_squared_value.(params.G[:,k])))
		params.e[k].variational_variance = 1. / (expected_value(params.⍵[k]) + exp_g_sq_sum)
		params.e[k].variational_mean = params.e[k].variational_variance * eta1
		ll += elbo(params.e[k])
	end

	# G
	for (t, drug) in enumerate(keys(m.data.results))
		exp_g_sum = sum(expected_value.(params.G[t,:]) .* expected_value.(params.e))
		for (k, kernel) in enumerate(m.precomputations["kernels"][drug])
			# G
			# NOTE: alter prior to allow correlations between cell lines, i.e. have a matrix prior on nu
			exp_e = expected_value(params.e[k])
			exp_g_sum -= expected_value(params.G[t,k]) * exp_e
			params.G[t,k].variational_covariance = Matrix(I, N[drug], N[drug]) .*
				1. / (
					2 * expected_value(params.ε[t]) *
					expected_squared_value(params.e[k]) +
					expected_value(params.ν[t])
				)

			eta1 = expected_value(params.ν[t]) .* kernel * expected_value(params.a[t]) +
				expected_value(params.ε[t]) * exp_e * (
					m.precomputations["targets"][drug] .-
					expected_value(params.b[t]) .-
					exp_g_sum
				)
			params.G[t,k].variational_mean = params.G[t,k].variational_covariance * eta1
			ll += elbo(params.G[t,k])
			exp_g_sum += expected_value(params.G[t,k]) * exp_e
		end
	end

	(dll, err) = data_likelihood(m, params)

	ll += dll
	ll, err
end

function data_likelihood(m::PredictionModel{BMTMKLModel}, params::BMTMKLModelParameters)
	ll = 0.
	mse = 0
	for (t, drug) in enumerate(keys(m.data.results))
		y_mean = sum(expected_value.(params.G[t,:]) .* expected_value.(params.e)) .+ expected_value(params.b[t])

		y_cov = 1. / expected_value(params.ε[t])

		actual_outcomes = m.precomputations["targets"][drug]
		y_mean_rescaled = y_mean .* m.data.results[drug].outcome_std .+ m.data.results[drug].outcome_mean
		outcomes_rescaled = actual_outcomes .* m.data.results[drug].outcome_std .+ m.data.results[drug].outcome_mean

		ll += Distributions.logpdf(Distributions.MvNormal(y_mean, y_cov .* Matrix(I, params.N[drug], params.N[drug])), actual_outcomes)

		mse += sum((outcomes_rescaled .- y_mean_rescaled).^2)
	end
	mse /= sum(collect(values(params.N)))
	(ll, mse)
end

function update_intermed_kernel_sum(m::PredictionModel{BMTMKLModel}, params::BMTMKLModelParameters)
	ret = Dict{Drug, Vector{Float64}}()
	for (t, drug) in enumerate(keys(m.data.results))
		g_kernel_sum = zeros(params.N[drug])
		for (k, kernel) in enumerate(m.precomputations["kernels"][drug])
			g_kernel_sum += kernel * expected_value(params.G[t,k])
		end
		ret[drug] = g_kernel_sum
	end
	ret
end
