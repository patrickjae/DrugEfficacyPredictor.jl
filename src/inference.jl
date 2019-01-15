import DrugEfficacyPredictor.Drug
import DrugEfficacyPredictor.expected_value
import DrugEfficacyPredictor.expected_squared_value
import DrugEfficacyPredictor.elbo

# using Plots, Rsvg

function data_likelihood(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, m::PredictionModel)
	ll = 0.
	mse = 0
	for (t, drug) in enumerate(keys(dep.experiment.results))
		y_mean = sum(expected_value.(m.G[t,:]) .* expected_value.(m.e)) .+ expected_value(m.b[t])

		y_cov = 1. / expected_value(m.ε[t])

		actual_outcomes = dep.targets[drug]
		y_mean_rescaled = y_mean .* dep.experiment.results[drug].outcome_std .+ dep.experiment.results[drug].outcome_mean
		outcomes_rescaled = actual_outcomes .* dep.experiment.results[drug].outcome_std .+ dep.experiment.results[drug].outcome_mean

		ll += Distributions.logpdf(Distributions.MvNormal(y_mean, y_cov.*Matrix(I, m.N[t], m.N[t])), actual_outcomes)

		mse += sum((outcomes_rescaled .- y_mean_rescaled).^2)
	end
	mse /= sum(m.N)
	(ll, mse)
end

# function gridsearch(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, dest_path::String="results")
# 	# gamma_dist_alphas = [1e-10, 1e-5, .1, 1., 10., 1e5, 1e10]
# 	# gamma_dist_betas = [1e-10, 1e-5, .1, 1., 10., 1e5, 1e10]
# 	gamma_dist_alphas = [250., 500., 750., 1000.]
# 	alpha_beta_ratios = [1., 5., 50., 100.]
# 	# gamma_dist_betas = [1e-10, .1, 1., 1e10]
# 	# gamma_dist_alphas = [1.]
# 	# gamma_dist_betas = [1.]
# 	# normal_means = [-1., 0., 1.]
# 	normal_vars = [.1, .5, 1., 2., 10., 20.]
# 	# gamma_dist_alphas = [1e-3]
# 	# gamma_dist_betas = [1e-3]
# 	normal_means = [1.]
# 	# normal_vars = [2.]
#
# 	mkpath(dest_path)
# 	log_message("created paths")
# 	# Plots.plotlyjs()
# 	all_configurations = Vector{Vector{Float64}}()
# 	for alpha in gamma_dist_alphas, ratio in alpha_beta_ratios, v in normal_vars
# 		push!(all_configurations, [alpha, alpha/ratio, 1., v])
# 	end
# 	log_message("created $(length(all_configurations)) parameter settings")
#
#     all_errors = String[]
#
#     inference_config = InferenceConfiguration()
#     inference_config.convergence_criterion = 1e-4
#     inference_config.max_iter = 200
#     inference_config.target_dir = dest_path
#
#     set_training_test_kernels(dep)
# 	log_message("starting inference")
# 	for i in 1:length(all_configurations)
# 		(alpha, beta, mu, v) = all_configurations[i]
# 		log_message("parameter setting $i, alpha=$alpha beta=$beta mu=$mu v=$v")
# 		mc = ModelConfiguration(alpha, beta, mu, v)
# 		mc.μ_b = 0.
# 		mc.μ_e = 1.
# 		try
# 			(lls, errs, test_errs, model, convergence) = parameter_inference(dep, inference_config = inference_config, model_config = mc)
# 			log_message("inference stats: likelihood=$(lls[end]), train_error=$(errs[end]) test_error=$(test_errs[end]) convergence=$convergence")
# 			push!(all_errors, "$alpha\tbeta\tmu\tv\t$(errs[end])\t$(test_errs[end])\t$(lls[end])\n")
# 			(predictions, ranks) = predict_outcomes(dep, model, collect(values(dep.experiment.cell_lines)))
# 			write_results(dep, inference_config.target_dir, "alpha_$(alpha)_beta_$(beta)_mean_$(mu)_var_$(v).txt", predictions, model)
# 		catch exc
# 			display(stacktrace(catch_backtrace()))
# 			log_message("Exception occurred, $exc alpha=$alpha mu=$mu variance=$v")
# 			break
# 		end
# 	end
#
# 	f = open(joinpath(dest_path, "errors.txt"), "w")
# 	for s in all_errors
# 		@printf(f, "%s", s)
# 	end
# 	close(f)
# end


function parameter_inference(dep::DrugEfficacyPredictor.DrugEfficacyPrediction; inference_config::InferenceConfiguration = InferenceConfiguration(), model_config::ModelConfiguration)
	model = DrugEfficacyPredictor.PredictionModel(dep.T, dep.K, dep.N, model_config = model_config)
	all_tasks = collect(keys(dep.experiment.results))
	kernel_products = Dict{DrugEfficacyPredictor.Drug, Matrix{Float64}}()
	for (t, d) in enumerate(all_tasks)
		views = dep.kernels[d]
		kp = zeros(model.N[t], model.N[t])
		for kernel in views
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
		ll, err = parameter_inference_step(dep, model, kernel_products, all_tasks)
		convergence = (old_ll - ll)/old_ll
		# convergence = (old_err - err)/old_err
		# break
		old_ll = ll
		old_err = err
		iter += 1
		push!(lls,ll)
		push!(errs,err)
		(test_err, _, _) = test(dep, model)
		push!(test_errs, test_err)
		if iter > inference_config.max_iter
			log_message("reached maxium number of iterations ($(inference_config.max_iter)), terminating")
			break
		end
	end
	(lls, errs, test_errs, model, convergence)
end

# actual variational inference algorithm
function parameter_inference_step(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, m::PredictionModel, kernel_products::Dict{DrugEfficacyPredictor.Drug, Matrix{Float64}}, all_tasks::Vector{DrugEfficacyPredictor.Drug})
    cell_lines = collect(values(dep.experiment.cell_lines))
	g_times_kernel = update_intermed_kernel_sum(dep, m)
	# updates for model parameters in turn
	ll = 0
	for (t, drug) in enumerate(all_tasks)
		exp_gte = sum(expected_value.(m.e) .* expected_value.(m.G[t,:]))
		# lambda
		a_square = expected_squared_value(m.a[t])
		for n in 1:m.N[t]
			m.λ[t][n].variational_a = m.λ[t][n].prior_a + .5
			m.λ[t][n].variational_b = m.λ[t][n].prior_b + .5*a_square[n,n]
		end
		ll += sum(elbo.(m.λ[t]))

		# a
		ν_expected = expected_value(m.ν[t])
		m.a[t].variational_covariance = inv(diagm(0 => expected_value.(m.λ[t])) + ν_expected.*kernel_products[drug])
		m.a[t].variational_mean = m.a[t].variational_covariance*(ν_expected .* g_times_kernel[drug])
		ll += elbo(m.a[t])

		# gamma
		m.ɣ[t].variational_a = m.ɣ[t].prior_a + .5
		m.ɣ[t].variational_b = m.ɣ[t].prior_b + .5*expected_squared_value(m.b[t])
		ll += elbo(m.ɣ[t])

		# b
		ε_expected = expected_value(m.ε[t])
		m.b[t].variational_variance = 1. / (expected_value(m.ɣ[t]) + ε_expected*m.N[t])
		m.b[t].variational_mean = m.b[t].variational_variance*ε_expected*sum(dep.targets[drug] .- exp_gte)
		ll += elbo(m.b[t])

		# nu
		m.ν[t].variational_a = m.ν[t].prior_a + .5*(m.K*m.N[t])
		m.ν[t].variational_b = m.ν[t].prior_b
						+ .5*tr(sum(expected_squared_value.(m.G[t,:])))
						- sum(transpose.(expected_value.(m.G[t,:])) .* dep.kernels[drug])*expected_value(m.a[t])
						+ .5*tr(kernel_products[drug]*expected_squared_value(m.a[t]))
		ll += elbo(m.ν[t])

		# epsilon
		m.ε[t].variational_a = m.ε[t].prior_a + .5*m.N[t]
		eps_beta_update = dot(dep.targets[drug], dep.targets[drug])
		eps_beta_update -= 2*dot(dep.targets[drug], exp_gte .+ expected_value(m.b[t]))
		eps_beta_update += sum(expected_squared_value.(m.e) .* tr.(expected_squared_value.(m.G[t,:])))
		eps_beta_update -= 2*sum(exp_gte .* expected_value(m.b[t]))
		eps_beta_update += m.N[t]*expected_squared_value(m.b[t])
		m.ε[t].variational_b = m.ε[t].prior_b + .5*eps_beta_update
		ll += elbo(m.ε[t])
	end
	for k in 1:m.K
		# omega
		m.⍵[k].variational_a = m.⍵[k].prior_a + .5
		m.⍵[k].variational_b = m.⍵[k].prior_b + .5*expected_squared_value(m.e[k])
		ll += elbo(m.⍵[k])

		# e
		eta1 = sum(expected_value.(m.ε) .* dot.([el .- expected_value(m.b[i]) for (i, el) in enumerate(collect(values(dep.targets)))], expected_value.(m.G[:,k])))
		exp_g_sq_sum = sum(expected_value.(m.ε) .* tr.(expected_squared_value.(m.G[:,k])))
		m.e[k].variational_variance = 1. / (expected_value(m.⍵[k]) + exp_g_sq_sum)
		m.e[k].variational_mean = m.e[k].variational_variance * eta1
		ll += elbo(m.e[k])
	end
	# G
	for (t, drug) in enumerate(all_tasks)
		exp_g_sum = sum(expected_value.(m.G[t,:]).*expected_value.(m.e))
		for (k, kernel) in enumerate(dep.kernels[drug])
			# G
			# NOTE: alter prior to allow correlations between cell lines, i.e. have a matrix prior on nu
			exp_e = expected_value(m.e[k])
			exp_g_sum -= expected_value(m.G[t,k])*exp_e
			m.G[t,k].variational_covariance = Matrix(I, m.N[t], m.N[t]) .* 1. / (2*expected_value(m.ε[t])*expected_squared_value(m.e[k]) + expected_value(m.ν[t]))

			eta1 = expected_value(m.ν[t]).*kernel*expected_value(m.a[t]) + expected_value(m.ε[t])*exp_e*(dep.targets[drug] .- expected_value(m.b[t]) .- exp_g_sum)
			m.G[t,k].variational_mean = m.G[t,k].variational_covariance*eta1
			ll += elbo(m.G[t,k])
			exp_g_sum += expected_value(m.G[t,k])*exp_e
		end
	end
	(dll, err) = data_likelihood(dep, m)
	ll += dll
	ll, err
end

function update_intermed_times_weights(dep::DrugEfficacyPredictor.DrugEfficacyPrediction)
	ret = Dict{Drug, Vector{Float64}}()
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_times_weight_sum = zeros(dep.model.N[t])
		for k in 1:length(dep.kernels[drug])
			g_times_weight_sum += expected_value(dep.model.G[t,k]).*expected_value(dep.model.e[k])
		end
		ret[drug] = g_times_weight_sum
	end
	ret
end

function update_intermed_kernel_sum(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, model::PredictionModel)
	ret = Dict{Drug, Vector{Float64}}()
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_kernel_sum = zeros(model.N[t])
		for (k, kernel) in enumerate(dep.kernels[drug])
			g_kernel_sum += kernel*expected_value(model.G[t,k])
		end
		ret[drug] = g_kernel_sum
	end
	ret
end
