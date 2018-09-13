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
		# @info "######### Drug $(drug.id) #########" target_ranking=sortperm(dep.targets[drug]) prediction_ranking=sortperm(y_mean) prediction_variance=y_cov gamma=expected_value(m.ɣ[t]) epsilon=expected_value(m.ε[t]) nu=expected_value(m.ν[t])
		ll += Distributions.logpdf(Distributions.MvNormal(y_mean, y_cov.*eye(m.N[t])), actual_outcomes)
		mse += sum((actual_outcomes .- y_mean).^2)/m.N[t]
	end
	# @info view_weights=expected_value.(m.e)
	(ll, mse)
end

function gridsearch(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, dest_path::String="results")
	gamma_dist_alphas = [1e-3, 1e-2, .1, 1., 10., 1e2, 1e3]
	gamma_dist_betas = [1e-3, 1e-2, .1, 1., 10., 1e2, 1e3]
	normal_means = [-1., 0., 1.]
	normal_vars = [1., 2., 5.]

	mkpath(joinpath(dest_path,"gridsearch_results"))
	mkpath(joinpath(dest_path,"gridsearch_charts"))
	@info "created paths"
	# Plots.plotlyjs()
	all_configurations = Vector{Vector{Float64}}()
	for alpha in gamma_dist_alphas, mu in normal_means, v in normal_vars
		push!(all_configurations, [alpha, mu, v])
	end
	@info "created $(length(all_configurations)) parameter settings"
    # (K, base_kernels, pathway_specific_kernels) = compute_all_kernels(dep.experiment, collect(values(dep.experiment.cell_lines)))
    
    # K = dep.model.K
    # base_kernels = dep.base_kernels
    # pathway_specific_kernels = dep.pathway_specific_kernels

    all_errors = Vector{String}(undef, length(all_configurations))

    # idx = Threads.Atomic{Int64}(1)
	f = open(joinpath(dest_path, "errors.txt"), "w")
	# for alpha in gamma_dist_alphas
	# 	for mu in normal_means, v in normal_vars
	@info "starting inference"
	Threads.@threads for i in 1:length(all_configurations)
		(alpha, mu, v) = all_configurations[i]
		@info "parameter setting $i" alpha mu v
		try
			(lls, errs, test_errs, model) = parameter_inference(dep, convergence_criterion=1e-4, min_iter=10, 
						⍺_ɣ=alpha, β_ɣ=1. / alpha,
						⍺_λ=alpha, β_λ=1. / alpha,
						⍺_ε=alpha, β_ε=1. / alpha,
						⍺_ν=alpha, β_ν=1. / alpha,
						⍺_⍵=alpha, β_⍵=1. / alpha,
						μ_b=mu, 𝜎_0=v,
						μ_e=mu, 𝜎_e=v,
						μ_a=mu, Σ_a=v,
						μ_g=mu, Σ_g=v)
			all_errors[i] = @sprintf("%f\t%f\t%f\t%f\t%f\t%f\n", alpha, mu, v, errs[end], test_errs[end], lls[end])
			# @printf(f, "%f\t%f\t%f\t%f\t%f\n", alpha, mu, v, errs[end], test_errs[end])
			# flush(f)
			# p = Plots.plot(title="alpha = $alpha, mu = $mu, variance = $v")
			# Plots.plot!(p, errs, label="training error")
			# Plots.plot!(p, test_errs, label="test error")
			# Plots.pdf(p, "gridsearch_charts/alpha_$(alpha)_mean_$(mu)_var_$(v).txt", p)
			predictions = predict_outcomes(dep, model, collect(values(dep.experiment.cell_lines)))
			write_prediction_file(joinpath(dest_path,"gridsearch_results","alpha_$(alpha)_mean_$(mu)_var_$(v).txt"), dep, predictions)
			# @info "Processed settings:" alpha mu variance=v training_error=errs[end] test_error=test_errs[end]
		catch exc
			# display(stacktrace(catch_backtrace()))
			# @warn "Exception occurred" exc alpha mu variance=v
		end
	end
	for s in all_errors
		@printf(f, "%s", s)
	end
	close(f)
end


function parameter_inference(dep::DrugEfficacyPredictor.DrugEfficacyPrediction; 
					convergence_criterion::Float64 = 1e-3, 
					min_iter::Int64 = 3, 
					⍺_ɣ::Float64=1., β_ɣ::Float64=1.,
					⍺_λ::Float64=1., β_λ::Float64=1.,
					⍺_ε::Float64=1., β_ε::Float64=1.,
					⍺_ν::Float64=1., β_ν::Float64=1.,
					⍺_⍵::Float64=1., β_⍵::Float64=1.,
					μ_b::Float64=0., 𝜎_0::Float64=20.,
					μ_e::Float64=0., 𝜎_e::Float64=20.,
					μ_a::Float64=0., Σ_a::Float64=20.,
					μ_g::Float64=0., Σ_g::Float64=20.)
	# Plots.plotly()
	model = DrugEfficacyPredictor.PredictionModel(dep.T, dep.K, dep.N,
					⍺_ɣ=⍺_ɣ, β_ɣ=β_ɣ,
					⍺_λ=⍺_λ, β_λ=β_λ,
					⍺_ε=⍺_ε, β_ε=β_ε,
					⍺_ν=⍺_ν, β_ν=β_ν,
					⍺_⍵=⍺_⍵, β_⍵=β_⍵,
					μ_b=μ_b, 𝜎_0=𝜎_0,
					μ_e=μ_e, 𝜎_e=𝜎_e,
					μ_a=μ_a, Σ_a=Σ_a,
					μ_g=μ_g, Σ_g=Σ_g)
	@info "created model..."
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
	@info "computed kernel products"
	old_ll = -1e100
	old_err = 1e100
	ll = 0.
	err = 0.
	convergence = 1
	err_convergence = 1
	iter = 1
	lls = Float64[]
	errs = Float64[]
	test_errs = Float64[]
	@info "entering optimization loop"
	while err_convergence > convergence_criterion || iter < min_iter
		ll, err = parameter_inference_step(dep, model, kernel_products, all_tasks)
		convergence = (old_ll - ll)/old_ll
		err_convergence = (old_err - err)/old_err
		# @info current_likelihood=ll current_error=err convergence error_convergence=err_convergence
		# break
		old_ll = ll
		old_err = err
		iter += 1
		push!(lls,ll)
		push!(errs,err)
		(test_err, _, _) = test(dep, model)
		push!(test_errs, test_err)
	end
	# println("highest log likelihood at iteration $(sortperm(lls, rev=true)[1])")
	# println("lowest error at iteration $(sortperm(errs)[1])")
	(lls, errs, test_errs, model)
end

# actual variational inference algorithm
function parameter_inference_step(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, m::PredictionModel, kernel_products::Dict{DrugEfficacyPredictor.Drug, Matrix{Float64}}, all_tasks::Vector{DrugEfficacyPredictor.Drug})
	@info "entering parameter inference single step"
    cell_lines = collect(values(dep.experiment.cell_lines))
	# @info "compute intermed kernel sums"
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
		# @info λ[t]=expected_value.(m.λ[t])
		ll += sum(elbo.(m.λ[t]))

		# a
		ν_expected = expected_value(m.ν[t])
		m.a[t].variational_covariance = inv(diagm(expected_value.(m.λ[t])) + ν_expected.*kernel_products[drug])
		m.a[t].variational_mean = m.a[t].variational_covariance*(ν_expected .* g_times_kernel[drug])
		# @info a[t]=expected_value(m.a[t])
		ll += elbo(m.a[t])

		# gamma
		m.ɣ[t].variational_a = m.ɣ[t].prior_a + .5
		m.ɣ[t].variational_b = m.ɣ[t].prior_b + .5*expected_squared_value(m.b[t])
		# @info ɣ[t]=expected_value(m.ɣ[t])
		ll += elbo(m.ɣ[t])

		# b
		ε_expected = expected_value(m.ε[t])
		m.b[t].variational_variance = 1. / (expected_value(m.ɣ[t]) + ε_expected*m.N[t])
		m.b[t].variational_mean = m.b[t].variational_variance*ε_expected*sum(dep.targets[drug] - exp_gte)
		# @info b[t]=expected_value(m.b[t])
		ll += elbo(m.b[t])

		# nu
		m.ν[t].variational_a = m.ν[t].prior_a + .5*(m.K*m.N[t])
		m.ν[t].variational_b = m.ν[t].prior_b 
						+ .5*tr(sum(expected_squared_value.(m.G[t,:]))) 
						- sum(transpose.(expected_value.(m.G[t,:])) .* dep.kernels[drug])*expected_value(m.a[t]) 
						+ .5*tr(kernel_products[drug]*expected_squared_value(m.a[t]))

		# sum_g_minus_kernel_a = sum(
		# 		tr.(expected_squared_value.(m.G[t,:])) 
		# 		.- dot.(expected_value.(m.G[t,:]), dep.kernels[drug] .* expected_value(m.a[t])) 
		# 		.+ tr.((*).(dep.kernels[drug], dep.kernels[drug]) .* expected_squared_value(m.a[t]))
		# 		)
		# sum_g_minus_kernel_a = 0
		# for k in 1:length(dep.kernels[drug])
		# 	exp_G = expected_value(m.G[t,k])
		# 	exp_G_sq = expected_squared_value(m.G[t,k])
		# 	exp_a = expected_value(m.a[t])
		# 	exp_a_sq = expected_squared_value(m.a[t])
		# 	sum_g_minus_kernel_a += tr(exp_G_sq) - 2*dot(exp_G, dep.kernels[drug][k]*exp_a) + tr(kernel_products[drug] * exp_a_sq)
		# end
		# m.ν[t].variational_b = m.ν[t].prior_b + .5*sum_g_minus_kernel_a
		# @info ν[t]=expected_value(m.ν[t])
		ll += elbo(m.ν[t])

		# epsilon
		m.ε[t].variational_a = m.ε[t].prior_a + .5*m.N[t]
		eps_beta_update = dot(dep.targets[drug], dep.targets[drug])
		eps_beta_update -= 2*dot(dep.targets[drug], exp_gte .+ expected_value(m.b[t]))
		eps_beta_update += sum(expected_squared_value.(m.e) .* tr.(expected_squared_value.(m.G[t,:])))
		eps_beta_update -= 2*sum(exp_gte .* expected_value(m.b[t]))
		eps_beta_update += m.N[t]*expected_squared_value(m.b[t])
		# exp_g_sq_sum = zeros(K, K)
		# for k in 1:m.K
		# 	exp_sq = expected_squared_value(m.G[t,k])
		# 	for i in 1:m.N[t], j in 1:m.N[t]
		# 		exp_g_sq_sum[i,j] += exp_sq[i,i]*exp_sq[j,j]
		# 	end
		# end
		# exp_g_sq_sum = zeros(m.N[t], m.N[t])			
		# [exp_g_sq_sum += expected_squared_value(m.G[t,k]) for k in 1:length(dep.kernels[drug])]
		# m.ε[t].variational_b = m.ε[t].prior_b + .5*(sum(expected_squared_value.(m.e))*tr(exp_g_sq_sum) + 2*sum(g_times_e[drug].*expected_value(m.b[t])) + expected_squared_value(m.b[t]))
		m.ε[t].variational_b = m.ε[t].prior_b + .5*eps_beta_update
		# @info ε[t]=expected_value(m.ε[t])
		ll += elbo(m.ε[t])
	end
	for k in 1:m.K
		# omega
		m.⍵[k].variational_a = m.⍵[k].prior_a + .5
		m.⍵[k].variational_b = m.⍵[k].prior_b + .5*expected_squared_value(m.e[k])
		ll += elbo(m.⍵[k])

		# e
		eta1 = sum(expected_value.(m.ε) .* dot.(collect(values(dep.targets)) .- expected_value.(m.b), expected_value.(m.G[:,k])))
		# exp_g_sq_sum = 0
		exp_g_sq_sum = sum(expected_value.(m.ε) .* tr.(expected_squared_value.(m.G[:,k])))
		# exp_g_times_y_minus_b = 0
		# for (t, drug) in enumerate(all_tasks)
		# 	# exp_g_sq_sum += tr(expected_squared_value(m.G[t,k]))
		# 	exp_g_times_y_minus_b += dot(expected_value(m.G[t,k]), dep.targets[drug] .- expected_value(m.b[t]))
		# end
		m.e[k].variational_variance = 1. / (expected_value(m.⍵[k]) + exp_g_sq_sum)
		m.e[k].variational_mean = m.e[k].variational_variance * eta1
		# @info e[k]=expected_value(m.e[k])
		ll += elbo(m.e[k])
	end
	# G
	for (t, drug) in enumerate(all_tasks)
		exp_g_sum = sum(expected_value.(m.G[t,:]).*expected_value.(m.e))
		# exp_g_sum = zeros(m.N[t])
		# [exp_g_sum += expected_value(m.G[t,k]) for k in 1:m.K]
		for (k, kernel) in enumerate(dep.kernels[drug])
			# G
			# NOTE: alter prior to allow correlations between cell lines, i.e. have a matrix prior on nu
			exp_e = expected_value(m.e[k])
			exp_g_sum -= expected_value(m.G[t,k])*exp_e
			m.G[t,k].variational_covariance = eye(m.N[t]) .* 1. / (2*expected_value(m.ε[t])*expected_squared_value(m.e[k]) + expected_value(m.ν[t]))

			eta1 = expected_value(m.ν[t]).*kernel*expected_value(m.a[t]) + expected_value(m.ε[t])*exp_e*(dep.targets[drug] .- expected_value(m.b[t]) - exp_g_sum)
			m.G[t,k].variational_mean = m.G[t,k].variational_covariance*eta1
			ll += elbo(m.G[t,k])
			exp_g_sum += expected_value(m.G[t,k])*exp_e
			# @info G[t,k]=expected_value(m.G[t,k])
		end
		#TODO: compute log likelihood E_q[ln p(gamma)] - E_q[ln q(gamma)]
	end
	(dll, err) = data_likelihood(dep, m)
	ll += dll
	ll, err
end

function update_intermed_times_weights(dep::DrugEfficacyPredictor.DrugEfficacyPrediction)
	ret = Dict{Drug, Vector{Float64}}()
	# @info "update_intermed_times_weights"
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_times_weight_sum = zeros(dep.model.N[t])
		for k in 1:length(dep.kernels[drug])
			g_times_weight_sum += expected_value(dep.model.G[t,k]).*expected_value(dep.model.e[k])
		end
		ret[drug] = g_times_weight_sum
		# info drug=g_times_weight_sum
	end
	ret
end

function update_intermed_kernel_sum(dep::DrugEfficacyPredictor.DrugEfficacyPrediction, model::PredictionModel)
	ret = Dict{Drug, Vector{Float64}}()
	# @info "update_intermed_kernel_sum"
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_kernel_sum = zeros(model.N[t])
		for (k, kernel) in enumerate(dep.kernels[drug])
			g_kernel_sum += kernel*expected_value(model.G[t,k])
		end
		ret[drug] = g_kernel_sum
		# @info drug=g_kernel_sum
	end
	ret
end