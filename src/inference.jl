import DrugEfficacyPredictor.DrugEfficacyPrediction
import DrugEfficacyPredictor.Drug
import DrugEfficacyPredictor.expected_value

function data_likelihood(dep::DrugEfficacyPredictor.DrugEfficacyPrediction)
	m = dep.model
	ll = 0.
	for (t, drug) in enumerate(keys(dep.experiment.results))
		y_mean = sum(expected_value.(m.G[t,:]) .* expected_value.(m.e) .+ expected_value(m.b[t]))
		y_cov = expected_value(m.ε[t]).*eye(m.N[t])

		ll += Distributions.logpdf(Distributions.MvNormal(y_mean, y_cov), values(dep.experiment.results[drug].normalized_outcome_values))
	end
	ll
end


function parameter_inference(dep::DrugEfficacyPrediction)
	old_ll = -1e100
	ll = 0
	convergence = 1
	while convergence > 1e-3
		ll = parameter_inference_step(dep)
		convergence = (old_ll - ll)/old_ll
		info("current likelihood: $ll, convergence: $convergence")
		old_ll = ll
	end
	ll
end

# actual variational inference algorithm
function parameter_inference_step(dep::DrugEfficacyPrediction)
	m = dep.model
	all_tasks = collect(keys(dep.experiment.results))

	kernel_products = Dict{Drug, Matrix{Float64}}()
    cell_lines = collect(values(dep.experiment.cell_lines))

	for (t, d) in enumerate(all_tasks)
		views = dep.kernels[d]
		kp = zeros(m.N[t], m.N[t])
		info("kernel product: $(size(kp))")
		for kernel in views
			info("kernel: $(size(kernel))")
			info("mult kernel: $(size(kernel'*kernel))")
			# info("length found idx: $(length(idx))")
			# kp[idx, idx] += kernel'*kernel
			kp += kernel'*kernel
		end
		kernel_products[d] = kp
	end
	g_times_kernel = update_intermed_kernel_sum(dep)
	g_times_e = update_intermed_times_weights(dep)
	# updates for model parameters in turn
	ll = 0
	for (t, drug) in enumerate(all_tasks)
		# lambda
		a_square = expected_squared_value(m.a[t])
		for n in 1:m.N[t]
			m.λ[t][n].variational_a = m.λ[t][n].prior_a + .5
			m.λ[t][n].variational_b = m.λ[t][n].prior_b + .5*a_square[n][n]
		end
		ll += sum(elbo.(m.λ[t]))

		# a
		ν_expected = expected_value(m.ν[t])
		m.a[t].variational_covariance = inv(diagm(expected_value.(m.λ[t])) + ν_expected.*kernel_products[drug])
		m.a[t].variational_mean = m.a[t].variational_covariance*(ν_expected .* g_times_kernel[drug])
		ll += elbo(m.a[t])

		# gamma
		m.ɣ[t].variational_a = m.ɣ[t].prior_a + .5
		m.ɣ[t].variational_b = m.ɣ[t].prior_b + .5*expected_squared_value(m.b[t])
		ll += elbo(m.ɣ[t])

		# b
		ε_expected = expected_value(m.ε[t])
		m.b[t].variational_variance = 1./(expected_value(m.ɣ[t]) + ε_expected*m.N[t])
		m.b[t].variational_mean = m.b[t].variational_variance*ε_expected*(dep.targets[drug] - g_times_e[drug])
		ll += elbo(m.b[t])

		# nu
		m.ν[t].variational_a = m.ν[t].prior_a + .5*(m.K*m.N[t])
		sum_g_minus_kernel_a = 0
		for k in 1:length(dep.kernels[drug])
			exp_G = expected_value(m.G[t,k])
			exp_G_sq = expected_squared_value(m.G[t,k])
			exp_a = expected_value(m.a[t])
			exp_a_sq = expected_squared_value(m.a[t])
			sum_g_minus_kernel_a += trace(exp_G_sq) - 2*dot(exp_G, dep.kernels[drug][k]*exp_a) + trace(kernel_products[drug] * exp_a_sq)
		end
		m.ν[t].variational_b = m.ν[t].prior_b + .5*sum_g_minus_kernel_a
		ll += elbo(m.ν[t])

		# epsilon
		m.ε[t].variational_a = m.ε[t].prior_a + .5*m.N[t]
		exp_g_sq_sum = zeros(m.N[t], m.N[t])
		[exp_g_sq_sum += expected_squared_value(m.G[t,k]) for k in 1:length(dep.kernels[drug])]
		m.ε[t].variational_b = m.ε[t].prior_b + .5*(trace(exp_g_sq_sum*diagm(expected_squared_value.(m.e))) + 2*sum(g_times_e[drug].*expected_value(m.b[t])) + expected_squared_value(m.b[t]))
		ll += elbo(m.ε[t])
	end
	for k in 1:m.K
		# omega
		m.⍵[k].variational_a = m.⍵[k].prior_a + .5
		m.⍵[k].variational_b = m.⍵[k].prior_b + .5*expected_squared_value(m.e[k])
		ll += elbo(m.⍵[k])

		# e
		exp_g_sq_sum = 0
		exp_g_times_y_minus_b = 0
		for (t, drug) in enumerate(all_tasks)
			exp_g_sq_sum += trace(expected_squared_value(m.G[t,k]))
			exp_g_times_y_minus_b += dot(expected_value(m.G[t,k]), dep.targets[drug] .- expected_value(m.b[t]))
		end
		m.e[k].variational_variance = 1./(expected_value(m.⍵[k]) + exp_g_sq_sum)
		m.e[k].variational_mean = m.e[k].variational_variance * exp_g_times_y_minus_b
		ll += elbo(m.e[k])
	end

	for (t, d) in enumerate(all_tasks)
		exp_g_sum = zeros(m.N[t])
		[exp_g_sum += expected_value(m.G[t,k]) for k in 1:length(all_tasks)]
		for (k, kernel) in enumerate(dep.kernesl[drug])
			# G
			# NOTE: alter prior to allow correlations between cell lines, i.e. have a matrix prior on nu
			exp_g_sum -= expected_value(m.G[t,k])
			m.G[t,k].variational_covariance = inv(diagm(expected_value(m.ε[t])*expected_squared_value(m.e[k]) + expected_value(m.ν[t])))
			eta1 = expected_value(m.ν[t]).*kernel*expected_value(m.a[t]) + expected_value(m.ε[t])*expected_value(m.e[k])*(dep.targets[drug] .- expected_value(m.b[t]) -.5*exp_g_sum)
			m.G[t,k].variational_mean = m.G[t,k].variational_covariance*eta1
			ll += elbo(m.G[t,k])
			exp_g_sum += expected_value(m.G[t,k])
		end
		#TODO: compute log likelihood E_q[ln p(gamma)] - E_q[ln q(gamma)]
	end
	ll += data_likelihood(dep)
	ll
end

function update_intermed_times_weights(dep::DrugEfficacyPrediction)
	ret = Dict{Drug, Vector{Float64}}()
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_times_weight_sum = zeros(dep.model.N[t])
		for k in 1:length(dep.kernels[drug])
			g_times_weight_sum += expected_value(dep.model.G[t,k]).*expected_value(dep.model.ε[t])
		end
		ret[drug] = g_times_weight_sum
	end
	ret
end

function update_intermed_kernel_sum(dep::DrugEfficacyPrediction)
	ret = Dict{Drug, Vector{Float64}}()
	info("G: $(size(dep.model.G))")
	for (t, drug) in enumerate(collect(keys(dep.experiment.results)))
		g_kernel_sum = zeros(dep.model.N[t])
		info("g kernel sum: $(size(g_kernel_sum))")
		for (k, kernel) in enumerate(dep.kernels[drug])
			info("G[t,k]: $(size(expected_value(dep.model.G[t,k]))), kernel: $(size(kernel))")
			g_kernel_sum += kernel*expected_value(dep.model.G[t,k])
		end
		ret[drug] = g_kernel_sum
	end
	ret
end