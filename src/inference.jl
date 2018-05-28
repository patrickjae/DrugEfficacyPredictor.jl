function data_likelihood(dep::DrugEfficacyPrediction)
	m = dep.model
	ll = 0.
	for (t, drug) in enumerate(keys(dep.experiment.results))
		y_mean = sum(expected_value.(m.G[t,:]) .* expected_value.(m.e) .+ expected_value(m.b[t]))
		y_cov = expected_value(m.ε[t]).*eye(m.N[t])

		ll += Distributions.logpdf(Distributions.MvNormal(y_mean, y_cov), values(dep.experiment.results[drug].normalized_outcome_values))
	end
	ll
end


# actual variational inference algorithm
function parameter_inference(dep::DrugEfficacyPrediction)
	m = dep.model
	all_tasks = collect(keys(dep.experiment.results))

	kernel_products = Dict{Drug, Matrix{Float64}}()
	for (t, d) in enumerate(all_tasks)
		views = dep.kernels[d]
		kp = zeros(m.N[t], m.N[t])
		[kp += kernel'*kernel for kernel in views]
		kernel_products[d] = kp
	end
	# updates for model parameters in turn
	ll = 0
	# lambda, a
	for (t, drug) in enumerate(all_tasks)
		aaT = expected_squared_value(m.a[t])
		bbT = expected_squared_value(m.b[t])
		for n in 1:m.N[t]
			m.λ[t][n].variational_a = m.λ[t][n].prior_a + .5
			m.λ[t][n].variational_b = m.λ[t][n].prior_b + .5*aaT[n,n]
		end
		ll += elbo.(m.λ[t])

		m.a[t].variational_covariance = inv(diagm(expected_value.(m.λ[t])) + expected_value(m.ν[t]).*kernel_products[drug])
		
		# m.a[t].variational_mean = m.a[t].variational_covariance*
	end
	# gamma
	for t in 1:m.T
		m.ɣ.variational_a = m.ɣ.prior_a + .5
		m.ɣ.variational_b = m.ɣ.prior_b + expected_squared_value(m.b[t])
		#TODO: compute log likelihood E_q[ln p(gamma)] - E_q[ln q(gamma)]
	end


end
