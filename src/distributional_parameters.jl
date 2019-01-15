using Distributions

abstract type ModelParameter end

mutable struct GammaParameter <: ModelParameter
	prior_a::Float64
	prior_b::Float64
	prior_distribution::Distributions.Gamma

	variational_a::Float64
	variational_b::Float64

	expectation_trace::Vector{Float64}

	var_name::String
	function GammaParameter(a::Float64, b::Float64; variable_name = "GammaParameter")
		gp = new()
		gp.prior_a = a
		gp.prior_b = b
		gp.prior_distribution = Distributions.Gamma(a, 1. / b)

		gp.variational_a = a
		gp.variational_b = b
		gp.expectation_trace = Vector{Float64}(undef, 0)
		gp.var_name = variable_name
		gp
	end
end

const VectorGammaParameter = Vector{GammaParameter}

mutable struct NormalParameter <: ModelParameter
	prior_mean::Float64
	prior_variance::Float64
	prior_distribution::Distributions.Normal

	variational_mean::Float64
	variational_variance::Float64

	expectation_trace::Vector{Float64}

	var_name::String
	function NormalParameter(m_0::Float64, s_0::Float64; variable_name = "NormalParameter")
		np = new()
		np.prior_mean = m_0
		np.prior_variance = s_0
		np.prior_distribution = Distributions.Normal(m_0, s_0)

		np.variational_mean = rand(np.prior_distribution)
		np.variational_variance = s_0
		np.expectation_trace = Vector{Float64}(undef, 0)
		np.var_name = variable_name
		np
	end
end

mutable struct MvNormalParameter <: ModelParameter
	prior_mean::Vector{Float64}
	prior_covariance::Matrix{Float64}
	prior_distribution::Distributions.MvNormal

	variational_mean::Vector{Float64}
	variational_covariance::Matrix{Float64}

	expectation_trace::Vector{Vector{Float64}}

	var_name::String
	function MvNormalParameter(m_0::Float64, s_0::T, D::Int64; variable_name = "MvNormalParameter") where {T<:Union{Float64, Vector{Float64}}}
		np = new()
		np.prior_mean = m_0.*ones(D)
		np.prior_covariance = s_0 .* Matrix(1.0I, D, D)
		np.prior_distribution = Distributions.MvNormal(np.prior_mean, np.prior_covariance)

		np.variational_mean = rand(np.prior_distribution)
		np.variational_covariance = s_0 .* Matrix(1.0I, D, D)
		np.expectation_trace = Vector{Vector{Float64}}()
		np.var_name = variable_name
		np
	end
end

set_variable_name(p::ModelParameter, name::String) = p.var_name = name

function check_params(p::GammaParameter)
	if p.variational_a < 0
		# @warn "paramater negative:" p.var_name  p.variational_a
		p.variational_a = 1e-7
	end
	if p.variational_b < 0
		# @warn "paramater negative:" p.var_name  p.variational_b
		p.variational_b = 1e-7
	end
end

function check_params(p::NormalParameter)
	if p.variational_variance ≤ 0
		# @warn "variance non-positive" p.var_name p.variational_variance
	end
	p.variational_variance = 1e-7
end

function check_params(p::MvNormalParameter)
	if !issymmetric(p.variational_covariance) || !isposdef(p.variational_covariance)
		p.variational_covariance = Array(Hermitian(p.variational_covariance))
	end
end

expected_value(p::GammaParameter) = p.variational_a/p.variational_b
expected_log_value(p::GammaParameter) = digamma(p.variational_a) - log(p.variational_b)

expected_value(p::Union{NormalParameter, MvNormalParameter}) = p.variational_mean
expected_squared_value(p::NormalParameter) = p.variational_mean^2 + p.variational_variance
expected_squared_value(p::MvNormalParameter) = p.variational_mean*p.variational_mean' + p.variational_covariance

function elbo(gp::GammaParameter)
	check_params(gp)
	q_dist = Distributions.Gamma(gp.variational_a, 1. / gp.variational_b)
	compute_elbo(q_dist, gp)
end

function elbo(np::NormalParameter)
	check_params(np)
	# p_dist = Distributions.Normal(np.prior_mean, np.prior_variance)
	q_dist = Distributions.Normal(np.variational_mean, np.variational_variance)
	compute_elbo(q_dist, np)
end

function elbo(np::MvNormalParameter)
	check_params(np)
	# p_dist = Distributions.MvNormal(np.prior_mean, np.prior_covariance)
	q_dist = Distributions.MvNormal(np.variational_mean, np.variational_covariance)
	compute_elbo(q_dist, np)
end


function compute_elbo(q_dist::Distributions.Distribution, param::ModelParameter)
	expected_param = expected_value(param)
	push!(param.expectation_trace, expected_param)
	Distributions.logpdf(param.prior_distribution, expected_param) - Distributions.logpdf(q_dist, expected_param)
end
