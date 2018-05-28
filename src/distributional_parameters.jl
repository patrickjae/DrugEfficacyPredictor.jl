abstract type ModelParameter end

mutable struct GammaParameter <: ModelParameter
	prior_a::Float64
	prior_b::Float64

	variational_a::Float64
	variational_b::Float64

	expectation_trace::Vector{Float64}

	function GammaParameter(a::Float64, b::Float64)
		gp = new()
		gp.prior_a = a
		gp.prior_b = b
		gp.variational_a = a
		gp.variational_b = b
		gp.expectation_trace = Float64[]
		gp
	end
end

const VectorGammaParameter = Vector{GammaParameter}

mutable struct NormalParameter <: ModelParameter
	prior_mean::Float64
	prior_variance::Float64

	variational_mean::Float64
	variational_variance::Float64

	mean_trace::Vector{Float64}

	function NormalParameter(m_0::Float64, s_0::Float64)
		np = new()
		np.prior_mean = m_0
		np.prior_variance = s_0
		np.variational_mean = m_0
		np.variational_variance = s_0
		np
	end
end

mutable struct MvNormalParameter <: ModelParameter
	prior_mean::Vector{Float64}
	prior_covariance::Matrix{Float64}

	variational_mean::Vector{Float64}
	variational_covariance::Matrix{Float64}

	mean_trace::Matrix{Float64}
	function MvNormalParameter(m_0::Float64, s_0::T, D::Int64) where {T<:Union{Float64, Vector{Float64}}}
		np = new()
		np.prior_mean = m_0.*ones(D)
		np.prior_covariance = s_0 .* eye(D)
		np.variational_mean = m_0 .* ones(D)
		np.variational_covariance = s_0 .* eye(D)
		np
	end
end

expected_value(p::GammaParameter) = p.variational_a/p.variational_b
expected_log_value(p::GammaParameter) = digamma(p.variational_a) - log(p.variational_b)

expected_value(p::Union{NormalParameter, MvNormalParameter}) = p.variational_mean
expected_squared_value(p::NormalParameter) = p.variational_mean^2 + p.variational_variance
expected_squared_value(p::MvNormalParameter) = p.variational_mean*p.variational_mean' + p.variational_covariance

function elbo(gp::GammaParameter)
	p_dist = Distributions.Gamma(gp.prior_a, 1./gp.prior_b)
	q_dist = Distributions.Gamma(gp.variational_a, 1./gp.variational_b)
	compute_elbo(p_dist, q_dist, gp)
end

function elbo(np::NormalParameter)
	p_dist = Distributions.Normal(np.prior_mean, np.prior_variance)
	q_dist = Distributions.Normal(np.variational_mean, np.variational_variance)
	compute_elbo(p_dist, q_dist, np)
end

function elbo(np::MvNormalParameter)
	p_dist = Distributions.MvNormal(np.prior_mean, np.prior_covariance)
	q_dist = Distributions.MvNormal(np.variational_mean, np.variational_covariance)
	compute_elbo(p_dist, q_dist, np)
end


function compute_elbo(p_dist::Distributions.Distribution, q_dist::Distributions.Distribution, param::ModelParameter)
	expected_param = expected_value(param)
	Distributions.logpdf(p_dist, expected_param) - Distributions.logpdf(q_dist, expected_param)
end
