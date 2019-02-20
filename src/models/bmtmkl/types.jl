abstract type BMTMKLModel <: Model end

function ModelConfiguration(
	pm::PredictionModel{BMTMKLModel},
	α::Float64 = 1000.,
	β::Float64 = 1000.,
	μ::Float64 = 1.,
	σ::Float64 = 2.
)
	mc = ModelConfiguration()
	# Gamma params for precision on bias term b
	mc.parameters["α_ɣ"] = α
	mc.parameters["β_ɣ"] = β
	# Gamma params for precision on weights a
	mc.parameters["α_λ"] = α
	mc.parameters["β_λ"] = β
	# Gamma params for precision on outcome y
	mc.parameters["α_ε"] = α
	mc.parameters["β_ε"] = β
	# Gamma params for precision on intermediate results g
	mc.parameters["α_ν"] = α
	mc.parameters["β_ν"] = β
	# Gamma params for precision on kernel weights e
	mc.parameters["α_⍵"] = α
	mc.parameters["β_⍵"] = β

	# normal params for mean on bias
	mc.parameters["μ_b"] = 0.
	mc.parameters["σ_b"] = σ
	# normal params for mean on kernel weights
	mc.parameters["μ_e"] = μ
	mc.parameters["σ_e"] = σ
	# mv normal params for weights
	mc.parameters["μ_a"] = μ
	mc.parameters["Σ_a"] = σ
	# mv normal weights for intermediate results
	mc.parameters["μ_g"] = μ
	mc.parameters["Σ_g"] = σ
	# get number of tasks (drugs)
	mc.parameters["T"] = length(pm.data.drugs)
	mc.parameters["N"] = OrderedDict{Drug, Int64}()
	# NOTE: N (the number of samples per drug) is set in the post_init! method for this type of model
	for drug in collect(keys(pm.data.results))
<<<<<<< HEAD
		mc.parameters["N"][drug] = length(filter((cl) -> !cl.in_test_set && haskey(pm.data.results[drug].outcome_values, cl), collect(values(pm.data.cell_lines))))
=======
		mc.parameters["N"][drug] = length(filter((cl) -> !cl.in_test_set && haskey(pm.data.results[drug].outcome_values, cl), collect(values(pm.data.cell_lines)))
>>>>>>> e6c64d1f42787631cd1c74ed639f565ccc3912bc
	end
	# recompute overall number of kernels (+4 accomodates for combined kernels)
	mc.parameters["K"] = length(pm.precomputations["base_kernels"]) + sum(map(length, collect(values(pm.precomputations["pathway_specific_kernels"])))) + 4
	mc
end

mutable struct BMTMKLModelParameters <: PredictionModelParameters
	# precision parameters, T-dimensional, one per drug/task
	# ɣ::IsometricPrecisionParameter #precision to b <-- could be a matrix introducing covariance between drugs
	ɣ::Vector{GammaParameter}
	# drug specific scalars, T-dimensional
	b::Vector{NormalParameter}

	# λ::Vector{IsometricPrecisionParameter} #precision to a <-- introduce covariance between cell lines
	λ::Vector{Vector{GammaParameter}}
	# drug specific vectors, TxN_t dimensional
	# a::Vector{Vector{NormalParameter}}
	a::Vector{MvNormalParameter}

	# ε::IsometricPrecisionParameter #precision to y <-- between drugs
	ε::Vector{GammaParameter}

	# ν::IsometricPrecisionParameter #precision to g <-- between drugs and views?
	ν::Vector{GammaParameter}
	# intermediate results, N_t dimensional, one per drug and view
	G::Matrix{MvNormalParameter}

	# precision parameter, K-dimensional, one per view
	# ⍵::IsometricPrecisionParameter
	⍵::Vector{GammaParameter}
	# view specific scalars, K-dimensional
	e::Vector{NormalParameter}

	# general model params
	T::Int64
	K::Int64
	N::OrderedDict{Drug, Int64}

	function BMTMKLModelParameters(mc::ModelConfiguration)
		p = new()
		p.T = mc.parameters["T"]
		p.N = mc.parameters["N"]
		p.K = mc.parameters["K"]
		p.ɣ = VectorGammaParameter(undef, p.T)
		p.b = Vector{NormalParameter}(undef, p.T)
		for t in 1:p.T
			p.ɣ[t] = GammaParameter(mc.parameters["α_ɣ"], mc.parameters["β_ɣ"])
			p.b[t] = NormalParameter(mc.parameters["μ_b"], mc.parameters["σ_b"])
			set_variable_name(p.ɣ[t], "ɣ[$t]")
			set_variable_name(p.b[t], "b[$t]")
		end

		p.λ = Vector{VectorGammaParameter}(undef, p.T)
		# p.a = Vector{Vector{NormalParameter}}(T)
		p.a = Vector{MvNormalParameter}(undef, p.T)
		for (t, drug) in enumerate(keys(p.N))
			p.λ[t] = VectorGammaParameter(undef, p.N[drug])
			p.a[t] = MvNormalParameter(mc.parameters["μ_a"], mc.parameters["Σ_a"], p.N[drug])
			# p.a[t] = Vector{NormalParameter}(N[t])
			set_variable_name(p.a[t], "a[$t]")
			for n in 1:p.N[drug]
				p.λ[t][n] = GammaParameter(mc.parameters["α_λ"], mc.parameters["β_λ"])
				set_variable_name(p.ɣ[t], "λ[$t][$n]")
				# p.a[t][n] = NormalParameter(m_0, s_0)
			end
		end

		p.ε = VectorGammaParameter(undef, p.T)
		for t in 1:p.T
			p.ε[t] = GammaParameter(mc.parameters["α_ε"], mc.parameters["β_ε"])
			set_variable_name(p.ε[t], "ε[$t]")
		end

		p.ν = VectorGammaParameter(undef, p.T)
		for t in 1:p.T
			p.ν[t] = GammaParameter(mc.parameters["α_ν"], mc.parameters["β_ν"])
			set_variable_name(p.ν[t], "ν[$t]")
		end
		p.G = Matrix{MvNormalParameter}(undef, p.T, p.K)
		for (t, drug) in enumerate(keys(p.N)), k in 1:p.K
			p.G[t,k] = MvNormalParameter(mc.parameters["μ_g"], mc.parameters["Σ_g"], p.N[drug])
			set_variable_name(p.G[t,k], "G[$t,$k]")
		end

		p.⍵ = VectorGammaParameter(undef, p.K)
		p.e = Vector{NormalParameter}(undef, p.K)
		for k in 1:p.K
			p.⍵[k] = GammaParameter(mc.parameters["α_⍵"], mc.parameters["β_⍵"])
			p.e[k] = NormalParameter(mc.parameters["μ_e"], mc.parameters["σ_e"])
			set_variable_name(p.⍵[k], "⍵[$k]")
			set_variable_name(p.e[k], "e[$k]")
		end
		p
	end

end
