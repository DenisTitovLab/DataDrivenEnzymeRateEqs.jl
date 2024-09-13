# using TestEnv
# TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test, BenchmarkTools

##
#=
TODO: Add tests for the following:

Use Symbolics.jl to test that the rate equation has all possible combinations of
substrates, products, and regulators with corresponding alphas
=#

##
#tests for rate function speed and allocations
substrates = [Symbol(:S, i) for i = 1:rand(1:3)]
products = [Symbol(:P, i) for i = 1:rand(1:3)]
inhibitors = [Symbol(:I, i) for i = 1:rand(1:2)]
regulators = [Symbol(:R, i) for i = 1:rand(1:5)]
Keq = 20_000.0
enzyme_parameters = (;
    substrates = substrates,
    products = products,
    inhibitors = inhibitors,
    regulators = regulators,
    Keq = Keq,
    oligomeric_state = rand(1:4),
    rate_equation_name = :rand_enz_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(Keq))
@test mean(benchmark_result.times) <= 800 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test mean(benchmark_result.times) <= 800 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt = NamedTuple{metab_names}((
    zeros(length(substrates))...,
    rand(length(products))...,
    rand(length(inhibitors))...,
    rand(length(regulators))...,
))
@test rand_enz_rate_equation(metabs_nt, params_nt, enzyme_parameters.Keq) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt = NamedTuple{metab_names}((
    rand(length(substrates))...,
    zeros(length(products))...,
    rand(length(inhibitors))...,
    rand(length(regulators))...,
))
@test rand_enz_rate_equation(metabs_nt, params_nt, Keq) > 0.0

#test Rate = Vmax_a when [Substrates] and [Activators] -> Inf and Vmax_a_rev when [Products] and [Activators] -> Inf
Vmax_a = 1.0
Vmax_i = rand()
params_vec = []
for param_name in propertynames(params_nt)
    if startswith(string(param_name), "K_a")
        push!(params_vec, 1e-3)
    elseif startswith(string(param_name), "alpha_")
        push!(params_vec, rand([0.0, 1.0]))
    else
        push!(params_vec, 1.0)
    end
end
params_nt = NamedTuple{param_names}(params_vec)
metabs_nt = NamedTuple{metab_names}((
    (1e12 .* ones(length(substrates)))...,
    zeros(length(products))...,
    zeros(length(inhibitors))...,
    1e12 .* ones(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax_a, rtol = 1e-2)
global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_params = [Symbol("K_a_", Symbol(:S, i)) for i = 1:length(substrates)]
product_params = [Symbol("K_a_", Symbol(:P, i)) for i = 1:length(products)]
for param in propertynames(params_nt)
    if param ∈ substrate_params
        global substrate_Ks *= params_nt[param]
    elseif param ∈ product_params
        global product_Ks *= params_nt[param]
    end
end
Vmax_a_rev = Vmax_a * product_Ks / (Keq * substrate_Ks)
metabs_nt = NamedTuple{metab_names}((
    (zeros(length(substrates)))...,
    1e12 .* ones(length(products))...,
    zeros(length(inhibitors))...,
    1e12 .* ones(length(regulators))...,
))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_a_rev, rtol = 1e-2)

#test Rate = Vmax_i when [Substrates] and [Allosteric Inhibitors] -> Inf and Vmax_i_rev when [Products] and [Allosteric Inhibitors] -> Inf
Vmax_a = 1.0
Vmax_i = rand()
params_vec = []
for param_name in propertynames(params_nt)
    if startswith(string(param_name), "K_i")
        push!(params_vec, 1e-3)
    elseif param_name == :Vmax_i
        push!(params_vec, Vmax_i)
    elseif startswith(string(param_name), "alpha_")
        push!(params_vec, rand([0.0, 1.0]))
    else
        push!(params_vec, 1.0)
    end
end
params_nt = NamedTuple{param_names}(params_vec)
metabs_nt = NamedTuple{metab_names}((
    (1e12 .* ones(length(substrates)))...,
    zeros(length(products))...,
    zeros(length(inhibitors))...,
    1e12 .* ones(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax_i, rtol = 1e-2)
global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_params = [Symbol("K_i_", Symbol(:S, i)) for i = 1:length(substrates)]
product_params = [Symbol("K_i_", Symbol(:P, i)) for i = 1:length(products)]
for param in propertynames(params_nt)
    if param ∈ substrate_params
        global substrate_Ks *= params_nt[param]
    elseif param ∈ product_params
        global product_Ks *= params_nt[param]
    end
end
Vmax_i_rev = Vmax_i * product_Ks / (Keq * substrate_Ks)
metabs_nt = NamedTuple{metab_names}((
    (zeros(length(substrates)))...,
    1e12 .* ones(length(products))...,
    zeros(length(inhibitors))...,
    1e12 .* ones(length(regulators))...,
))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_i_rev, rtol = 1e-2)

#test Rate = 0 when [Inhibitors] -> Inf
Vmax_a = 1.0
Vmax_i = 1.0
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}((
    rand(length(substrates))...,
    rand(length(products))...,
    1e12 .* ones(length(inhibitors))...,
    rand(length(regulators))...,
))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), 0.0, atol = 1e-10)
