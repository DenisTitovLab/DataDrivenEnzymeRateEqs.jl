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
substrates = [Symbol(:S, i) for i in 1:rand(1:3)]
products = [Symbol(:P, i) for i in 1:rand(1:3)]
regulators = [Symbol(:R, i) for i in 1:rand(0:2)]
Keq = 20_000.0
enzyme_parameters = (; substrates=substrates, products=products, regulators=regulators, Keq=Keq, rate_equation_name=:rand_enz_rate_equation)
metab_names, param_names = @derive_general_qssa_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(Keq))
@test mean(benchmark_result.times) <= 1000 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test mean(benchmark_result.times) <= 1200 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt = NamedTuple{metab_names}((zeros(length(substrates))..., rand(length(products))..., rand(length(regulators))...))
@test rand_enz_rate_equation(metabs_nt, params_nt, enzyme_parameters.Keq) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt = NamedTuple{metab_names}((rand(length(substrates))..., zeros(length(products))..., rand(length(regulators))...))
@test rand_enz_rate_equation(metabs_nt, params_nt, Keq) > 0.0

#test Rate = Vmax when [Substrates] -> Inf and Vmax_rev when [Products] -> Inf
Vmax = 1.0
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(((1e6 .* ones(length(substrates)))..., zeros(length(products))..., zeros(length(regulators))...))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax, rtol=1e-2)
global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_Ks_names = Symbol("K", (Symbol("_", :S, i) for i in 1:length(substrates))...)
product_Ks_names = Symbol("K", (Symbol("_", :P, i) for i in 1:length(products))...)
Vmax_rev = Vmax * (params_nt[product_Ks_names])^length(products) / (Keq * (params_nt[substrate_Ks_names])^length(substrates) )
metabs_nt = NamedTuple{metab_names}(((zeros(length(substrates)))..., 1e6 .* ones(length(products))..., zeros(length(regulators))...))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_rev, rtol=1e-2)
