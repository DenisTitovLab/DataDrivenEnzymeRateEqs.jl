# using TestEnv
# TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test, BenchmarkTools

##
substrates = [Symbol(:S, i) for i in 1:rand(1:3)]
products = [Symbol(:P, i) for i in 1:rand(1:3)]
Keq = 20_000.0
enzyme_parameters = (; substrates=substrates, products=products, Keq=Keq, rate_equation_name=:rand_enz_rate_equation)

enzyme_parameters
metab_names, param_names = @derive_general_rapid_equilibrium_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(Keq))
@test mean(benchmark_result.times) <= 100 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test mean(benchmark_result.times) <= 150 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt = NamedTuple{metab_names}((zeros(length(substrates))..., rand(length(products))...))
@test rand_enz_rate_equation(metabs_nt, params_nt, Keq) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt = NamedTuple{metab_names}((rand(length(substrates))..., zeros(length(products))...))
@test rand_enz_rate_equation(metabs_nt, params_nt, Keq) > 0.0

#test Rate = Vmax when [Substrates] -> Inf and Vmax_rev when [Products] -> Inf
Vmax = 1.0
metabs_nt = NamedTuple{metab_names}(((1e12 .* rand(length(substrates)))..., zeros(length(products))...))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax, rtol=1e-6)

propertynames(params_nt)[3]

global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_params = [Symbol("K_", Symbol(:S, i)) for i in 1:length(substrates)]
product_params = [Symbol("K_", Symbol(:P, i)) for i in 1:length(products)]
for param in propertynames(params_nt)
    if param ∈ substrate_params
        global substrate_Ks *= params_nt[param]
    elseif param ∈ product_params
        global product_Ks *= params_nt[param]
    end
end

Vmax_rev = Vmax * product_Ks / (Keq * substrate_Ks)
metabs_nt = NamedTuple{metab_names}(((zeros(length(substrates)))..., 1e12 .* rand(length(products))...))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_rev, rtol=1e-6)
