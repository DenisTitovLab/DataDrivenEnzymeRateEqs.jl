using TestEnv
TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test, BenchmarkTools

##
substrates = [Symbol(:S, i) for i in 1:rand(1:3)]
products = [Symbol(:P, i) for i in 1:rand(1:3)]
enzyme_parameters = (; substrates=substrates, products=products, Keq=20_000.0, rate_equation_name=:rand_enz_rate_equation)

enzyme_parameters
metab_names, param_names = @derive_general_rapid_equilibrium_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(20000.0))
@test mean(benchmark_result.times) <= 100 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, 20000.0)
@test mean(benchmark_result.times) <= 150 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1.0e-3, ATP=1.0e-3, F16BP=1.0e-3, Phenylalanine=1.0e-3)
@test pkm2_rate_equation(metabs_nt, params_nt, 20000.0) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt =
    (PEP=1.0e-3, ADP=1.0e-3, Pyruvate=0.0, ATP=0.0, F16BP=1.0e-3, Phenylalanine=1.0e-3)
@test pkm2_rate_equation(metabs_nt, params_nt, 20000.0) > 0.0

#test Rate = Vmax_a/i when [Substrates] -> Inf and Vmax_a/i_rev when [Products] -> Inf
Vmax_a = 1.0
Vmax_i = 1.0
metabs_nt =
    (PEP=1e12, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), Vmax_a, rtol=1e-6)
metabs_nt =
    (PEP=1e12, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=0.0, Phenylalanine=1.0e12)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), Vmax_i, rtol=1e-6)
Vmax_a_rev = params_nt.Vmax_a * params_nt.K_a_ATP_cat1 * params_nt.K_a_Pyruvate_cat2 / (20000.0 * params_nt.K_a_PEP_cat2 * params_nt.K_a_ADP_cat1)
Vmax_i_rev = params_nt.Vmax_i * params_nt.K_i_ATP_cat1 * params_nt.K_i_Pyruvate_cat2 / (20000.0 * params_nt.K_i_PEP_cat2 * params_nt.K_i_ADP_cat1)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e24, ATP=1e24, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_a_rev, rtol=1e-6)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e12, ATP=1e12, F16BP=0.0, Phenylalanine=1e12)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_i_rev, rtol=1e-6)
