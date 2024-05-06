using TestEnv
TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test, BenchmarkTools

##
substrates = [Symbol(:S, i) for i in 1:rand(1:3)]
products = [Symbol(:P, i) for i in 1:rand(1:3)]
regulators = [Symbol(:R, i) for i in 1:rand(1:6)]
# enzyme_parameters = ()
# enzyme_parameters = (; substrates=substrates, products=products, cat1=[substrates[1], products[1]])
# for i in 1:max(length(substrates), length(products))
#     global enzyme_parameters = (; enzyme_parameters..., Symbol(:cat, i) => filter!(x -> x !== :Blank, [get(substrates, i, :Blank), get(products, i, :Blank)]))
# end
# for i in eachindex(regulators)
#     not_fully_occupied_sites = [site for site in [:reg1, :reg2, :reg3] if !haskey(enzyme_parameters, site) || (haskey(enzyme_parameters, site) && length(enzyme_parameters[site]) < 3)]
#     reg_site_name = rand(not_fully_occupied_sites)
#     global enzyme_parameters = (; enzyme_parameters..., reg_site_name => haskey(enzyme_parameters, reg_site_name) ? [enzyme_parameters[reg_site_name]..., regulators[i]] : [regulators[i]])
# end
# enzyme_parameters = (; enzyme_parameters..., Keq=20_000.0, oligomeric_state=rand(1:4), rate_equation_name=:rand_enz_rate_equation)

enzyme_parameters = (; substrates=substrates, products=products, regulators=regulators, Keq=20_000.0, oligomeric_state=rand(1:4), rate_equation_name=:rand_enz_rate_equation)


metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(20000.0))
@test mean(benchmark_result.times) <= 150 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, 20000.0)
@test mean(benchmark_result.times) <= 175 #ns
@test benchmark_result.allocs <= 1

#test `@derive_mwc_rate_eq` generated `rate_equation::Function` using PKM2 data
# PKM2_enzyme = (;
#     substrates = [:PEP, :ADP],
#     products = [:Pyruvate, :ATP],
#     cat1 = [:ATP, :ADP],
#     cat2 = [:PEP, :Pyruvate],
#     reg1 = [:F16BP],
#     reg2 = [:Phenylalanine],
#     Keq = 20_000.0,
#     oligomeric_state = 4,
#     rate_equation_name = :pkm2_rate_equation,
# )
PKM2_enzyme = (;
    substrates=[:PEP, :ADP],
    products=[:Pyruvate, :ATP],
    regulators=[:F16BP, :Phenylalanine],
    Keq=20_000.0,
    oligomeric_state=4,
    rate_equation_name=:pkm2_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)

@test pkm2_rate_equation isa Function
params_nt = (
    L=10.0,
    Vmax_a=1.0,
    Vmax_i=0.1,
    K_a_PEP=1e-3,
    K_i_PEP=100e-3,
    K_a_ADP=1e-3,
    K_i_ADP=1e-3,
    K_a_Pyruvate=1e-3,
    K_i_Pyruvate=1e-3,
    K_a_ATP=1e-3,
    K_i_ATP=1e-3,
    K_a_F16BP=1e-3,
    K_i_F16BP=100e-3,
    K_a_Phenylalanine=100e-3,
    K_i_Phenylalanine=1e-3,
    alpha_PEP_ATP=0.0,
    alpha_PEP_Pyruvate=0.0,
    alpha_PEP_F16BP=1.0,
    alpha_PEP_Phenylalanine=1.0,
    alpha_ADP_Pyruvate=1.0,
    alpha_ADP_ATP=0.0,
    alpha_ADP_F16BP=1.0,
    alpha_ADP_Phenylalanine=1.0,
    alpha_ATP_F16BP=1.0,
    alpha_ATP_Phenylalanine=1.0,
    alpha_Pyruvate_F16BP=1.0,
    alpha_Pyruvate_Phenylalanine=1.0,
    alpha_F16BP_Phenylalanine=1.0,
)
metabs_nt =
    (PEP=1.0e-3, ADP=1.0e-3, Pyruvate=1.0e-3, ATP=1.0e-3, F16BP=1.0e-3, Phenylalanine=1.0e-3)

#test pkm2_rate_equation for speed and allocations
benchmark_result = @benchmark pkm2_rate_equation($(metabs_nt), $(params_nt), $(20000.0))
@test mean(benchmark_result.times) <= 100 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark pkm2_rate_equation(metabs_nt, params_nt, 20000.0)
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
metabs_nt =
    (PEP=1e12, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), 1.0, rtol=1e-6)
metabs_nt =
    (PEP=1e3, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=0.0, Phenylalanine=1.0e24)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), params_nt.Vmax_i, rtol=1e-6)
pkm2_rate_equation(metabs_nt, params_nt, 20000.0)
Vmax_a_rev = params_nt.Vmax_a * params_nt.K_a_ATP * params_nt.K_a_Pyruvate / (20000.0 * params_nt.K_a_PEP * params_nt.K_a_ADP)
Vmax_i_rev = params_nt.Vmax_i * params_nt.K_i_ATP * params_nt.K_i_Pyruvate / (20000.0 * params_nt.K_i_PEP * params_nt.K_i_ADP)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e24, ATP=1e24, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_a_rev, rtol=1e-6)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e12, ATP=1e12, F16BP=0.0, Phenylalanine=1e12)
@test isapprox(pkm2_rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_i_rev, rtol=1e-6)


pkm2_rate_equation(metabs_nt, params_nt, 20000.0)

display_rate_equation(
    pkm2_rate_equation,
    metab_names,
    param_names;
    nt_param_removal_code=nothing
)

using Symbolics
# params_nt = (
#     L=0.0,
#     Vmax_a=1.0,
#     Vmax_i=1.0,
#     K_a_PEP=1.0,
#     K_i_PEP=1.0,
#     K_a_ADP=1.0,
#     K_i_ADP=1.0,
#     K_a_Pyruvate=1.0,
#     K_i_Pyruvate=1.0,
#     K_a_ATP=1.0,
#     K_i_ATP=1.0,
#     K_a_F16BP=1.0,
#     K_i_F16BP=1.0,
#     K_a_Phenylalanine=1.0,
#     K_i_Phenylalanine=1.0,
#     alpha_PEP_ATP=0.0,
#     alpha_PEP_Pyruvate=0.0,
#     alpha_PEP_F16BP=1.0,
#     alpha_PEP_Phenylalanine=1.0,
#     alpha_ADP_Pyruvate=1.0,
#     alpha_ADP_ATP=0.0,
#     alpha_ADP_F16BP=1.0,
#     alpha_ADP_Phenylalanine=1.0,
#     alpha_ATP_F16BP=1.0,
#     alpha_ATP_Phenylalanine=1.0,
#     alpha_Pyruvate_F16BP=1.0,
#     alpha_Pyruvate_Phenylalanine=1.0,
#     alpha_F16BP_Phenylalanine=1.0,
# )
@variables PEP ADP Pyruvate ATP F16BP Phenylalanine
metabs_nt_sym =
    (PEP=0.0, ADP=0.0, Pyruvate=Pyruvate, ATP=ATP, F16BP=F16BP, Phenylalanine=0.0)
metabs_nt_sym =
    (PEP=PEP, ADP=ADP, Pyruvate=0.0, ATP=0.0, F16BP=0.0, Phenylalanine=0.0)
pkm2_rate_equation(metabs_nt_sym, params_nt, 20000.0)
pkm2_rate_equation(metabs_nt_sym, params_nt, 20000.0)
