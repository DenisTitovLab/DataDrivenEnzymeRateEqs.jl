using Pkg
package_path = "/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl"
Pkg.activate(package_path)
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools
# include("rate_equation_selection.jl")

file_path = joinpath(package_path, "test/Data_for_tests/PKM2_data.csv")
data = CSV.read(file_path, DataFrame)

# enzyme_parameters = (; 
# substrates=[:PEP,:ADP], 
# products=[:Pyruvate, :ATP], 
# cat1=[:PEP, :Pyruvate],
# cat2 = [:ADP, :ATP], 
# reg1=[:F16BP], reg2=[:Phenylalanine], 
# Keq=20_000, oligomeric_state=4,
# rate_equation_name=:derived_rate_equation)

PKM2_enzyme = (;
    substrates=[:PEP, :ADP],
    products=[:Pyruvate, :ATP],
    regulators=[:F16BP, :Phenylalanine],
    Keq=20_000.0,
    oligomeric_state=4,
    rate_equation_name=:pkm2_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)
pkm2_rate_equation_no_Keq(metabs, p) = pkm2_rate_equation(metabs, p, 20000.0)

# metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
# derived_rate_equation_no_Keq(nt_metabs, nt_params) = derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(pkm2_rate_equation_no_Keq,
 data,
  metab_names, 
  param_names, 
  (7, 15), 
  true;
  n_reps_opt=1, # n repeats optimization
  maxiter_opt=5,# n iteration opt algorithm
  model_selection_method = "current_subsets_filtering",
  p_val_threshold = .4,   # pval threshould for choosing best n params
  train_loss = "likelihood",
  test_loss = "sse"
  )




# data_gen_rate_equation_Keq = 1.0
# mwc_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) = (1 / params.K_a_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P)
# mwc_alternative_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) = (1 / params.K_a_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P + metabs.S * metabs.P / (params.K_a_S * params.K_a_P))

# data_gen_param_names = (:Vmax_a, :K_a_S, :K_a_P)
# metab_names = (:S, :P)
# params = (Vmax=10.0, K_a_S=1e-3, K_a_P=5e-3)
# #create DataFrame of simulated data
# num_datapoints = 80
# num_figures = 4
# S_concs = Float64[]
# P_concs = Float64[]
# sources = String[]

# for i in 1:num_figures
#     if i < num_figures ÷ 2
#         for S in range(0, rand(1:10) * params.K_a_S, rand(num_datapoints÷2:num_datapoints*2))
#             push!(S_concs, S)
#             push!(P_concs, 0.0)
#             push!(sources, "Figure$i")
#         end
#     else
#         for P in range(0, rand(1:10) * params.K_a_P, rand(num_datapoints÷2:num_datapoints*2))
#             push!(S_concs, 0.0)
#             push!(P_concs, P)
#             push!(sources, "Figure$i")
#         end
#     end
# end
# data = DataFrame(S=S_concs, P=P_concs, source=sources)
# noise_sd = 0.2
# data.Rate = [mwc_data_gen_rate_equation(row, params, data_gen_rate_equation_Keq) * (1 + noise_sd * randn()) for row in eachrow(data)]

# enzyme_parameters = (; substrates=[:S,], products=[:P], regulators=[], Keq=1.0, oligomeric_state=1, rate_equation_name=:mwc_derived_rate_equation)

# metab_names, derived_param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
# mwc_derived_rate_equation_no_Keq(nt_metabs, nt_params) = mwc_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
# selection_result = @time data_driven_rate_equation_selection(
#     mwc_derived_rate_equation_no_Keq, 
#     data, metab_names, derived_param_names, (3, 7), 
#     true,
#     model_selection_method = "current_subsets_filtering",
#     p_val_threshold = .4,   # pval threshould for choosing best n params
#     train_loss = "likelihood",
#     test_loss = "sse"
#     )

