using TestEnv
TestEnv.activate()

using EnzymeFitting, Test
using BenchmarkTools, CMAEvolutionStrategy, DataFrames, CSV, Statistics

@derive_mwc_rate_eq(substrates = [:PEP, :ADP],
    products = [:Pyruvate, :ATP], reg1 = [:F16BP], reg2 = [:Phenylalanine], Keq = 20_000.0)
rate_equation(metabs, p) = rate_equation(metabs, p, 20000.0)
#Load and process data
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/PKM2_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source = PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

# "Names of parameters. Make sure it matches exactly allocation of p in rate_equation()"
param_names = (
    :L,
    :Vmax_a,
    :Vmax_i,
    :K_a_PEP_cat,
    :K_i_PEP_cat,
    :K_a_ADP_cat,
    :K_i_ADP_cat,
    :K_a_Pyruvate_cat,
    :K_i_Pyruvate_cat,
    :K_a_ATP_cat,
    :K_i_ATP_cat,
    :K_a_F16BP_reg1,
    :K_i_F16BP_reg1,
    :K_a_Phenylalanine_reg2,
    :K_i_Phenylalanine_reg2,
    :alpha_PEP_ATP,
    :alpha_ADP_Pyruvate,
)
# "Names of PKM substrate, products and regulators. Ensure that columns in data.csv file have the same exact spelling for metabolites"
metab_names = (:PEP, :ADP, :Pyruvate, :ATP, :F16BP, :Phenylalanine)

data = PKM2_data_for_fit
data.fig_num = vcat(
    [
        i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for i = 1:length(unique(data.source))
    ]...,
)
# Convert DF to NamedTuple for better type stability / speed
rate_data_nt = Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])

# Make a vector containing indexes of points corresponding to each figure
fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]
kinetic_params = [rand() for i = 1:length(param_names)]
benchmark_result = @benchmark EnzymeFitting.loss_rate_equation(kinetic_params, rate_equation, rate_data_nt, param_names, fig_point_indexes)
@test mean(benchmark_result.times) <= 100_000 #ns
benchmark_result = @benchmark EnzymeFitting.loss_rate_equation($(kinetic_params), rate_equation, $(rate_data_nt), $(param_names), $(fig_point_indexes))
@test mean(benchmark_result.times) <= 100_000 #ns

fit_result = @time "fit_rate_function() on PKM2 data" fit_rate_equation(rate_equation, data, metab_names, param_names; n_iter=20)
@test isapprox(fit_result.loss, 0.08946088323758938, rtol=1e-3)
