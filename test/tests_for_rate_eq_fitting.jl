# using TestEnv
# TestEnv.activate()

#test `@derive_general_mwc_rate_eq`, `loss_rate_equation` and `fit_rate_equation` on real PKM2 data
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools

PKM2_enzyme = (;
    substrates = [:PEP, :ADP],
    products = [:Pyruvate, :ATP],
    cat1 = [:ATP, :ADP],
    cat2 = [:PEP, :Pyruvate],
    reg1 = [:F16BP],
    reg2 = [:Phenylalanine],
    Keq = 20_000.0,
    oligomeric_state = 4,
    rate_equation_name = :rate_equation,
)
@derive_general_mwc_rate_eq(PKM2_enzyme)



rate_equation(metabs, p) = rate_equation(metabs, p, 20000.0)
#Load and process data
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/PKM2_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source = PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig

param_names = (
    :L,
    :Vmax_a,
    :Vmax_i,
    :K_a_PEP_cat2,
    :K_i_PEP_cat2,
    :K_a_ADP_cat1,
    :K_i_ADP_cat1,
    :K_a_Pyruvate_cat2,
    :K_i_Pyruvate_cat2,
    :K_a_ATP_cat1,
    :K_i_ATP_cat1,
    :K_a_F16BP_reg1,
    :K_i_F16BP_reg1,
    :K_a_Phenylalanine_reg2,
    :K_i_Phenylalanine_reg2,
    :alpha_PEP_ATP,
    :alpha_ADP_Pyruvate,
)
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
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation(kinetic_params, rate_equation, rate_data_nt, param_names, fig_point_indexes)
@test mean(benchmark_result.times) <= 100_000 #ns
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation($(kinetic_params), rate_equation, $(rate_data_nt), $(param_names), $(fig_point_indexes))
@test mean(benchmark_result.times) <= 100_000 #ns

#TODO: make fake data with noise and known params and ensure known params are recovered
fit_result = fit_rate_equation(rate_equation, data, metab_names, param_names; n_iter=20)
@test isapprox(fit_result.train_loss, 0.08946088323758938, rtol=1e-3)
@test fit_result.params isa NamedTuple{param_names}{NTuple{length(param_names),Float64}}

##
#test the ability of `fit_rate_equation` to recover parameters used to generated data for an arbitrary enzyme
#=
TODO: delete PKM2 example above after making below to use more complex test_rate_equation with 1-3 S, P and R and
randomly generated parameters and data around K values.
Add an option to fit real Vmax values instead of fixing Vmax=1.0.
=#
test_rate_equation(metabs, params) = params.Vmax * (metabs.S / params.K_S) / (1 + metabs.S / params.K_S)

param_names = (:Vmax, :K_S)
metab_names = (:S,)
params = (Vmax=10.0, K_S=1.0)
data = DataFrame(S=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
noise_sd = 0.2
data.Rate = [test_rate_equation(row, params) * (1 + noise_sd * randn()) for row in eachrow(data)]
data.source = ["Figure1" for i in 1:nrow(data)]
fit_result = fit_rate_equation(test_rate_equation, data, metab_names, param_names; n_iter=20)
@test isapprox(fit_result.params.K_S, params.K_S, rtol=3 * noise_sd)
