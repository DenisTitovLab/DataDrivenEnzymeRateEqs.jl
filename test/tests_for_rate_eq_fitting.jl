# using TestEnv
# TestEnv.activate()

#test `@derive_general_mwc_rate_eq`, `loss_rate_equation` and `fit_rate_equation` on real PKM2 data
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools

#Test the loss_rate_equation speed using PKM2 enzyme data
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
#Load and process data
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/PKM2_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source = PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig
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
num_alphas = sum([1 for param_name in param_names if occursin("alpha", string(param_name))])
kinetic_params = [[rand() for i = 1:length(param_names)-num_alphas]..., [rand([0, 1]) for i = 1:num_alphas]...]
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation(kinetic_params, pkm2_rate_equation_no_Keq, rate_data_nt, param_names, fig_point_indexes)
@test mean(benchmark_result.times) <= 150_000 #ns
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation($(kinetic_params), pkm2_rate_equation_no_Keq, $(rate_data_nt), $(param_names), $(fig_point_indexes))
@test mean(benchmark_result.times) <= 150_000 #ns


#Test the loss_rate_equation speed using LDH enzyme data
LDH_enzyme = (;
    substrates=[:Pyruvate, :NADH],
    products=[:Lactate, :NAD],
    regulators=[],
    Keq=20_000.0,
    rate_equation_name=:ldh_rate_equation,
)
metab_names, param_names = @derive_general_qssa_rate_eq(LDH_enzyme)
ldh_rate_equation_no_Keq(metabs, p) = ldh_rate_equation(metabs, p, 20000.0)
#Load and process data
LDH_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/LDH_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
LDH_data_for_fit.source = LDH_data_for_fit.Article .* "_" .* LDH_data_for_fit.Fig
data = LDH_data_for_fit
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
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation(kinetic_params, ldh_rate_equation_no_Keq, rate_data_nt, param_names, fig_point_indexes)
@test mean(benchmark_result.times) <= 150_000 #ns
benchmark_result = @benchmark DataDrivenEnzymeRateEqs.loss_rate_equation($(kinetic_params), ldh_rate_equation_no_Keq, $(rate_data_nt), $(param_names), $(fig_point_indexes))
@test mean(benchmark_result.times) <= 150_000 #ns

# Test the ability of `fit_rate_equation` to recover parameters used to generated data for an arbitrary enzyme
#=
TODO: make below to use more complex test_rate_equation with 1-3 S, P and R and
randomly generated parameters and data around K values.
Add an option to fit real Vmax values instead of fixing Vmax=1.0.
=#
test_rate_equation_Keq = 1.0
test_rate_equation(metabs, params) = params.Vmax * (metabs.S / params.K_S - (1 / test_rate_equation_Keq) * metabs.P / params.K_P) / (1 + metabs.S / params.K_S + metabs.P / params.K_P)
param_names = (:Vmax, :K_S, :K_P)
metab_names = (:S, :P)
params = (Vmax=10.0, K_S=1.0, K_P=5.0)
#create DataFrame of simulated data
num_datapoints = 10
num_figures = 8
S_concs = Float64[]
P_concs = Float64[]
sources = String[]
for i in 1:num_figures
    if i < num_figures ÷ 2
        for S in range(0, 10 * params.K_S, num_datapoints)
            push!(S_concs, S)
            push!(P_concs, 0.0)
            push!(sources, "Figure$i")
        end
    else
        for P in range(0, 10 * params.K_P, num_datapoints)
            push!(S_concs, 0.0)
            push!(P_concs, P)
            push!(sources, "Figure$i")
        end
    end
end
data = DataFrame(S=S_concs, P=P_concs, source=sources)
noise_sd = 0.2
data.Rate = [test_rate_equation(row, params) * (1 + noise_sd * randn()) for row in eachrow(data)]
filter!(row -> row.Rate > 0, data)
fit_result = fit_rate_equation(test_rate_equation, data, metab_names, param_names; n_iter=20)

@test isapprox(fit_result.params.K_S, params.K_S, rtol=3 * noise_sd)
@test isapprox(fit_result.params.K_P, params.K_P, rtol=3 * noise_sd)
