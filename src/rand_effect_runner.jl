using DataDrivenEnzymeRateEqs, Test, MixedModels
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
# using BenchmarkTools
include("rate_equation_fitting.jl")


"Loss function used for fitting that calculate log of ratio of rate equation predicting of rate and rate data"
function loss_likelihood_rate_equation(
    params,
    rate_equation::Function,
    rate_data_nt::NamedTuple,
    param_names::Tuple{Symbol,Vararg{Symbol}},
    fig_point_indexes::Vector{Vector{Int}};
    rescale_params_from_0_10_scale = true,
    nt_param_removal_code = nothing,
)
    if rescale_params_from_0_10_scale
        kinetic_params = param_rescaling(params, param_names)
    else
        !rescale_params_from_0_10_scale
        kinetic_params = params
    end
    if !isnothing(nt_param_removal_code)
        kinetic_params .=
            param_subset_select(kinetic_params, param_names, nt_param_removal_code)
    end

    #precalculate log_pred_vs_data_ratios for all points as it is expensive and reuse it for weights and loss
    kinetic_params_nt =
        NamedTuple{param_names}(ntuple(i -> kinetic_params[i], Val(length(kinetic_params))))
    log_pred_vs_data_ratios =
        log_ratio_predict_vs_data(rate_equation, rate_data_nt, kinetic_params_nt)

    df = DataFrame(figure=rate_data_nt.fig_num, log_ratios=log_pred_vs_data_ratios)
    model_formula = @formula(-log_ratios ~ (1|figure)) # - since it suppose to be log(actual/pred)
    model = fit!(LinearMixedModel(model_formula, df))
    
    return -loglikelihood(model) 
end



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
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "/home/ec2-user/code/DataDrivenEnzymeRateEqs.jl/test/Data_for_tests/PKM2_data.csv"), DataFrame)
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

loss_liklihood_val = loss_likelihood_rate_equation(kinetic_params, pkm2_rate_equation_no_Keq, rate_data_nt, param_names, fig_point_indexes)
may = 5
# benchmark_result = @benchmark loss_likelihood_rate_equation(kinetic_params, pkm2_rate_equation_no_Keq, rate_data_nt, param_names, fig_point_indexes)
# @test mean(benchmark_result.times) <= 150_000 #ns
# benchmark_result = @benchmark loss_likelihood_rate_equation($(kinetic_params), pkm2_rate_equation_no_Keq, $(rate_data_nt), $(param_names), $(fig_point_indexes))
# @test mean(benchmark_result.times) <= 150_000 #ns


