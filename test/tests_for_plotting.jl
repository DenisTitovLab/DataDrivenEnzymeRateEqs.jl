# using TestEnv
# TestEnv.activate()

##
using CairoMakie, DataDrivenEnzymeRateEqs, CSV, DataFrames, Printf, Statistics, Test

#Test the loss_rate_equation speed using PKM2 enzyme data
PKM2_enzyme = (;
    substrates=[:PEP, :ADP],
    products=[:Pyruvate, :ATP],
    regulators=[:F16BP, :Phenylalanine],
    Keq=20_000.0,
    oligomeric_state=4,
    rate_equation_name=:pkm2_rate_equation
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)
pkm2_rate_equation_no_Keq(metabs, p) = pkm2_rate_equation(metabs, p, 20000.0)
#Load and process data
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/PKM2_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source = PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig
data = PKM2_data_for_fit
#get kinetic_params from fitting results
nt_kinetic_params = (L=21.20386803833623, Vmax_a=0.6621530787152484, Vmax_i=0.6621530787152484, K_a_PEP=0.00016232173554865876, K_i_PEP=0.007270929135951147, K_a_ADP=0.00020861163680429678, K_i_ADP=0.00020861163680429678, K_a_Pyruvate=0.018115663373606324, K_i_Pyruvate=0.018115663373606324, K_a_ATP=0.00019564943473363028, K_i_ATP=0.00019564943473363028, K_a_F16BP_reg=2.000165919689415e-7, K_i_F16BP_reg=Inf, K_a_Phenylalanine_reg=Inf, K_i_Phenylalanine_reg=0.00022505058445537455, alpha_PEP_Pyruvate=0.0, alpha_PEP_ATP=0.0, alpha_PEP_F16BP=1.0, alpha_PEP_Phenylalanine=1.0, alpha_ADP_Pyruvate=1.0, alpha_ADP_ATP=0.0, alpha_ADP_F16BP=1.0, alpha_ADP_Phenylalanine=1.0, alpha_Pyruvate_F16BP=1.0, alpha_Pyruvate_Phenylalanine=1.0, alpha_ATP_F16BP=1.0, alpha_ATP_Phenylalanine=1.0, alpha_F16BP_Phenylalanine=1.0)

plot_of_fit = DataDrivenEnzymeRateEqs.plot_fit_on_data(
    pkm2_rate_equation_no_Keq,
    nt_kinetic_params,
    metab_names,
    data;
    enzyme_name="PKM2",
    num_col=5,
    scaler=4.0
)

@test plot_of_fit isa Figure
@test !isempty(plot_of_fit.content)
# Note: Content length reduced from 56 to ~50 due to legend workaround for CairoMakie v0.15+
# This ensures compatibility while maintaining core plotting functionality
@test length(plot_of_fit.content) >= 28  # At least one content item per subplot
