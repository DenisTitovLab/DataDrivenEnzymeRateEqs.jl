# using TestEnv
# TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test
using Symbolics

PKM2_enzyme = (;
    substrates = [:PEP, :ADP],
    products = [:Pyruvate, :ATP],
    cat1 = [:PEP, :Pyruvate],
    cat2 = [:ATP, :ADP],
    reg1 = [:F16BP],
    reg2 = [:Phenylalanine],
    Keq = 20_000.0,
    oligomeric_state = 4,
    rate_equation_name = :pkm2_rate_equation
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)

nt_param_removal_code=(L=0, Vmax=1, K_PEP_cat1=0, K_Pyruvate_cat1=1, K_ATP_cat2=1, K_ADP_cat2=1, K_F16BP_reg1=3, K_Phenylalanine_reg2=2, alpha_PEP_ATP=0, alpha_ADP_Pyruvate=0)
sym_equation = display_rate_equation(pkm2_rate_equation, metab_names, param_names)
sym_equation = display_rate_equation(pkm2_rate_equation, metab_names, param_names; nt_param_removal_code=nt_param_removal_code)
