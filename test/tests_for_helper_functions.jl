# WORK IN PROGRESS

# using TestEnv
# TestEnv.activate()

##

# using DataDrivenEnzymeRateEqs, Test
# using Symbolics

# PKM2_enzyme = (;
#     substrates = [:PEP, :ADP],
#     products = [:Pyruvate, :ATP],
#     cat1 = [:ATP, :ADP],
#     cat2 = [:PEP, :Pyruvate],
#     reg1 = [:F16BP],
#     reg2 = [:Phenylalanine],
#     Keq = 20_000.0,
#     oligomeric_state = 4,
#     rate_equation_name = :pkm_rate_equation,
# )
# metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)



# vect=[]
# for x in metab_names
#     push!(vect, (@variables $(x))...)
# end
# nt_metab_names = NamedTuple{metab_names}(vect)
# vect=[]
# for x in param_names
#     push!(vect, (@variables $(x))...)
# end
# nt_param_names = NamedTuple{param_names}(vect)
# @variables Keq
# sym_rate_eqn = pkm_rate_equation(nt_metab_names, nt_param_names, Keq)
