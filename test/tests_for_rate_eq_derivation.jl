using EnzymeFitting, Test, BenchmarkTools

@derive_mwc_rate_eq(substrates = [:PEP, :ADP],
    products = [:Pyruvate, :ATP], reg1 = [:F16BP], reg2 = [:Phenylalanine], Keq = 20_000.0)
#test `@derive_mwc_rate_eq` generated `rate_equation::Function`
@test rate_equation isa Function
params_nt = (
    L=10.0,
    Vmax_a=1.0,
    Vmax_i=1.0,
    K_a_PEP_cat=1e-3,
    K_i_PEP_cat=100e-3,
    K_a_ADP_cat=1e-3,
    K_i_ADP_cat=1e-3,
    K_a_Pyruvate_cat=1e-3,
    K_i_Pyruvate_cat=1e-3,
    K_a_ATP_cat=1e-3,
    K_i_ATP_cat=1e-3,
    K_a_F16BP_reg1=1e-3,
    K_i_F16BP_reg1=100e-3,
    K_a_Phenylalanine_reg2=100e-3,
    K_i_Phenylalanine_reg2=1e-3,
    alpha_PEP_ATP=1.0,
    alpha_ADP_Pyruvate=1.0,
)
metabs_nt =
    (PEP=1.0e-3, ADP=1.0e-3, Pyruvate=1.0e-3, ATP=1.0e-3, F16BP=1.0e-3, Phenylalanine=1.0e-3)

#test rate_equation for speed and allocations
benchmark_result = @benchmark rate_equation($(metabs_nt), $(params_nt), $(20000.0))
@test mean(benchmark_result.times) <= 100 #ns
@test benchmark_result.allocs == 0

benchmark_result = @benchmark rate_equation(metabs_nt, params_nt, 20000.0)
@test mean(benchmark_result.times) <= 200 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1.0e-3, ATP=1.0e-3, F16BP=1.0e-3, Phenylalanine=1.0e-3)
@test rate_equation(metabs_nt, params_nt, 20000.0) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt =
    (PEP=1.0e-3, ADP=1.0e-3, Pyruvate=0.0, ATP=0.0, F16BP=1.0e-3, Phenylalanine=1.0e-3)
@test rate_equation(metabs_nt, params_nt, 20000.0) > 0.0

#test Rate = Vmax_a/i when [Substrates] -> Inf and Vmax_a/i_rev when [Products] -> Inf
Vmax_a = 1.0
Vmax_i = 1.0
metabs_nt =
    (PEP=1e12, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(rate_equation(metabs_nt, params_nt, 20000.0), Vmax_a, rtol = 1e-6)
metabs_nt =
    (PEP=1e12, ADP=1e12, Pyruvate=0.0, ATP=0.0, F16BP=0.0, Phenylalanine=1.0e12)
@test isapprox(rate_equation(metabs_nt, params_nt, 20000.0), Vmax_i, rtol = 1e-6)
Vmax_a_rev = params_nt.Vmax_a * params_nt.K_a_ATP_cat * params_nt.K_a_Pyruvate_cat / (20000.0 * params_nt.K_a_PEP_cat * params_nt.K_a_ADP_cat)
Vmax_i_rev = params_nt.Vmax_i * params_nt.K_i_ATP_cat * params_nt.K_i_Pyruvate_cat / (20000.0 * params_nt.K_i_PEP_cat * params_nt.K_i_ADP_cat)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e24, ATP=1e24, F16BP=1.0e12, Phenylalanine=0.0)
@test isapprox(rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_a_rev, rtol = 1e-6)
metabs_nt =
    (PEP=0.0, ADP=0.0, Pyruvate=1e12, ATP=1e12, F16BP=0.0, Phenylalanine=1e12)
@test isapprox(rate_equation(metabs_nt, params_nt, 20000.0), -Vmax_i_rev, rtol = 1e-6)
