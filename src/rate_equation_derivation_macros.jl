#=
CODE FOR RATE EQUATION DERIVATION
=#

function general_mwc_rate_equation(
        S1::T,
        S2::T,
        P1::T,
        P2::T,
        R1_reg1::T,
        R2_reg1::T,
        R3_reg1::T,
        R1_reg2::T,
        R2_reg2::T,
        R3_reg2::T,
        R1_reg3::T,
        R2_reg3::T,
        R3_reg3::T,
        L::T,
        Vmax_a::T,
        Vmax_i::T,
        K_a_S1_cat::T,
        K_i_S1_cat::T,
        K_a_S2_cat::T,
        K_i_S2_cat::T,
        K_a_P1_cat::T,
        K_i_P1_cat::T,
        K_a_P2_cat::T,
        K_i_P2_cat::T,
        K_a_R1_reg1::T,
        K_i_R1_reg1::T,
        K_a_R2_reg1::T,
        K_i_R2_reg1::T,
        K_a_R3_reg1::T,
        K_i_R3_reg1::T,
        K_a_R1_reg2::T,
        K_i_R1_reg2::T,
        K_a_R2_reg2::T,
        K_i_R2_reg2::T,
        K_a_R3_reg2::T,
        K_i_R3_reg2::T,
        K_a_R1_reg3::T,
        K_i_R1_reg3::T,
        K_a_R2_reg3::T,
        K_i_R2_reg3::T,
        K_a_R3_reg3::T,
        K_i_R3_reg3::T,
        alpha_S1_P2::T,
        alpha_S2_P1::T,
        Keq::T
) where {T <: Float64}
    Vmax_a = 1.0
    Vmax_a_rev = ifelse(
        (K_a_P2_cat * K_a_P1_cat) != Inf,
        Vmax_a * K_a_P2_cat * K_a_P1_cat / (Keq * K_a_S1_cat * K_a_S2_cat),
        0.0
    )
    Vmax_i_rev = ifelse(
        (K_i_P2_cat * K_i_P1_cat) != Inf,
        Vmax_i * K_i_P2_cat * K_i_P1_cat / (Keq * K_i_S1_cat * K_i_S2_cat),
        0.0
    )

    Z_a_cat = (
        1 +
        (S1 / K_a_S1_cat) +
        (S2 / K_a_S2_cat) +
        (P1 / K_a_P1_cat) +
        (P2 / K_a_P2_cat) +
        (S1 / K_a_S1_cat) * (S2 / K_a_S2_cat) +
        (P1 / K_a_P1_cat) * (P2 / K_a_P2_cat) +
        alpha_S1_P2 * (S1 / K_a_S1_cat) * (P2 / K_a_P2_cat) +
        alpha_S2_P1 * (P1 / K_a_P1_cat) * (S2 / K_a_S2_cat)
    )
    Z_i_cat = (
        1 +
        (S1 / K_i_S1_cat) +
        (S2 / K_i_S2_cat) +
        (P1 / K_i_P1_cat) +
        (P2 / K_i_P2_cat) +
        (S1 / K_i_S1_cat) * (S2 / K_i_S2_cat) +
        (P1 / K_i_P1_cat) * (P2 / K_i_P2_cat) +
        alpha_S1_P2 * (S1 / K_i_S1_cat) * (P2 / K_i_P2_cat) +
        alpha_S2_P1 * (P1 / K_i_P1_cat) * (S2 / K_i_S2_cat)
    )
    Z_a_reg = (
        (1 + R1_reg1 / K_a_R1_reg1 + R2_reg1 / K_a_R2_reg1 + R3_reg1 / K_a_R3_reg1) *
        (1 + R1_reg2 / K_a_R1_reg2 + R2_reg2 / K_a_R2_reg2 + R3_reg2 / K_a_R3_reg2) *
        (1 + R1_reg3 / K_a_R1_reg3 + R2_reg3 / K_a_R2_reg3 + R3_reg3 / K_a_R3_reg3)
    )
    Z_i_reg = (
        (1 + R1_reg1 / K_i_R1_reg1 + R2_reg1 / K_i_R2_reg1 + R3_reg1 / K_i_R3_reg1) *
        (1 + R1_reg2 / K_i_R1_reg2 + R2_reg2 / K_i_R2_reg2 + R3_reg2 / K_i_R3_reg2) *
        (1 + R1_reg3 / K_i_R1_reg3 + R2_reg3 / K_i_R2_reg3 + R3_reg3 / K_i_R3_reg3)
    )

    Rate = (
        (
            Vmax_a * (S1 / K_a_S1_cat) * (S2 / K_a_S2_cat) -
            Vmax_a_rev * (P1 / K_a_P1_cat) * (P2 / K_a_P2_cat)
        ) *
        (Z_a_cat^3) *
        (Z_a_reg^4) +
        L *
        (
            Vmax_i * (S1 / K_i_S1_cat) * (S2 / K_i_S2_cat) -
            Vmax_i_rev * (P1 / K_i_P1_cat) * (P2 / K_i_P2_cat)
        ) *
        (Z_i_cat^3) *
        (Z_i_reg^4)
    ) / ((Z_a_cat^4) * (Z_a_reg^4) + L * (Z_i_cat^4) * (Z_i_reg^4))

    return Rate
end

# propertynames(lowered_code)
# Base.method_argnames(methods(general_mwc_rate_equation)[1])[2:end]

macro derive_rate_eq(metabs_and_regulators_kwargs...)
    expected_kwargs = [:substrates, :products, :reg1, :reg2, :reg3, :Keq]
    enz = NamedTuple()
    for expr in metabs_and_regulators_kwargs
        expr.args[1] ∈ expected_kwargs || error("invalid keyword: ", expr.args[1])
        enz = merge(enz, (; expr.args[1] => eval(expr)))
    end
    enz1 = NamedTuple()
    for field in propertynames(enz)
        if field == :substrates
            for (i, metab_name) in enumerate(enz[field])
                enz1 = merge(enz1, (; Symbol(:S, i) => enz[field][i]))
            end
        elseif field == :products
            for (i, metab_name) in enumerate(enz[field])
                enz1 = merge(enz1, (; Symbol(:P, i) => enz[field][i]))
            end
        elseif field == :reg1
            for (i, metab_name) in enumerate(enz[field])
                enz1 = merge(enz1, (; Symbol(:R, i, "_", field) => enz[field][i]))
            end
        elseif field == :reg2
            for (i, metab_name) in enumerate(enz[field])
                enz1 = merge(enz1, (; Symbol(:R, i, "_", field) => enz[field][i]))
            end
        elseif field == :reg3
            for (i, metab_name) in enumerate(enz[field])
                enz1 = merge(enz1, (; Symbol(:R, i, "_", field) => enz[field][i]))
            end
        end
    end
    enz_nt_keys = [:S1, :S2, :P1, :P2, :R1_reg1, :R2_reg1, :R3_reg1, :R1_reg2,
        :R2_reg2, :R3_reg2, :R1_reg3, :R2_reg3, :R3_reg3]
    enz2 = NamedTuple()
    for key in enz_nt_keys
        if haskey(enz1, key)
            enz2 = merge(enz2, (; key => enz1[key]))
        else
            enz2 = merge(enz2, (; key => nothing))
        end
    end
    println(enz2)
    qualified_name = esc(GlobalRef(Main, :rate_equation))
    return :(function $(qualified_name)(metabs, params, Keq)
        general_mwc_rate_equation(
            $(enz2.S1 isa Symbol) ? metabs.$(enz2.S1) : 1.0,
            $(enz2.S2 isa Symbol) ? metabs.$(enz2.S2) : 1.0,
            $(enz2.P1 isa Symbol) ? metabs.$(enz2.P1) : 1.0,
            $(enz2.P2 isa Symbol) ? metabs.$(enz2.P2) : 1.0,
            $(enz2.R1_reg1 isa Symbol) ? metabs.$(enz2.R1_reg1) : 0.0,
            $(enz2.R2_reg1 isa Symbol) ? metabs.$(enz2.R2_reg1) : 0.0,
            $(enz2.R3_reg1 isa Symbol) ? metabs.$(enz2.R3_reg1) : 0.0,
            $(enz2.R1_reg2 isa Symbol) ? metabs.$(enz2.R1_reg2) : 0.0,
            $(enz2.R2_reg2 isa Symbol) ? metabs.$(enz2.R2_reg2) : 0.0,
            $(enz2.R3_reg2 isa Symbol) ? metabs.$(enz2.R3_reg2) : 0.0,
            $(enz2.R1_reg3 isa Symbol) ? metabs.$(enz2.R1_reg3) : 0.0,
            $(enz2.R2_reg3 isa Symbol) ? metabs.$(enz2.R2_reg3) : 0.0,
            $(enz2.R3_reg3 isa Symbol) ? metabs.$(enz2.R3_reg3) : 0.0,
            params.L,
            params.Vmax_a,
            params.Vmax_i,
            $(enz2.S1 isa Symbol) ? params.$(Symbol("K_a_", enz2.S1, "_cat")) : 1.0,
            $(enz2.S1 isa Symbol) ? params.$(Symbol("K_i_", enz2.S1, "_cat")) : 1.0,
            $(enz2.S2 isa Symbol) ? params.$(Symbol("K_a_", enz2.S2, "_cat")) : 1.0,
            $(enz2.S2 isa Symbol) ? params.$(Symbol("K_i_", enz2.S2, "_cat")) : 1.0,
            $(enz2.P1 isa Symbol) ? params.$(Symbol("K_a_", enz2.P1, "_cat")) : 1.0,
            $(enz2.P1 isa Symbol) ? params.$(Symbol("K_i_", enz2.P1, "_cat")) : 1.0,
            $(enz2.P2 isa Symbol) ? params.$(Symbol("K_a_", enz2.P2, "_cat")) : 1.0,
            $(enz2.P2 isa Symbol) ? params.$(Symbol("K_i_", enz2.P2, "_cat")) : 1.0,
            $(enz2.R1_reg1 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R1_reg1, "_reg1")) : Inf,
            $(enz2.R1_reg1 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R1_reg1, "_reg1")) : Inf,
            $(enz2.R2_reg1 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R2_reg1, "_reg1")) : Inf,
            $(enz2.R2_reg1 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R2_reg1, "_reg1")) : Inf,
            $(enz2.R3_reg1 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R3_reg1, "_reg1")) : Inf,
            $(enz2.R3_reg1 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R3_reg1, "_reg1")) : Inf,
            $(enz2.R1_reg2 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R1_reg2, "_reg2")) : Inf,
            $(enz2.R1_reg2 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R1_reg2, "_reg2")) : Inf,
            $(enz2.R2_reg2 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R2_reg2, "_reg2")) : Inf,
            $(enz2.R2_reg2 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R2_reg2, "_reg2")) : Inf,
            $(enz2.R3_reg2 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R3_reg2, "_reg2")) : Inf,
            $(enz2.R3_reg2 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R3_reg2, "_reg2")) : Inf,
            $(enz2.R1_reg3 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R1_reg3, "_reg3")) : Inf,
            $(enz2.R1_reg3 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R1_reg3, "_reg3")) : Inf,
            $(enz2.R2_reg3 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R2_reg3, "_reg3")) : Inf,
            $(enz2.R2_reg3 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R2_reg3, "_reg3")) : Inf,
            $(enz2.R3_reg3 isa Symbol) ?
            params.$(Symbol("K_a_", enz2.R3_reg3, "_reg3")) : Inf,
            $(enz2.R3_reg3 isa Symbol) ?
            params.$(Symbol("K_i_", enz2.R3_reg3, "_reg3")) : Inf,
            params.alpha_PEP_ATP,
            params.alpha_ADP_Pyruvate,
            Keq
        )
    end
    )
end


# metabs_nt =
#     (PEP = 1.0e-3, ADP = 1.0e-3, Pyruvate = 1.0e-3, ATP = 1.0e-3, F16BP = 1.0e-3, Phenylalanine = 1.0e-3)

# params_nt = (
#     L = 1.0,
#     Vmax_a = 1.0,
#     Vmax_i = 1.0,
#     K_a_PEP_cat = 1e-3,
#     K_i_PEP_cat = 100e-3,
#     K_a_ADP_cat = 1e-3,
#     K_i_ADP_cat = 1e-3,
#     K_a_Pyruvate_cat = 1e-3,
#     K_i_Pyruvate_cat = 1e-3,
#     K_a_ATP_cat = 1e-3,
#     K_i_ATP_cat = 1e-3,
#     K_a_F16BP_reg1 = 1e-3,
#     K_i_F16BP_reg1 = 100e-3,
#     K_a_Phenylalanine_reg2 = 100e-3,
#     K_i_Phenylalanine_reg2 = 1e-3,
#     alpha_PEP_ATP = 1.0,
#     alpha_ADP_Pyruvate = 1.0,
# )
# @derive_rate_eq(substrates=[:PEP, :ADP],
#     products=[:Pyruvate, :ATP], reg1=[:F16BP], reg2=[:Phenylalanine],Keq=20_000.0)
# @code_warntype rate_equation(metabs_nt, params_nt, 20000.0)
# using BenchmarkTools
# @benchmark rate_equation(metabs_nt, params_nt, 20000.0)
# @benchmark rate_equation($(metabs_nt), $(params_nt), 20000.0)
# rate_equation(metabs_nt, params_nt, 20000.0)
