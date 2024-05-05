#=
CODE FOR MWC RATE EQUATION DERIVATION
=#
using Distributed

@doc raw"""
    derive_general_mwc_rate_eq(metabs_and_regulators_kwargs...)

Derive a function that calculates the rate of a reaction using the general MWC rate equation given the list of substrates, products, and regulators that bind to specific cat or reg sites.

The general MWC rate equation is given by:

```math
Rate = \\frac{{V_{max}^a \\prod_{i=1}^{n} \\left(\\frac{S_i}{K_{a, i}}\\right) - V_{max}^a_{rev} \\prod_{i=1}^{n} \\left(\\frac{P_i}{K_{a, i}}\\right) \\cdot Z_{a, cat}^{n-1} \\cdot Z_{a, reg}^n + L \\left(V_{max}^i \\prod_{i=1}^{n} \\left(\\frac{S_i}{K_{i, i}}\\right) - V_{max}^i_{rev} \\prod_{i=1}^{n} \\left(\\frac{P_i}{K_{i, i}}\\right)\\right) \\cdot Z_{i, cat}^{n-1} \\cdot Z_{i, reg}^n}}{Z_{a, cat}^n \\cdot Z_{a, reg}^n + L \\cdot Z_{i, cat}^n \\cdot Z_{i, reg}^n}
```

where:
- ``V_{max}^a`` is the maximum rate of the forward reaction
- ``V_{max rev}^a`` is the maximum rate of the reverse reaction
- ``V_{max}^i`` is the maximum rate of the forward reaction
- ``V_{max rev}^i`` is the maximum rate of the reverse reaction
- ``S_i`` is the concentration of the ``i^{th}`` substrate
- ``P_i`` is the concentration of the ``i^{th}`` product
- ``K_{a, i}`` is the Michaelis constant for the ``i^{th}`` substrate
- ``K_{i, i}`` is the Michaelis constant for the ``i^{th}`` product
- ``Z_{a, cat}`` is the allosteric factor for the catalytic site
- ``Z_{i, cat}`` is the allosteric factor for the catalytic site
- ``Z_{a, reg}`` is the allosteric factor for the regulatory site
- ``Z_{i, reg}`` is the allosteric factor for the regulatory site
- ``L`` is the ratio of inactive to active enzyme conformations in the absence of ligands
- ``n`` is the oligomeric state of the enzyme

# Arguments
- `metabs_and_regulators_kwargs...`: keyword arguments that specify the substrates, products, catalytic sites, regulatory sites, and other parameters of the reaction.

# Returns
- A function that calculates the rate of the reaction using the general MWC rate equation
- A tuple of the names of the metabolites and parameters used in the rate equation
"""
macro derive_general_mwc_rate_eq(metabs_and_regulators_kwargs)

    processed_input = getfield(__module__, Symbol(metabs_and_regulators_kwargs))
    expected_input_kwargs =
        [:substrates, :products, :regulators, :Keq, :oligomeric_state, :rate_equation_name]
    for kwarg in keys(processed_input)
        kwarg ∈ expected_input_kwargs || error(
            "invalid keyword: ",
            kwarg,
            ". The only supported keywords are: ",
            expected_input_kwargs,
        )
    end
    @assert 0 < length(processed_input.substrates) <= 3 "At least 1 and no more that 3 substrates are supported"
    @assert 0 < length(processed_input.products) <= 3 "At least 1 and no more that 3 products are supported"
    @assert 0 < length(processed_input.regulators) <= 6 "no more that 6 regulators are supported"

    metab_names = generate_metab_names(processed_input)
    param_names = generate_param_names(processed_input)

    enz = NamedTuple()
    for (i, substrate) in enumerate(processed_input[:substrates])
        enz = merge(enz, (; Symbol(:S, i) => substrate))
    end
    for (i, product) in enumerate(processed_input[:products])
        enz = merge(enz, (; Symbol(:P, i) => product))
    end
    for (i, regulator) in enumerate(processed_input[:regulators])
        enz = merge(enz, (; Symbol(:R, i) => regulator))
    end

    #TODO: use Base.method_argnames(methods(general_mwc_rate_equation)[1])[2:end] to get args
    mwc_rate_eq_args = [:S1, :S2, :S3, :P1, :P2, :P3, :R1, :R2, :R3, :R4, :R5, :R6]
    missing_keys = filter(x -> !haskey(enz, x), mwc_rate_eq_args)
    for key in missing_keys
        enz = merge(enz, (; key => nothing))
    end
    # qualified_name = esc(GlobalRef(Main, :rate_equation))
    function_name = (
        hasproperty(processed_input, :rate_equation_name) ?
        esc(processed_input.rate_equation_name) : esc(:rate_equation)
    )
    oligomeric_state = esc(processed_input.oligomeric_state)
    return quote
        @inline function $(function_name)(metabs, params, Keq)
            general_mwc_rate_equation(
                $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 1.0,
                $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 1.0,
                $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 1.0,
                $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 1.0,
                $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 1.0,
                $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 1.0,
                params.L,
                params.Vmax_a,
                params.Vmax_i,
                $(enz.S1 isa Symbol) ? params.$(Symbol("K_a_", enz.S1)) : 1.0,
                $(enz.S1 isa Symbol) ? params.$(Symbol("K_i_", enz.S1)) : 1.0,
                $(enz.S2 isa Symbol) ? params.$(Symbol("K_a_", enz.S2)) : 1.0,
                $(enz.S2 isa Symbol) ? params.$(Symbol("K_i_", enz.S2)) : 1.0,
                $(enz.S3 isa Symbol) ? params.$(Symbol("K_a_", enz.S3)) : 1.0,
                $(enz.S3 isa Symbol) ? params.$(Symbol("K_i_", enz.S3)) : 1.0,
                $(enz.P1 isa Symbol) ? params.$(Symbol("K_a_", enz.P1)) : 1.0,
                $(enz.P1 isa Symbol) ? params.$(Symbol("K_i_", enz.P1)) : 1.0,
                $(enz.P2 isa Symbol) ? params.$(Symbol("K_a_", enz.P2)) : 1.0,
                $(enz.P2 isa Symbol) ? params.$(Symbol("K_i_", enz.P2)) : 1.0,
                $(enz.P3 isa Symbol) ? params.$(Symbol("K_a_", enz.P3)) : 1.0,
                $(enz.P3 isa Symbol) ? params.$(Symbol("K_i_", enz.P3)) : 1.0,
                #Z_a_cat
                calculate_z_cat(
                    $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 0.0,
                    $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 0.0,
                    $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 0.0,
                    $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 0.0,
                    $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 0.0,
                    $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 0.0,
                    $(enz.S1 isa Symbol) ? params.$(Symbol("K_a_", enz.S1)) : Inf,
                    $(enz.S2 isa Symbol) ? params.$(Symbol("K_a_", enz.S2)) : Inf,
                    $(enz.S3 isa Symbol) ? params.$(Symbol("K_a_", enz.S3)) : Inf,
                    $(enz.P1 isa Symbol) ? params.$(Symbol("K_a_", enz.P1)) : Inf,
                    $(enz.P2 isa Symbol) ? params.$(Symbol("K_a_", enz.P2)) : Inf,
                    $(enz.P3 isa Symbol) ? params.$(Symbol("K_a_", enz.P3)) : Inf,
                    $(enz.S1 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P2)) : 0.0,
                    $(enz.S1 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P3)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P1)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P3)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P1)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P2)) : 0.0,
                ),
                #Z_i_cat
                calculate_z_cat(
                    $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 0.0,
                    $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 0.0,
                    $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 0.0,
                    $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 0.0,
                    $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 0.0,
                    $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 0.0,
                    $(enz.S1 isa Symbol) ? params.$(Symbol("K_i_", enz.S1)) : Inf,
                    $(enz.S2 isa Symbol) ? params.$(Symbol("K_i_", enz.S2)) : Inf,
                    $(enz.S3 isa Symbol) ? params.$(Symbol("K_i_", enz.S3)) : Inf,
                    $(enz.P1 isa Symbol) ? params.$(Symbol("K_i_", enz.P1)) : Inf,
                    $(enz.P2 isa Symbol) ? params.$(Symbol("K_i_", enz.P2)) : Inf,
                    $(enz.P3 isa Symbol) ? params.$(Symbol("K_i_", enz.P3)) : Inf,
                    $(enz.S1 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P2)) : 0.0,
                    $(enz.S1 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P3)) : 0.0,
                    $(enz.S1 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P1)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P3)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P1)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P2)) : 0.0,
                ),
                #Z_a_reg
                calculate_z_reg(
                    $(enz.R1 isa Symbol) ? metabs.$(enz.R1) : 0.0,
                    $(enz.R2 isa Symbol) ? metabs.$(enz.R2) : 0.0,
                    $(enz.R3 isa Symbol) ? metabs.$(enz.R3) : 0.0,
                    $(enz.R4 isa Symbol) ? metabs.$(enz.R4) : 0.0,
                    $(enz.R5 isa Symbol) ? metabs.$(enz.R5) : 0.0,
                    $(enz.R6 isa Symbol) ? metabs.$(enz.R6) : 0.0,
                    $(enz.R1 isa Symbol) ? params.$(Symbol("K_a_", enz.R1)) : Inf,
                    $(enz.R2 isa Symbol) ? params.$(Symbol("K_a_", enz.R2)) : Inf,
                    $(enz.R3 isa Symbol) ? params.$(Symbol("K_a_", enz.R3)) : Inf,
                    $(enz.R4 isa Symbol) ? params.$(Symbol("K_a_", enz.R4)) : Inf,
                    $(enz.R5 isa Symbol) ? params.$(Symbol("K_a_", enz.R5)) : Inf,
                    $(enz.R6 isa Symbol) ? params.$(Symbol("K_a_", enz.R6)) : Inf,
                ),
                #Z_i_reg
                calculate_z_reg(
                    $(enz.R1 isa Symbol) ? metabs.$(enz.R1) : 0.0,
                    $(enz.R2 isa Symbol) ? metabs.$(enz.R2) : 0.0,
                    $(enz.R3 isa Symbol) ? metabs.$(enz.R3) : 0.0,
                    $(enz.R4 isa Symbol) ? metabs.$(enz.R4) : 0.0,
                    $(enz.R5 isa Symbol) ? metabs.$(enz.R5) : 0.0,
                    $(enz.R6 isa Symbol) ? metabs.$(enz.R6) : 0.0,
                    $(enz.R1_reg3 isa Symbol) ? metabs.$(enz.R1_reg3) : 0.0,
                    $(enz.R2_reg3 isa Symbol) ? metabs.$(enz.R2_reg3) : 0.0,
                    $(enz.R3_reg3 isa Symbol) ? metabs.$(enz.R3_reg3) : 0.0,
                    $(enz.R1 isa Symbol) ? params.$(Symbol("K_i_", enz.R1)) : Inf,
                    $(enz.R2 isa Symbol) ? params.$(Symbol("K_i_", enz.R2)) : Inf,
                    $(enz.R3 isa Symbol) ? params.$(Symbol("K_i_", enz.R3)) : Inf,
                    $(enz.R4 isa Symbol) ? params.$(Symbol("K_i_", enz.R4)) : Inf,
                    $(enz.R5 isa Symbol) ? params.$(Symbol("K_i_", enz.R5)) : Inf,
                    $(enz.R6 isa Symbol) ? params.$(Symbol("K_i_", enz.R6)) : Inf,
                ),
                Keq,
                $(oligomeric_state),
            )
        end

        metab_names = $(esc(metab_names))
        param_names = $(esc(param_names))

        metab_names, param_names
    end
end

@inline function general_mwc_rate_equation(
    S1,
    S2,
    S3,
    P1,
    P2,
    P3,
    L,
    Vmax_a,
    Vmax_i,
    K_a_S1,
    K_i_S1,
    K_a_S2,
    K_i_S2,
    K_a_S3,
    K_i_S3,
    K_a_P1,
    K_i_P1,
    K_a_P2,
    K_i_P2,
    K_a_P3,
    K_i_P3,
    Z_a_cat,
    Z_i_cat,
    Z_a_reg,
    Z_i_reg,
    Keq,
    n,
)
    Vmax_a = 1.0
    Vmax_a_rev = ifelse(
        !isinf(K_a_P1 * K_a_P2 * K_a_P3),
        Vmax_a * K_a_P1 * K_a_P2 * K_a_P3 / (Keq * K_a_S1 * K_a_S2 * K_a_S3),
        0.0,
    )
    Vmax_i_rev = ifelse(
        !isinf(K_i_P1 * K_i_P2 * K_i_P3),
        Vmax_i * K_i_P1 * K_i_P2 * K_i_P3 / (Keq * K_i_S1 * K_i_S2 * K_i_S3),
        0.0,
    )
    Rate =
        (
            (
                Vmax_a * (S1 / K_a_S1) * (S2 / K_a_S2) * (S3 / K_a_S3) -
                Vmax_a_rev * (P1 / K_a_P1) * (P2 / K_a_P2) * (P3 / K_a_P3)
            ) *
            (Z_a_cat^(n - 1)) *
            (Z_a_reg^n) +
            L *
            (
                Vmax_i * (S1 / K_i_S1) * (S2 / K_i_S2) * (S3 / K_i_S3) -
                Vmax_i_rev * (P1 / K_i_P1) * (P2 / K_i_P2) * (P3 / K_i_P3)
            ) *
            (Z_i_cat^(n - 1)) *
            (Z_i_reg^n)
        ) / ((Z_a_cat^n) * (Z_a_reg^n) + L * (Z_i_cat^n) * (Z_i_reg^n))

    return Rate
end

@inline function calculate_z_reg(
    R1,
    R2,
    R3,
    R4,
    R5,
    R6,
    K_R1,
    K_R2,
    K_R3,
    K_R4,
    K_R5,
    K_R6,
    delta_R1_R2,
    delta_R1_R3,
    delta_R1_R4,
    delta_R1_R5,
    delta_R1_R6,
    delta_R2_R3,
    delta_R2_R4,
    delta_R2_R5,
    delta_R2_R6,
    delta_R3_R4,
    delta_R3_R5,
    delta_R3_R6,
    delta_R4_R5,
    delta_R4_R6,
    delta_R5_R6,
)
    # Z_reg = (
    #     (1 + R1 / K_R1) *
    #     (1 + R2 / K_R2) *
    #     (1 + R3 / K_R3) *
    #     (1 + R4 / K_R4) *
    #     (1 + R5 / K_R5) *
    #     (1 + R6 / K_R6)
    # )
    Z_reg = (
        (
            1 +
            R1 / K_R1 +
            delta_R1_R2 * R2 / K_R2 +
            delta_R1_R3 * R3 / K_R3 +
            delta_R1_R4 * R4 / K_R4 +
            delta_R1_R5 * R5 / K_R5 +
            delta_R1_R6 * R6 / K_R6
        ) *
        (
            1 +
            (1 - delta_R1_R2) * (
                R2 / K_R2 +
                delta_R2_R3 * R3 / K_R3 +
                delta_R2_R4 * R4 / K_R4 +
                delta_R2_R5 * R5 / K_R5 +
                delta_R2_R6 * R6 / K_R6
            )
        ) *
        (
            1 +
            (1 - delta_R1_R3) *
            (1 - delta_R2_R3) *
            (
                R3 / K_R3 +
                delta_R3_R4 * R4 / K_R4 +
                delta_R3_R5 * R5 / K_R5 +
                delta_R3_R6 * R6 / K_R6
            )
        ) *
        (
            1 +
            (1 - delta_R1_R4) *
            (1 - delta_R2_R4) *
            (1 - delta_R3_R4) *
            (R4 / K_R4 + delta_R4_R5 * R5 / K_R5 + delta_R4_R6 * R6 / K_R6)
        ) *
        (
            1 +
            (1 - delta_R1_R5) *
            (1 - delta_R2_R5) *
            (1 - delta_R3_R5) *
            (1 - delta_R4_R5) *
            (R5 / K_R5 + delta_R5_R6 * R6 / K_R6)
        ) *
        (
            1 +
            (1 - delta_R1_R6) *
            (1 - delta_R2_R6) *
            (1 - delta_R3_R6) *
            (1 - delta_R4_R6) *
            (1 - delta_R5_R6) *
            (R6 / K_R6)
        )
    )
    return Z_reg
end

@inline function calculate_z_cat(
    S1,
    S2,
    S3,
    P1,
    P2,
    P3,
    K_S1,
    K_S2,
    K_S3,
    K_P1,
    K_P2,
    K_P3,
    alpha_S1_P1,
    alpha_S1_P2,
    alpha_S1_P3,
    alpha_S2_P1,
    alpha_S2_P2,
    alpha_S2_P3,
    alpha_S3_P1,
    alpha_S3_P2,
    alpha_S3_P3,
)
    Z_cat = (
        1 +
        (((1 + S1 / K_S1) * (1 + S2 / K_S2) * (1 + S3 / K_S3)) - 1) +
        (((1 + P1 / K_P1) * (1 + P2 / K_P2) * (1 + P3 / K_P3)) - 1) +
        (S1 / K_S1) * (
            alpha_S1_P1 * (P1 / K_P1) +
            alpha_S1_P2 * (P2 / K_P2) +
            alpha_S1_P3 * (P3 / K_P3)
        ) +
        (S2 / K_S2) * (
            alpha_S2_P1 * (P1 / K_P1) +
            alpha_S2_P2 * (P2 / K_P2) +
            alpha_S2_P3 * (P3 / K_P3)
        ) +
        (S3 / K_S3) * (
            alpha_S3_P1 * (P1 / K_P1) +
            alpha_S3_P2 * (P2 / K_P2) +
            alpha_S3_P3 * (P3 / K_P3)
        ) +
        (S1 / K_S1) *
        (S2 / K_S2) *
        (
            alpha_S1_P1 * alpha_S2_P1 * (P1 / K_P1) +
            alpha_S1_P2 * alpha_S2_P2 * (P2 / K_P2) +
            alpha_S1_P3 * alpha_S2_P3 * (P3 / K_P3)
        ) +
        (S1 / K_S1) *
        (S3 / K_S3) *
        (
            alpha_S1_P1 * alpha_S3_P1 * (P1 / K_P1) +
            alpha_S1_P2 * alpha_S3_P2 * (P2 / K_P2) +
            alpha_S1_P3 * alpha_S3_P3 * (P3 / K_P3)
        ) +
        (S2 / K_S2) *
        (S3 / K_S3) *
        (
            alpha_S2_P1 * alpha_S3_P1 * (P1 / K_P1) +
            alpha_S2_P2 * alpha_S3_P2 * (P2 / K_P2) +
            alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3)
        ) +
        (P1 / K_P1) *
        (P2 / K_P2) *
        (
            alpha_S1_P1 * alpha_S2_P1 * (S1 / K_S1) +
            alpha_S1_P2 * alpha_S2_P2 * (S2 / K_S2) +
            alpha_S1_P3 * alpha_S2_P3 * (S3 / K_S3)
        ) +
        (P1 / K_P1) *
        (P3 / K_P3) *
        (
            alpha_S1_P1 * alpha_S3_P1 * (S1 / K_S1) +
            alpha_S1_P2 * alpha_S3_P2 * (S2 / K_S2) +
            alpha_S1_P3 * alpha_S3_P3 * (S3 / K_S3)
        ) +
        (P2 / K_P2) *
        (P3 / K_P3) *
        (
            alpha_S2_P1 * alpha_S3_P1 * (S1 / K_S1) +
            alpha_S2_P2 * alpha_S3_P2 * (S2 / K_S2) +
            alpha_S2_P3 * alpha_S3_P3 * (S3 / K_S3)
        )
    )
    return Z_cat
end

"Generate the names of the parameters for the rate equation using the same input as @derive_general_mwc_rate_eq"
function generate_param_names(processed_input)
    param_names = (:L, :Vmax_a, :Vmax_i)
    for metab in [
        processed_input[:substrates]...,
        processed_input[:products]...,
        processed_input[:regulators]...,
    ]
        param_names = (
            param_names...,
            Symbol("K_a_", metab, "_", site),
            Symbol("K_i_", metab, "_", site),
        )
    end
    for substrate in processed_input[:substrates]
        for product in processed_input[:products]
            param_names = (param_names..., Symbol("alpha_", substrate, "_", product))
        end
    end
    for (i, regulator1) in enumerate(processed_input[:regulators])[1:(end-1)]
        for regulator2 in processed_input[:regulators][(i+1):end]
            param_names = (param_names..., Symbol("delta_", regulator1, "_", regulator2))
        end
    end
    return param_names
end

"Generate the names of the metabolites for the rate equation using the same input as @derive_general_mwc_rate_eq"
function generate_metab_names(processed_input)
    metab_names = ()
    for field in keys(processed_input)
        if field ∈ [:substrates, :products, :regulators]
            metab_names = (metab_names..., processed_input[field]...)
        end
    end
    return Tuple(unique(metab_names))
end
