#=
CODE FOR RATE EQUATION DERIVATION
=#
using Distributed

@doc raw"""
    derive_general_rapid_equilibrium_rate_eq(metabs_and_regulators_kwargs...)

Derive a function that calculates the rate of a reaction using the general rate equation using rapid equilibrium approximation given the list of substrates and products that bind to specific catalytic sites.

The general rapid equilibrium rate equation is given by:

```math
Rate = \\frac{{V_{max} \\prod_{i=1}^{n} \\left(\\frac{S_i}{K_{S_i}}\\right) - V_{max}_{rev} \\prod_{i=1}^{n} \\left(\\frac{P_i}{K_{P_i}}\\right) - V_{max rev} \\prod_{i=1}^{n} \\left(\\frac{P_i}{K_{P_i}}\\right)\\right)}{Z_{cat}}
```

where:
- ``V_{max}`` is the maximum rate of the forward reaction
- ``V_{max rev}`` is the maximum rate of the reverse reaction
- ``S_i`` is the concentration of the ``i^{th}`` substrate
- ``P_i`` is the concentration of the ``i^{th}`` product
- ``K_{S_i}`` is the binding constant of the ``i^{th}`` substrate to enzyme
- ``K_{S_i}`` is the binding constant of the ``i^{th}`` product to enzyme
- ``Z_{cat}`` is the allosteric factor for the catalytic site

# Arguments
- `metabs_and_regulators_kwargs...`: keyword arguments that specify the substrates, products, catalytic sites, regulatory sites, and other parameters of the reaction.

# Returns
- A function that calculates the rate of the reaction using the general rapid_equilibrium rate equation
- A tuple of the names of the metabolites and parameters used in the rate equation
"""
macro derive_general_rapid_equilibrium_rate_eq(metabs_and_regulators_kwargs)

    processed_input = getfield(__module__, Symbol(metabs_and_regulators_kwargs))
    expected_input_kwargs = [:substrates, :products, :Keq, :rate_equation_name]
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

    metab_names = rapid_equilibrium_metab_names(processed_input)
    param_names = rapid_equilibrium_param_names(processed_input)

    enz = NamedTuple()
    for (i, substrate) in enumerate(processed_input[:substrates])
        enz = merge(enz, (; Symbol(:S, i) => substrate))
    end
    for (i, product) in enumerate(processed_input[:products])
        enz = merge(enz, (; Symbol(:P, i) => product))
    end

    rapid_equilibrium_rate_eq_args = [:S1, :S2, :S3, :P1, :P2, :P3]
    missing_keys = filter(x -> !haskey(enz, x), rapid_equilibrium_rate_eq_args)
    for key in missing_keys
        enz = merge(enz, (; key => nothing))
    end
    # qualified_name = esc(GlobalRef(Main, :rate_equation))
    function_name = (
        hasproperty(processed_input, :rate_equation_name) ?
        esc(processed_input.rate_equation_name) : esc(:rate_equation)
    )
    return quote
        @inline function $(function_name)(metabs, params, Keq)
            general_rapid_equilibrium_rate_equation(
                $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 1.0,
                $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 1.0,
                $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 1.0,
                $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 1.0,
                $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 1.0,
                $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 1.0,
                params.Vmax,
                $(enz.S1 isa Symbol) ? params.$(Symbol("K_", enz.S1)) : 1.0,
                $(enz.S2 isa Symbol) ? params.$(Symbol("K_", enz.S2)) : 1.0,
                $(enz.S3 isa Symbol) ? params.$(Symbol("K_", enz.S3)) : 1.0,
                $(enz.P1 isa Symbol) ? params.$(Symbol("K_", enz.P1)) : 1.0,
                $(enz.P2 isa Symbol) ? params.$(Symbol("K_", enz.P2)) : 1.0,
                $(enz.P3 isa Symbol) ? params.$(Symbol("K_", enz.P3)) : 1.0,
                #Z_cat
                calculate_z_cat_rapid_equilibrium(
                    $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 0.0,
                    $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 0.0,
                    $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 0.0,
                    $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 0.0,
                    $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 0.0,
                    $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 0.0,
                    $(enz.S1 isa Symbol) ? params.$(Symbol("K_", enz.S1)) : Inf,
                    $(enz.S2 isa Symbol) ? params.$(Symbol("K_", enz.S2)) : Inf,
                    $(enz.S3 isa Symbol) ? params.$(Symbol("K_", enz.S3)) : Inf,
                    $(enz.P1 isa Symbol) ? params.$(Symbol("K_", enz.P1)) : Inf,
                    $(enz.P2 isa Symbol) ? params.$(Symbol("K_", enz.P2)) : Inf,
                    $(enz.P3 isa Symbol) ? params.$(Symbol("K_", enz.P3)) : Inf,
                    $(enz.S1 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P1)) : 0.0,
                    $(enz.S1 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P2)) : 0.0,
                    $(enz.S1 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P3)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P1)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P2)) : 0.0,
                    $(enz.S2 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P3)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P1)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P2)) : 0.0,
                    $(enz.S3 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P3)) : 0.0,
                ),
                Keq,
            )
        end

        metab_names = $(esc(metab_names))
        param_names = $(esc(param_names))

        metab_names, param_names
    end
end

@inline function general_rapid_equilibrium_rate_equation(
    S1,
    S2,
    S3,
    P1,
    P2,
    P3,
    Vmax,
    K_S1,
    K_S2,
    K_S3,
    K_P1,
    K_P2,
    K_P3,
    Z_cat,
    Keq,
)
    Vmax = 1.0
    Vmax_rev = ifelse(
        !isinf(K_P1 * K_P2 * K_P3),
        Vmax * K_P1 * K_P2 * K_P3 / (Keq * K_S1 * K_S2 * K_S3),
        0.0,
    )
    Rate =
        (
            Vmax * (S1 / K_S1) * (S2 / K_S2) * (S3 / K_S3) -
            Vmax_rev * (P1 / K_P1) * (P2 / K_P2) * (P3 / K_P3)
        ) / Z_cat

    return Rate
end

@inline function calculate_z_cat_rapid_equilibrium(
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

"Generate the names of the parameters for the rate equation using the same input as @derive_general_rapid_equilibrium_rate_eq"
function rapid_equilibrium_param_names(processed_input)
    param_names = (:Vmax,)
    for metab in [processed_input[:substrates]..., processed_input[:products]...]
        param_names = (param_names..., Symbol("K_", metab))
    end
    for substrate in processed_input[:substrates]
        for product in processed_input[:products]
            param_names = (param_names..., Symbol("alpha_", substrate, "_", product))
        end
    end

    return param_names
end

"Generate the names of the metabolites for the rate equation using the same input as @derive_general_rapid_equilibrium_rate_eq"
function rapid_equilibrium_metab_names(processed_input)
    metab_names = ()
    for site in keys(processed_input)
        if site ∈ [:substrates, :products]
            metab_names = (metab_names..., processed_input[site]...)
        end
    end
    return Tuple(unique(metab_names))
end
