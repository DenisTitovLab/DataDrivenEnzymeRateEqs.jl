#=
CODE FOR MWC RATE EQUATION DERIVATION
=#
using Distributed

@doc raw"""
    derive_general_mwc_rate_eq(metabs_and_regulators_kwargs...)

Derive a function that calculates the rate of a reaction using the general MWC rate equation given the list of substrates, products, and regulators that bind to specific cat or reg sites.

The general MWC rate equation is given by:

```math
Rate = \frac{{V_{max}^a \prod_{i=1}^{n} \left(\frac{S_i}{K_{a, i}}\right) - V_{max, rev}^a \prod_{i=1}^{n} \left(\frac{P_i}{K_{a, i}}\right) \cdot Z_{a, cat}^{n-1} \cdot Z_{a, reg}^n + L \left(V_{max}^i \prod_{i=1}^{n} \left(\frac{S_i}{K_{i, i}}\right) - V_{max, rev}^i \prod_{i=1}^{n} \left(\frac{P_i}{K_{i, i}}\right)\right) \cdot Z_{i, cat}^{n-1} \cdot Z_{i, reg}^n}}{Z_{a, cat}^n \cdot Z_{a, reg}^n + L \cdot Z_{i, cat}^n \cdot Z_{i, reg}^n}
```

where:
- ``V_{max}^a`` is the maximum rate of the forward reaction
- ``V_{max, rev}^a`` is the maximum rate of the reverse reaction
- ``V_{max}^i`` is the maximum rate of the forward reaction
- ``V_{max, rev}^i`` is the maximum rate of the reverse reaction
- ``S_i`` is the concentration of the ``i^{th}`` substrate
- ``P_i`` is the concentration of the ``i^{th}`` product
- ``I_i`` is the concentration of the ``i^{th}`` catalytic site inhibitor
- ``R_i`` is the concentration of the ``i^{th}`` allosteric regulator
- ``K_{a, X}`` is the binding constant of the ``X`` metabolite for active MWC state
- ``K_{i, X}`` is the binding constant of the ``X`` metabolite for inactive MWC state
- ``Z_{a, cat}`` is the allosteric factor for the catalytic site in the active MWC state
- ``Z_{i, cat}`` is the allosteric factor for the catalytic site in the inactive MWC state
- ``Z_{a, reg}`` is the allosteric factor for the regulatory site in the active MWC state
- ``Z_{i, reg}`` is the allosteric factor for the regulatory site in the inactive MWC state
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
    expected_input_kwargs = [
        :substrates,
        :products,
        :inhibitors,
        :regulators,
        :Keq,
        :oligomeric_state,
        :rate_equation_name,
    ]
    if !haskey(processed_input, :regulators)
        processed_input.regulators = []
    elseif !haskey(processed_input, :inhibitors)
        processed_input.inhibitors = []
    end
    for kwarg in keys(processed_input)
        kwarg ∈ expected_input_kwargs || error(
            "invalid keyword: ",
            kwarg,
            ". The only supported keywords are: ",
            expected_input_kwargs,
        )
        @assert haskey(processed_input, kwarg) "missing keyword: $kwarg"
    end

    @assert 0 < length(processed_input.substrates) <= 3 "At least 1 and no more that 3 substrates are supported"
    @assert 0 < length(processed_input.products) <= 3 "At least 1 and no more that 3 products are supported"
    @assert 0 <= length(processed_input.inhibitors) <= 6 "no more that 6 inhibitors are supported"
    @assert 0 <= length(processed_input.regulators) <= 6 "no more that 6 regulators are supported"

    metab_names = generate_metab_names(processed_input)
    param_names = generate_param_names(processed_input)

    enz = NamedTuple()
    for (i, substrate) in enumerate(processed_input[:substrates])
        enz = merge(enz, (; Symbol(:S, i) => substrate))
    end
    for (i, product) in enumerate(processed_input[:products])
        enz = merge(enz, (; Symbol(:P, i) => product))
    end
    for (i, inhibitor) in enumerate(processed_input[:inhibitors])
        enz = merge(enz, (; Symbol(:I, i) => inhibitor))
    end
    for (i, regulator) in enumerate(processed_input[:regulators])
        enz = merge(enz, (; Symbol(:R, i) => regulator))
    end

    mwc_rate_eq_args = [
        :S1,
        :S2,
        :S3,
        :P1,
        :P2,
        :P3,
        :I1,
        :I2,
        :I3,
        :I4,
        :I5,
        :I6,
        :R1,
        :R2,
        :R3,
        :R4,
        :R5,
        :R6,
    ]
    missing_keys = filter(x -> !haskey(enz, x), mwc_rate_eq_args)
    for key in missing_keys
        enz = merge(enz, (; key => nothing))
    end
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
                    $(enz.I1 isa Symbol) ? metabs.$(enz.I1) : 0.0,
                    $(enz.I2 isa Symbol) ? metabs.$(enz.I2) : 0.0,
                    $(enz.I3 isa Symbol) ? metabs.$(enz.I3) : 0.0,
                    $(enz.I4 isa Symbol) ? metabs.$(enz.I4) : 0.0,
                    $(enz.I5 isa Symbol) ? metabs.$(enz.I5) : 0.0,
                    $(enz.I6 isa Symbol) ? metabs.$(enz.I6) : 0.0,
                    $(enz.S1 isa Symbol) ? params.$(Symbol("K_a_", enz.S1)) : Inf,
                    $(enz.S2 isa Symbol) ? params.$(Symbol("K_a_", enz.S2)) : Inf,
                    $(enz.S3 isa Symbol) ? params.$(Symbol("K_a_", enz.S3)) : Inf,
                    $(enz.P1 isa Symbol) ? params.$(Symbol("K_a_", enz.P1)) : Inf,
                    $(enz.P2 isa Symbol) ? params.$(Symbol("K_a_", enz.P2)) : Inf,
                    $(enz.P3 isa Symbol) ? params.$(Symbol("K_a_", enz.P3)) : Inf,
                    $(enz.I1 isa Symbol) ? params.$(Symbol("K_a_", enz.I1, "_cat")) : Inf,
                    $(enz.I2 isa Symbol) ? params.$(Symbol("K_a_", enz.I2, "_cat")) : Inf,
                    $(enz.I3 isa Symbol) ? params.$(Symbol("K_a_", enz.I3, "_cat")) : Inf,
                    $(enz.I4 isa Symbol) ? params.$(Symbol("K_a_", enz.I4, "_cat")) : Inf,
                    $(enz.I5 isa Symbol) ? params.$(Symbol("K_a_", enz.I5, "_cat")) : Inf,
                    $(enz.I6 isa Symbol) ? params.$(Symbol("K_a_", enz.I6, "_cat")) : Inf,
                    $(enz.S1 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I4)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I5)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I6)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I4)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I5)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I6)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I4)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I5)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I6)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I1)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I2)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I3)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I4)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I5)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I6)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I1)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I2)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I3)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I4)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I5)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I6)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I1)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I2)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I3)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I4)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I5)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I6)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I2)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I3)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I4)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I5)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I6)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I3)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I4)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I5)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I6)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I4)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I5)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I6)) : 1.0,
                    $(enz.I4 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I4, "_", enz.I5)) : 1.0,
                    $(enz.I4 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I4, "_", enz.I6)) : 1.0,
                    $(enz.I5 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I5, "_", enz.I6)) : 1.0,
                ),
                #Z_i_cat
                calculate_z_cat(
                    $(enz.S1 isa Symbol) ? metabs.$(enz.S1) : 0.0,
                    $(enz.S2 isa Symbol) ? metabs.$(enz.S2) : 0.0,
                    $(enz.S3 isa Symbol) ? metabs.$(enz.S3) : 0.0,
                    $(enz.P1 isa Symbol) ? metabs.$(enz.P1) : 0.0,
                    $(enz.P2 isa Symbol) ? metabs.$(enz.P2) : 0.0,
                    $(enz.P3 isa Symbol) ? metabs.$(enz.P3) : 0.0,
                    $(enz.I1 isa Symbol) ? metabs.$(enz.I1) : 0.0,
                    $(enz.I2 isa Symbol) ? metabs.$(enz.I2) : 0.0,
                    $(enz.I3 isa Symbol) ? metabs.$(enz.I3) : 0.0,
                    $(enz.I4 isa Symbol) ? metabs.$(enz.I4) : 0.0,
                    $(enz.I5 isa Symbol) ? metabs.$(enz.I5) : 0.0,
                    $(enz.I6 isa Symbol) ? metabs.$(enz.I6) : 0.0,
                    $(enz.S1 isa Symbol) ? params.$(Symbol("K_i_", enz.S1)) : Inf,
                    $(enz.S2 isa Symbol) ? params.$(Symbol("K_i_", enz.S2)) : Inf,
                    $(enz.S3 isa Symbol) ? params.$(Symbol("K_i_", enz.S3)) : Inf,
                    $(enz.P1 isa Symbol) ? params.$(Symbol("K_i_", enz.P1)) : Inf,
                    $(enz.P2 isa Symbol) ? params.$(Symbol("K_i_", enz.P2)) : Inf,
                    $(enz.P3 isa Symbol) ? params.$(Symbol("K_i_", enz.P3)) : Inf,
                    $(enz.I1 isa Symbol) ? params.$(Symbol("K_i_", enz.I1, "_cat")) : Inf,
                    $(enz.I2 isa Symbol) ? params.$(Symbol("K_i_", enz.I2, "_cat")) : Inf,
                    $(enz.I3 isa Symbol) ? params.$(Symbol("K_i_", enz.I3, "_cat")) : Inf,
                    $(enz.I4 isa Symbol) ? params.$(Symbol("K_i_", enz.I4, "_cat")) : Inf,
                    $(enz.I5 isa Symbol) ? params.$(Symbol("K_i_", enz.I5, "_cat")) : Inf,
                    $(enz.I6 isa Symbol) ? params.$(Symbol("K_i_", enz.I6, "_cat")) : Inf,
                    $(enz.S1 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.P3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I4)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I5)) : 1.0,
                    $(enz.S1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.I6)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.P3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I4)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I5)) : 1.0,
                    $(enz.S2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.I6)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.P3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.P3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I4)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I5)) : 1.0,
                    $(enz.S3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.I6)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I1)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I2)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I3)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I4)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I5)) : 1.0,
                    $(enz.P1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.I6)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I1)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I2)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I3)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I4)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I5)) : 1.0,
                    $(enz.P2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.I6)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I1)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I2)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I3)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I4)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I5)) : 1.0,
                    $(enz.P3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.I6)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I2)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I3)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I4)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I5)) : 1.0,
                    $(enz.I1 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I1, "_", enz.I6)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I3)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I4)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I5)) : 1.0,
                    $(enz.I2 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I2, "_", enz.I6)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I4)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I5)) : 1.0,
                    $(enz.I3 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I3, "_", enz.I6)) : 1.0,
                    $(enz.I4 isa Symbol && enz.I5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I4, "_", enz.I5)) : 1.0,
                    $(enz.I4 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I4, "_", enz.I6)) : 1.0,
                    $(enz.I5 isa Symbol && enz.I6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.I5, "_", enz.I6)) : 1.0,
                ),
                #Z_a_reg
                calculate_z_reg(
                    $(enz.R1 isa Symbol) ? metabs.$(enz.R1) : 0.0,
                    $(enz.R2 isa Symbol) ? metabs.$(enz.R2) : 0.0,
                    $(enz.R3 isa Symbol) ? metabs.$(enz.R3) : 0.0,
                    $(enz.R4 isa Symbol) ? metabs.$(enz.R4) : 0.0,
                    $(enz.R5 isa Symbol) ? metabs.$(enz.R5) : 0.0,
                    $(enz.R6 isa Symbol) ? metabs.$(enz.R6) : 0.0,
                    $(enz.R1 isa Symbol) ? params.$(Symbol("K_a_", enz.R1, "_reg")) : Inf,
                    $(enz.R2 isa Symbol) ? params.$(Symbol("K_a_", enz.R2, "_reg")) : Inf,
                    $(enz.R3 isa Symbol) ? params.$(Symbol("K_a_", enz.R3, "_reg")) : Inf,
                    $(enz.R4 isa Symbol) ? params.$(Symbol("K_a_", enz.R4, "_reg")) : Inf,
                    $(enz.R5 isa Symbol) ? params.$(Symbol("K_a_", enz.R5, "_reg")) : Inf,
                    $(enz.R6 isa Symbol) ? params.$(Symbol("K_a_", enz.R6, "_reg")) : Inf,
                    $(enz.S1 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R4)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R5)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R6)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R4)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R5)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R6)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R4)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R5)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R6)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R1)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R2)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R3)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R4)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R5)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R6)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R1)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R2)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R3)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R4)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R5)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R6)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R1)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R2)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R3)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R4)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R5)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R6)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R2)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R3)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R4)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R5)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R6)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R3)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R4)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R5)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R6)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R4)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R5)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R6)) : 1.0,
                    $(enz.R4 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R4, "_", enz.R5)) : 1.0,
                    $(enz.R4 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R4, "_", enz.R6)) : 1.0,
                    $(enz.R5 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R5, "_", enz.R6)) : 1.0,
                ),
                #Z_i_reg
                calculate_z_reg(
                    $(enz.R1 isa Symbol) ? metabs.$(enz.R1) : 0.0,
                    $(enz.R2 isa Symbol) ? metabs.$(enz.R2) : 0.0,
                    $(enz.R3 isa Symbol) ? metabs.$(enz.R3) : 0.0,
                    $(enz.R4 isa Symbol) ? metabs.$(enz.R4) : 0.0,
                    $(enz.R5 isa Symbol) ? metabs.$(enz.R5) : 0.0,
                    $(enz.R6 isa Symbol) ? metabs.$(enz.R6) : 0.0,
                    $(enz.R1 isa Symbol) ? params.$(Symbol("K_i_", enz.R1, "_reg")) : Inf,
                    $(enz.R2 isa Symbol) ? params.$(Symbol("K_i_", enz.R2, "_reg")) : Inf,
                    $(enz.R3 isa Symbol) ? params.$(Symbol("K_i_", enz.R3, "_reg")) : Inf,
                    $(enz.R4 isa Symbol) ? params.$(Symbol("K_i_", enz.R4, "_reg")) : Inf,
                    $(enz.R5 isa Symbol) ? params.$(Symbol("K_i_", enz.R5, "_reg")) : Inf,
                    $(enz.R6 isa Symbol) ? params.$(Symbol("K_i_", enz.R6, "_reg")) : Inf,
                    $(enz.S1 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R1)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R2)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R3)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R4)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R5)) : 1.0,
                    $(enz.S1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1, "_", enz.R6)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R1)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R2)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R3)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R4)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R5)) : 1.0,
                    $(enz.S2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2, "_", enz.R6)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R1)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R2)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R3)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R4)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R5)) : 1.0,
                    $(enz.S3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3, "_", enz.R6)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R1)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R2)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R3)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R4)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R5)) : 1.0,
                    $(enz.P1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P1, "_", enz.R6)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R1)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R2)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R3)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R4)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R5)) : 1.0,
                    $(enz.P2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P2, "_", enz.R6)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R1)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R2)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R3)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R4)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R5)) : 1.0,
                    $(enz.P3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.P3, "_", enz.R6)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R2)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R3)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R4)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R5)) : 1.0,
                    $(enz.R1 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R1, "_", enz.R6)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R3)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R4)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R5)) : 1.0,
                    $(enz.R2 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R2, "_", enz.R6)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R4)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R5)) : 1.0,
                    $(enz.R3 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R3, "_", enz.R6)) : 1.0,
                    $(enz.R4 isa Symbol && enz.R5 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R4, "_", enz.R5)) : 1.0,
                    $(enz.R4 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R4, "_", enz.R6)) : 1.0,
                    $(enz.R5 isa Symbol && enz.R6 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.R5, "_", enz.R6)) : 1.0,
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

@inline function calculate_z_cat(
    S1,
    S2,
    S3,
    P1,
    P2,
    P3,
    I1,
    I2,
    I3,
    I4,
    K_S1,
    K_S2,
    K_S3,
    K_P1,
    K_P2,
    K_P3,
    K_I1,
    K_I2,
    K_I3,
    K_I4,
    alpha_S1_P1,
    alpha_S1_P2,
    alpha_S1_P3,
    alpha_S1_I1,
    alpha_S1_I2,
    alpha_S1_I3,
    alpha_S1_I4,
    alpha_S2_P1,
    alpha_S2_P2,
    alpha_S2_P3,
    alpha_S2_I1,
    alpha_S2_I2,
    alpha_S2_I3,
    alpha_S2_I4,
    alpha_S3_P1,
    alpha_S3_P2,
    alpha_S3_P3,
    alpha_S3_I1,
    alpha_S3_I2,
    alpha_S3_I3,
    alpha_S3_I4,
    alpha_P1_I1,
    alpha_P1_I2,
    alpha_P1_I3,
    alpha_P1_I4,
    alpha_P2_I1,
    alpha_P2_I2,
    alpha_P2_I3,
    alpha_P2_I4,
    alpha_P3_I1,
    alpha_P3_I2,
    alpha_P3_I3,
    alpha_P3_I4,
    alpha_I1_I2,
    alpha_I1_I3,
    alpha_I1_I4,
    alpha_I2_I3,
    alpha_I2_I4,
    alpha_I3_I4,
)
    one_metab_bound = (
        S1 / K_S1 +
        S2 / K_S2 +
        S3 / K_S3 +
        P1 / K_P1 +
        P2 / K_P2 +
        P3 / K_P3 +
        I1 / K_I1 +
        I2 / K_I2 +
        I3 / K_I3 +
        I4 / K_I4
    )
    two_metab_bound =
        (S1 / K_S1) * (S2 / K_S2) +
        (S1 / K_S1) * (S3 / K_S3) +
        (S2 / K_S2) * (S3 / K_S3) +
        (P1 / K_P1) * (P2 / K_P2) +
        (P1 / K_P1) * (P3 / K_P3) +
        (P2 / K_P2) * (P3 / K_P3) +
        (S1 / K_S1) * (
            alpha_S1_P1 * (P1 / K_P1) +
            alpha_S1_P2 * (P2 / K_P2) +
            alpha_S1_P3 * (P3 / K_P3) +
            alpha_S1_I1 * (I1 / K_I1) +
            alpha_S1_I2 * (I2 / K_I2) +
            alpha_S1_I3 * (I3 / K_I3) +
            alpha_S1_I4 * (I4 / K_I4)
        ) +
        (S2 / K_S2) * (
            alpha_S2_P1 * (P1 / K_P1) +
            alpha_S2_P2 * (P2 / K_P2) +
            alpha_S2_P3 * (P3 / K_P3) +
            alpha_S2_I1 * (I1 / K_I1) +
            alpha_S2_I2 * (I2 / K_I2) +
            alpha_S2_I3 * (I3 / K_I3) +
            alpha_S2_I4 * (I4 / K_I4)
        ) +
        (S3 / K_S3) * (
            alpha_S3_P1 * (P1 / K_P1) +
            alpha_S3_P2 * (P2 / K_P2) +
            alpha_S3_P3 * (P3 / K_P3) +
            alpha_S3_I1 * (I1 / K_I1) +
            alpha_S3_I2 * (I2 / K_I2) +
            alpha_S3_I3 * (I3 / K_I3) +
            alpha_S3_I4 * (I4 / K_I4)
        ) +
        (P1 / K_P1) * (
            alpha_P1_I1 * (I1 / K_I1) +
            alpha_P1_I2 * (I2 / K_I2) +
            alpha_P1_I3 * (I3 / K_I3) +
            alpha_P1_I4 * (I4 / K_I4)
        ) +
        (P2 / K_P2) * (
            alpha_P2_I1 * (I1 / K_I1) +
            alpha_P2_I2 * (I2 / K_I2) +
            alpha_P2_I3 * (I3 / K_I3) +
            alpha_P2_I4 * (I4 / K_I4)
        ) +
        (P3 / K_P3) * (
            alpha_P3_I1 * (I1 / K_I1) +
            alpha_P3_I2 * (I2 / K_I2) +
            alpha_P3_I3 * (I3 / K_I3) +
            alpha_P3_I4 * (I4 / K_I4)
        ) +
        (I1 / K_I1) * (
            alpha_I1_I2 * (I2 / K_I2) +
            alpha_I1_I3 * (I3 / K_I3) +
            alpha_I1_I4 * (I4 / K_I4)
        ) +
        (I2 / K_I2) * (alpha_I2_I3 * (I3 / K_I3) + alpha_I2_I4 * (I4 / K_I4)) +
        (I3 / K_I3) * (alpha_I3_I4 * (I4 / K_I4))

    three_metab_bound =
        (S1 / K_S1) * (
            (S2 / K_S2) * (
                (S3 / K_S3) +
                alpha_S1_P1 * alpha_S2_P1 * (P1 / K_P1) +
                alpha_S1_P2 * alpha_S2_P2 * (P2 / K_P2) +
                alpha_S1_P3 * alpha_S2_P3 * (P3 / K_P3) +
                alpha_S1_I1 * alpha_S2_I1 * (I1 / K_I1) +
                alpha_S1_I2 * alpha_S2_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_S2_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_S2_I4 * (I4 / K_I4)
            ) +
            (S3 / K_S3) * (
                alpha_S1_P1 * alpha_S3_P1 * (P1 / K_P1) +
                alpha_S1_P2 * alpha_S3_P2 * (P2 / K_P2) +
                alpha_S1_P3 * alpha_S3_P3 * (P3 / K_P3) +
                alpha_S1_I1 * alpha_S3_I1 * (I1 / K_I1) +
                alpha_S1_I2 * alpha_S3_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_S3_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_S3_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_P1 * (P1 / K_P1)) * (
                alpha_S1_P2 * (P2 / K_P2) +
                alpha_S1_P3 * (P3 / K_P3) +
                alpha_S1_I1 * alpha_P1_I1 * (I1 / K_I1) +
                alpha_S1_I2 * alpha_P1_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_P1_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_P1_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_P2 * (P2 / K_P2)) * (
                alpha_S1_P3 * (P3 / K_P3) +
                alpha_S1_I1 * alpha_P2_I1 * (I1 / K_I1) +
                alpha_S1_I2 * alpha_P2_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_P2_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_P2_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_P3 * (P3 / K_P3)) * (
                alpha_S1_I1 * alpha_P3_I1 * (I1 / K_I1) +
                alpha_S1_I2 * alpha_P3_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_P3_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_P3_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_I1 * (I1 / K_I1)) * (
                alpha_S1_I2 * alpha_I1_I2 * (I2 / K_I2) +
                alpha_S1_I3 * alpha_I1_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_I1_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_I2 * (I2 / K_I2)) * (
                alpha_S1_I3 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_S1_I4 * alpha_I2_I4 * (I4 / K_I4)
            ) +
            (alpha_S1_I3 * (I3 / K_I3)) * (alpha_S1_I4 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (S2 / K_S2) * (
            (S3 / K_S3) * (
                alpha_S2_P1 * alpha_S3_P1 * (P1 / K_P1) +
                alpha_S2_P2 * alpha_S3_P2 * (P2 / K_P2) +
                alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3) +
                alpha_S2_I1 * alpha_S3_I1 * (I1 / K_I1) +
                alpha_S2_I2 * alpha_S3_I2 * (I2 / K_I2) +
                alpha_S2_I3 * alpha_S3_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_S3_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_P1 * (P1 / K_P1)) * (
                alpha_S2_P2 * (P2 / K_P2) +
                alpha_S2_P3 * (P3 / K_P3) +
                alpha_S2_I1 * alpha_P1_I1 * (I1 / K_I1) +
                alpha_S2_I2 * alpha_P1_I2 * (I2 / K_I2) +
                alpha_S2_I3 * alpha_P1_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_P1_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_P2 * (P2 / K_P2)) * (
                alpha_S2_P3 * (P3 / K_P3) +
                alpha_S2_I1 * alpha_P2_I1 * (I1 / K_I1) +
                alpha_S2_I2 * alpha_P2_I2 * (I2 / K_I2) +
                alpha_S2_I3 * alpha_P2_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_P2_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_P3 * (P3 / K_P3)) * (
                alpha_S2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                alpha_S2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                alpha_S2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_P3_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_I1 * (I1 / K_I1)) * (
                alpha_S2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                alpha_S2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_I1_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_I2 * (I2 / K_I2)) * (
                alpha_S2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_S2_I4 * alpha_I2_I4 * (I4 / K_I4)
            ) +
            (alpha_S2_I3 * (I3 / K_I3)) * (alpha_S2_I4 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (S3 / K_S3) * (
            (alpha_S3_P1 * (P1 / K_P1)) * (
                alpha_S3_P2 * (P2 / K_P2) +
                alpha_S3_P3 * (P3 / K_P3) +
                alpha_S3_I1 * alpha_P1_I1 * (I1 / K_I1) +
                alpha_S3_I2 * alpha_P1_I2 * (I2 / K_I2) +
                alpha_S3_I3 * alpha_P1_I3 * (I3 / K_I3) +
                alpha_S3_I4 * alpha_P1_I4 * (I4 / K_I4)
            ) +
            (alpha_S3_P2 * (P2 / K_P2)) * (
                alpha_S3_P3 * (P3 / K_P3) +
                alpha_S3_I1 * alpha_P2_I1 * (I1 / K_I1) +
                alpha_S3_I2 * alpha_P2_I2 * (I2 / K_I2) +
                alpha_S3_I3 * alpha_P2_I3 * (I3 / K_I3) +
                alpha_S3_I4 * alpha_P2_I4 * (I4 / K_I4)
            ) +
            (alpha_S3_P3 * (P3 / K_P3)) * (
                alpha_S3_I1 * alpha_P3_I1 * (I1 / K_I1) +
                alpha_S3_I2 * alpha_P3_I2 * (I2 / K_I2) +
                alpha_S3_I3 * alpha_P3_I3 * (I3 / K_I3) +
                alpha_S3_I4 * alpha_P3_I4 * (I4 / K_I4)
            ) +
            (alpha_S3_I1 * (I1 / K_I1)) * (
                alpha_S3_I2 * alpha_S1_I2 * (I2 / K_I2) +
                alpha_S3_I3 * alpha_S1_I3 * (I3 / K_I3) +
                alpha_S3_I4 * alpha_S1_I4 * (I4 / K_I4)
            ) +
            (alpha_S3_I2 * (I2 / K_I2)) * (
                alpha_S3_I3 * alpha_S1_I3 * (I3 / K_I3) +
                alpha_S3_I4 * alpha_S1_I4 * (I4 / K_I4)
            ) +
            (alpha_S3_I3 * (I3 / K_I3)) * (alpha_S3_I4 * alpha_S1_I4 * (I4 / K_I4))
        ) +
        (P1 / K_P1) * (
            (P2 / K_P2) * (
                (P3 / K_P3) +
                alpha_P1_I1 * alpha_P2_I1 * (I1 / K_I1) +
                alpha_P1_I2 * alpha_P2_I2 * (I2 / K_I2) +
                alpha_P1_I3 * alpha_P2_I3 * (I3 / K_I3) +
                alpha_P1_I4 * alpha_P2_I4 * (I4 / K_I4)
            ) +
            (P3 / K_P3) * (
                alpha_P1_I1 * alpha_P3_I1 * (I1 / K_I1) +
                alpha_P1_I2 * alpha_P3_I2 * (I2 / K_I2) +
                alpha_P1_I3 * alpha_P3_I3 * (I3 / K_I3) +
                alpha_P1_I4 * alpha_P3_I4 * (I4 / K_I4)
            ) +
            (alpha_P1_I1 * (I1 / K_I1)) * (
                alpha_P2_I1 * alpha_I1_I2 * (I2 / K_I2) +
                alpha_P2_I2 * alpha_I1_I3 * (I3 / K_I3) +
                alpha_P2_I3 * alpha_I1_I4 * (I4 / K_I4)
            ) +
            (alpha_P1_I2 * (I2 / K_I2)) * (
                alpha_P2_I1 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_P2_I2 * alpha_I2_I4 * (I4 / K_I4)
            ) +
            (alpha_P1_I3 * (I3 / K_I3)) * (alpha_P2_I1 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (P2 / K_P2) * (
            (P3 / K_P3) * (
                alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                alpha_P2_I4 * alpha_P3_I4 * (I4 / K_I4)
            ) +
            (alpha_P2_I1 * (I1 / K_I1)) * (
                alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                alpha_P2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                alpha_P2_I4 * alpha_I1_I4 * (I4 / K_I4)
            ) +
            (alpha_P2_I2 * (I2 / K_I2)) * (
                alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_P2_I4 * alpha_I2_I4 * (I4 / K_I4)
            ) +
            (alpha_P2_I3 * (I3 / K_I3)) * (alpha_P2_I4 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (P3 / K_P3) * (
            (alpha_P3_I1 * (I1 / K_I1)) * (
                alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4)
            ) +
            (alpha_P3_I2 * (I2 / K_I2)) * (
                alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4)
            ) +
            (alpha_P3_I3 * (I3 / K_I3)) * (alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (I1 / K_I1) * (
            (alpha_I1_I2 * (I2 / K_I2)) * (
                alpha_I1_I3 * alpha_I2_I3 * (I3 / K_I3) +
                alpha_I1_I4 * alpha_I2_I4 * (I4 / K_I4)
            ) + (alpha_I1_I3 * (I3 / K_I3)) * (alpha_I1_I4 * alpha_I3_I4 * (I4 / K_I4))
        ) +
        (I2 / K_I2) *
        ((alpha_I2_I3 * (I3 / K_I3)) * (alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4)))

    four_metab_bound =
        (S1 / K_S1) * (
            (S2 / K_S2) * (
                (S3 / K_S3) * (
                    alpha_S1_P1 * alpha_S2_P1 * alpha_S3_P1 * (P1 / K_P1) +
                    alpha_S1_P2 * alpha_S2_P2 * alpha_S3_P2 * (P2 / K_P2) +
                    alpha_S1_P3 * alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3) +
                    alpha_S1_I1 * alpha_S2_I1 * alpha_S3_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S2_I2 * alpha_S3_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S2_I3 * alpha_S3_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_S3_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P1 * alpha_S2_P1 * (P1 / K_P1)) * (
                    alpha_S1_P2 * alpha_S2_P2 * (P2 / K_P2) +
                    alpha_S1_P3 * alpha_S2_P3 * (P3 / K_P3) +
                    alpha_S1_I1 * alpha_S2_I1 * alpha_P1_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S2_I2 * alpha_P1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S2_I3 * alpha_P1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_P1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P2 * alpha_S2_P2 * (P2 / K_P2)) * (
                    alpha_S1_P3 * alpha_S2_P3 * (P3 / K_P3) +
                    alpha_S1_I1 * alpha_S2_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S2_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S2_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P3 * alpha_S2_P3 * (P3 / K_P3)) * (
                    alpha_S1_I1 * alpha_S2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I1 * alpha_S2_I1 * (I1 / K_I1)) * (
                    alpha_S1_I2 * alpha_S2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I2 * alpha_S2_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_S2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S2_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * alpha_S2_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_S2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (S3 / K_S3) * (
                (alpha_S1_P1 * alpha_S3_P1 * (P1 / K_P1)) * (
                    alpha_S1_P2 * alpha_S3_P2 * (P2 / K_P2) +
                    alpha_S1_P3 * alpha_S3_P3 * (P3 / K_P3) +
                    alpha_S1_I1 * alpha_S3_I1 * alpha_P1_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S3_I2 * alpha_P1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S3_I3 * alpha_P1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S3_I4 * alpha_P1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P2 * alpha_S3_P2 * (P2 / K_P2)) * (
                    alpha_S1_P3 * alpha_S3_P3 * (P3 / K_P3) +
                    alpha_S1_I1 * alpha_S3_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S3_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S3_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S3_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P3 * alpha_S3_P3 * (P3 / K_P3)) * (
                    alpha_S1_I1 * alpha_S3_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_S3_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S3_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S3_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I1 * alpha_S3_I1 * (I1 / K_I1)) * (
                    alpha_S1_I2 * alpha_S3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_S3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S3_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I2 * alpha_S3_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_S3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_S3_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * alpha_S3_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_S3_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S1_P1 * (P1 / K_P1)) * (
                (alpha_S1_P2 * (P2 / K_P2)) * (
                    (alpha_S1_P3 * (P3 / K_P3)) +
                    alpha_S1_I1 * alpha_P1_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_P1_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P1_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P1_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_P3 * (P3 / K_P3)) * (
                    alpha_S1_I1 * alpha_P1_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_P1_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P1_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P1_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I1 * (I1 / K_I1)) * (
                    alpha_S1_I2 * alpha_P1_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P1_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P1_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_P1_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P1_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_P1_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S1_P2 * (P2 / K_P2)) * (
                (alpha_S1_P3 * (P3 / K_P3)) * (
                    alpha_S1_I1 * alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S1_I2 * alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P2_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I1 * alpha_P2_I1 * (I1 / K_I1)) * (
                    alpha_S1_I2 * alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I2 * alpha_P2_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P2_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * alpha_P2_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_P2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S1_P3 * (P3 / K_P3)) * (
                (alpha_S1_I1 * alpha_P3_I1 * (I1 / K_I1)) * (
                    alpha_S1_I2 * alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S1_I3 * alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I2 * alpha_P3_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * alpha_P3_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S1_I1 * (I1 / K_I1)) * (
                (alpha_S1_I2 * alpha_I1_I2 * (I2 / K_I2)) * (
                    alpha_S1_I3 * alpha_I2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S1_I4 * alpha_I2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S1_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S1_I2 * (I2 / K_I2)) * (
                (alpha_S1_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S1_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4))
            )
        ) +
        (S2 / K_S2) * (
            (S3 / K_S3) * (
                (alpha_S2_P1 * alpha_S3_P1 * (P1 / K_P1)) * (
                    alpha_S2_P2 * alpha_S3_P2 * (P2 / K_P2) +
                    alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3) +
                    alpha_S2_I1 * alpha_S3_I1 * alpha_P1_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_S3_I2 * alpha_P1_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_S3_I3 * alpha_P1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_S3_I4 * alpha_P1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_P2 * alpha_S3_P2 * (P2 / K_P2)) * (
                    alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3) +
                    alpha_S2_I1 * alpha_S3_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_S3_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_S3_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_S3_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_P3 * alpha_S3_P3 * (P3 / K_P3)) * (
                    alpha_S2_I1 * alpha_S3_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_S3_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_S3_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_S3_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I1 * alpha_S3_I1 * (I1 / K_I1)) * (
                    alpha_S2_I2 * alpha_S3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_S3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_S3_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I2 * alpha_S3_I2 * (I2 / K_I2)) * (
                    alpha_S2_I3 * alpha_S3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_S3_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I3 * alpha_S3_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_S3_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S2_P1 * (P1 / K_P1)) * (
                (alpha_S2_P2 * (P2 / K_P2)) * (
                    (alpha_S2_P3 * (P3 / K_P3)) +
                    alpha_S2_I1 * alpha_P1_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_P1_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P1_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P1_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_P3 * (P3 / K_P3)) * (
                    alpha_S2_I1 * alpha_P1_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_P1_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P1_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P1_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I1 * (I1 / K_I1)) * (
                    alpha_S2_I2 * alpha_P1_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P1_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P1_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I2 * (I2 / K_I2)) * (
                    alpha_S2_I3 * alpha_P1_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P1_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_P1_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S2_P2 * (P2 / K_P2)) * (
                (alpha_S2_P3 * (P3 / K_P3)) * (
                    alpha_S2_I1 * alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S2_I2 * alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P2_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I1 * alpha_P2_I1 * (I1 / K_I1)) * (
                    alpha_S2_I2 * alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I2 * alpha_P2_I2 * (I2 / K_I2)) * (
                    alpha_S2_I3 * alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P2_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I3 * alpha_P2_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_P2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S2_P3 * (P3 / K_P3)) * (
                (alpha_S2_I1 * alpha_P3_I1 * (I1 / K_I1)) * (
                    alpha_S2_I2 * alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S2_I3 * alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I2 * alpha_P3_I2 * (I2 / K_I2)) * (
                    alpha_S2_I3 * alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I3 * alpha_P3_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S2_I1 * (I1 / K_I1)) * (
                (alpha_S2_I2 * alpha_I1_I2 * (I2 / K_I2)) * (
                    alpha_S2_I3 * alpha_I2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S2_I4 * alpha_I2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S2_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S2_I2 * (I2 / K_I2)) * (
                (alpha_S2_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S2_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4))
            )
        ) +
        (S3 / K_S3) * (
            (alpha_S3_P1 * (P1 / K_P1)) * (
                (alpha_S3_P2 * (P2 / K_P2)) * (
                    (alpha_S3_P3 * (P3 / K_P3)) +
                    alpha_S3_I1 * alpha_P1_I1 * alpha_P2_I1 * (I1 / K_I1) +
                    alpha_S3_I2 * alpha_P1_I2 * alpha_P2_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P1_I3 * alpha_P2_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P1_I4 * alpha_P2_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_P3 * (P3 / K_P3)) * (
                    alpha_S3_I1 * alpha_P1_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S3_I2 * alpha_P1_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P1_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P1_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I1 * (I1 / K_I1)) * (
                    alpha_S3_I2 * alpha_P1_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P1_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P1_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I2 * (I2 / K_I2)) * (
                    alpha_S3_I3 * alpha_P1_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P1_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I3 * (I3 / K_I3)) *
                (alpha_S3_I4 * alpha_P1_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S3_P2 * (P2 / K_P2)) * (
                (alpha_S3_P3 * (P3 / K_P3)) * (
                    alpha_S3_I1 * alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1) +
                    alpha_S3_I2 * alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P2_I4 * alpha_P3_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I1 * alpha_P2_I1 * (I1 / K_I1)) * (
                    alpha_S3_I2 * alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I2 * alpha_P2_I2 * (I2 / K_I2)) * (
                    alpha_S3_I3 * alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P2_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I3 * alpha_P2_I3 * (I3 / K_I3)) *
                (alpha_S3_I4 * alpha_P2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S3_P3 * (P3 / K_P3)) * (
                (alpha_S3_I1 * alpha_P3_I1 * (I1 / K_I1)) * (
                    alpha_S3_I2 * alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2) +
                    alpha_S3_I3 * alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I2 * alpha_P3_I2 * (I2 / K_I2)) * (
                    alpha_S3_I3 * alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I3 * alpha_P3_I3 * (I3 / K_I3)) *
                (alpha_S3_I4 * alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S3_I1 * (I1 / K_I1)) * (
                (alpha_S3_I2 * alpha_I1_I2 * (I2 / K_I2)) * (
                    alpha_S3_I3 * alpha_I2_I3 * alpha_I1_I3 * (I3 / K_I3) +
                    alpha_S3_I4 * alpha_I2_I4 * alpha_I1_I4 * (I4 / K_I4)
                ) +
                (alpha_S3_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S3_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4))
            ) +
            (alpha_S3_I2 * (I2 / K_I2)) * (
                (alpha_S3_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                (alpha_S3_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4))
            )
        ) +
        (P1 / K_P1) * (
            (P2 / K_P2) * (
                (P3 / K_P3) * (
                    (alpha_P1_I1 * alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1)) +
                    (alpha_P1_I2 * alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2)) +
                    (alpha_P1_I3 * alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_P2_I4 * alpha_P3_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I1 * alpha_P2_I1 * (I1 / K_I1)) * (
                    (alpha_P1_I2 * alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2)) +
                    (alpha_P1_I3 * alpha_P2_I3 * alpha_I1_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_P2_I4 * alpha_I1_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I2 * alpha_P2_I2 * (I2 / K_I2)) * (
                    (alpha_P1_I3 * alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_P2_I4 * alpha_I2_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I3 * alpha_P2_I3 * (I3 / K_I3)) *
                ((alpha_P1_I4 * alpha_P2_I4 * alpha_I3_I4 * (I4 / K_I4)))
            ) +
            (P3 / K_P3) * (
                (alpha_P1_I1 * alpha_P3_I1 * (I1 / K_I1)) * (
                    (alpha_P1_I2 * alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2)) +
                    (alpha_P1_I3 * alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I2 * alpha_P3_I2 * (I2 / K_I2)) * (
                    (alpha_P1_I3 * alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I3 * alpha_P3_I3 * (I3 / K_I3)) *
                ((alpha_P1_I4 * alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4)))
            ) +
            (alpha_P1_I1 * (I1 / K_I1)) * (
                (alpha_P1_I2 * alpha_I1_I2 * (I2 / K_I2)) * (
                    (alpha_P1_I3 * alpha_I2_I3 * alpha_I1_I3 * (I3 / K_I3)) +
                    (alpha_P1_I4 * alpha_I2_I4 * alpha_I1_I4 * (I4 / K_I4))
                ) +
                (alpha_P1_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                ((alpha_P1_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4)))
            ) +
            (alpha_P1_I2 * (I2 / K_I2)) * (
                (alpha_P1_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                ((alpha_P1_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4)))
            )
        ) +
        (P2 / K_P2) * (
            (P3 / K_P3) * (
                (alpha_P2_I1 * alpha_P3_I1 * (I1 / K_I1)) * (
                    (alpha_P2_I2 * alpha_P3_I2 * alpha_I1_I2 * (I2 / K_I2)) +
                    (alpha_P2_I3 * alpha_P3_I3 * alpha_I1_I3 * (I3 / K_I3)) +
                    (alpha_P2_I4 * alpha_P3_I4 * alpha_I1_I4 * (I4 / K_I4))
                ) +
                (alpha_P2_I2 * alpha_P3_I2 * (I2 / K_I2)) * (
                    (alpha_P2_I3 * alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3)) +
                    (alpha_P2_I4 * alpha_P3_I4 * alpha_I2_I4 * (I4 / K_I4))
                ) +
                (alpha_P2_I3 * alpha_P3_I3 * (I3 / K_I3)) *
                ((alpha_P2_I4 * alpha_P3_I4 * alpha_I3_I4 * (I4 / K_I4)))
            ) +
            (alpha_P2_I1 * (I1 / K_I1)) * (
                (alpha_P2_I2 * alpha_I1_I2 * (I2 / K_I2)) * (
                    (alpha_P2_I3 * alpha_I2_I3 * alpha_I1_I3 * (I3 / K_I3)) +
                    (alpha_P2_I4 * alpha_I2_I4 * alpha_I1_I4 * (I4 / K_I4))
                ) +
                (alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                ((alpha_P2_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4)))
            ) +
            (alpha_P2_I2 * (I2 / K_I2)) * (
                (alpha_P2_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                ((alpha_P2_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4)))
            )
        ) +
        (P3 / K_P3) * (
            (alpha_P3_I1 * (I1 / K_I1)) * (
                (alpha_P3_I2 * alpha_I1_I2 *(I2 / K_I2)) * (
                    (alpha_P3_I3 * alpha_I1_I3 * alpha_I2_I3 * (I3 / K_I3)) +
                    (alpha_P3_I4 * alpha_I1_I4 * alpha_I2_I4 * (I4 / K_I4))
                )
            ) +
            (alpha_P3_I2 * (I2 / K_I2)) * (
                (alpha_P3_I3 * alpha_I2_I3 * (I3 / K_I3)) *
                ((alpha_P3_I4 * alpha_I3_I4 * alpha_I2_I4 * (I4 / K_I4)))
            )
        ) +
        (I1 / K_I1) * (
            alpha_I1_I2 *
            (I2 / K_I2) *
            (alpha_I1_I3 * alpha_I2_I3 * (I3 / K_I3)) *
            (alpha_I1_I4 * alpha_I2_I4 * alpha_I3_I4 * (I4 / K_I4))
        )

    return 1 + one_metab_bound + two_metab_bound + three_metab_bound + four_metab_bound
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
    alpha_R1_R2,
    alpha_R1_R3,
    alpha_R1_R4,
    alpha_R1_R5,
    alpha_R1_R6,
    alpha_R2_R3,
    alpha_R2_R4,
    alpha_R2_R5,
    alpha_R2_R6,
    alpha_R3_R4,
    alpha_R3_R5,
    alpha_R3_R6,
    alpha_R4_R5,
    alpha_R4_R6,
    alpha_R5_R6,
)

    one_metab_bound =
        (R1 / K_R1) + (R2 / K_R2) + (R3 / K_R3) + (R4 / K_R4) + (R5 / K_R5) + (R6 / K_R6)
    two_metab_bound =
        (R1 / K_R1) * (
            alpha_R1_R2 * (R2 / K_R2) +
            alpha_R1_R3 * (R3 / K_R3) +
            alpha_R1_R4 * (R4 / K_R4) +
            alpha_R1_R5 * (R5 / K_R5) +
            alpha_R1_R6 * (R6 / K_R6)
        ) +
        (R2 / K_R2) * (
            alpha_R2_R3 * (R3 / K_R3) +
            alpha_R2_R4 * (R4 / K_R4) +
            alpha_R2_R5 * (R5 / K_R5) +
            alpha_R2_R6 * (R6 / K_R6)
        ) +
        (R3 / K_R3) * (
            alpha_R3_R4 * (R4 / K_R4) +
            alpha_R3_R5 * (R3 / K_R3) * (R5 / K_R5) +
            alpha_R3_R6 * (R3 / K_R3) * (R6 / K_R6)
        ) +
        (R4 / K_R4) *
        (alpha_R4_R5 * (R5 / K_R5) + alpha_R4_R6 * (R4 / K_R4) * (R6 / K_R6)) +
        (R5 / K_R5) * (alpha_R5_R6 * (R6 / K_R6))

    three_metab_bound =
        (R1 / K_R1) * (
            (alpha_R1_R2 * (R2 / K_R2)) * (
                alpha_R1_R3 * alpha_R2_R3 * (R3 / K_R3) +
                alpha_R1_R4 * alpha_R2_R4 * (R4 / K_R4) +
                alpha_R1_R5 * alpha_R2_R5 * (R5 / K_R5) +
                alpha_R1_R6 * alpha_R2_R6 * (R6 / K_R6)
            ) +
            (alpha_R1_R3 * (R3 / K_R3)) * (
                alpha_R1_R4 * alpha_R3_R4 * (R4 / K_R4) +
                alpha_R1_R5 * alpha_R3_R5 * (R5 / K_R5) +
                alpha_R1_R6 * alpha_R3_R6 * (R6 / K_R6)
            ) +
            (alpha_R1_R4 * (R4 / K_R4)) * (
                alpha_R1_R5 * alpha_R4_R5 * (R5 / K_R5) +
                alpha_R1_R6 * alpha_R4_R6 * (R6 / K_R6)
            ) +
            (alpha_R1_R5 * (R5 / K_R5)) * (alpha_R1_R6 * alpha_R5_R6 * (R6 / K_R6))
        ) +
        (R2 / K_R2) * (
            (alpha_R2_R3 * (R3 / K_R3)) * (
                alpha_R2_R4 * alpha_R3_R4 * (R4 / K_R4) +
                alpha_R2_R5 * alpha_R3_R5 * (R5 / K_R5) +
                alpha_R2_R6 * alpha_R3_R6 * (R6 / K_R6)
            ) +
            (alpha_R2_R4 * (R4 / K_R4)) * (
                alpha_R2_R5 * alpha_R4_R5 * (R5 / K_R5) +
                alpha_R2_R6 * alpha_R4_R6 * (R6 / K_R6)
            ) +
            (alpha_R2_R5 * (R5 / K_R5)) * (alpha_R2_R6 * alpha_R5_R6 * (R6 / K_R6))
        ) +
        (R3 / K_R3) * (
            (alpha_R3_R4 * (R4 / K_R4)) * (
                alpha_R3_R5 * alpha_R4_R5 * (R5 / K_R5) +
                alpha_R3_R6 * alpha_R4_R6 * (R6 / K_R6)
            ) + (alpha_R3_R5 * (R5 / K_R5)) * (alpha_R3_R6 * alpha_R5_R6 * (R6 / K_R6))
        ) +
        (R4 / K_R4) *
        ((alpha_R4_R5 * (R5 / K_R5)) * (alpha_R4_R6 * alpha_R5_R6 * (R6 / K_R6)))

    four_metab_bound =
        (R1 / K_R1) * (
            (alpha_R1_R2 * (R2 / K_R2)) * (
                (alpha_R1_R3 * alpha_R2_R3 * (R3 / K_R3)) * (
                    alpha_R1_R4 * alpha_R2_R4 * alpha_R3_R4 * (R4 / K_R4) +
                    alpha_R1_R5 * alpha_R2_R5 * alpha_R3_R5 * (R5 / K_R5) +
                    alpha_R1_R6 * alpha_R2_R6 * alpha_R3_R6 * (R6 / K_R6)
                ) +
                (alpha_R1_R4 * alpha_R2_R4 * (R4 / K_R4)) * (
                    alpha_R1_R5 * alpha_R2_R5 * alpha_R4_R5 * (R5 / K_R5) +
                    alpha_R1_R6 * alpha_R2_R6 * alpha_R4_R6 * (R6 / K_R6)
                ) +
                (alpha_R1_R5 * alpha_R2_R5 * (R5 / K_R5)) *
                (alpha_R1_R6 * alpha_R2_R6 * alpha_R5_R6 * (R6 / K_R6))
            ) +
            (alpha_R1_R3 * (R3 / K_R3)) * (
                (alpha_R1_R4 * alpha_R3_R4 * (R4 / K_R4)) * (
                    alpha_R1_R5 * alpha_R3_R5 * alpha_R4_R5 * (R5 / K_R5) +
                    alpha_R1_R6 * alpha_R3_R6 * alpha_R4_R6 * (R6 / K_R6)
                ) +
                (alpha_R1_R5 * alpha_R3_R5 * (R5 / K_R5)) *
                (alpha_R1_R6 * alpha_R3_R6 * alpha_R5_R6 * (R6 / K_R6))
            ) +
            (alpha_R1_R4 * (R4 / K_R4)) * (
                (alpha_R1_R5 * alpha_R4_R5 * (R5 / K_R5)) *
                (alpha_R1_R6 * alpha_R4_R6 * alpha_R5_R6 * (R6 / K_R6))
            )
        ) +
        (R2 / K_R2) * (
            (alpha_R2_R3 * (R3 / K_R3)) * (
                (alpha_R2_R4 * alpha_R3_R4 * (R4 / K_R4)) * (
                    alpha_R2_R5 * alpha_R3_R5 * alpha_R4_R5 * (R5 / K_R5) +
                    alpha_R2_R6 * alpha_R3_R6 * alpha_R4_R6 * (R6 / K_R6)
                ) +
                (alpha_R2_R5 * alpha_R3_R5 * (R5 / K_R5)) *
                (alpha_R2_R6 * alpha_R3_R6 * alpha_R5_R6 * (R6 / K_R6))
            ) +
            (alpha_R2_R4 * (R4 / K_R4)) * (
                (alpha_R2_R5 * alpha_R4_R5 * (R5 / K_R5)) *
                (alpha_R2_R6 * alpha_R4_R6 * alpha_R5_R6 * (R6 / K_R6))
            )
        ) +
        (R3 / K_R3) * (
            alpha_R3_R4 *
            (R4 / K_R4) *
            (
                (alpha_R3_R5 * alpha_R4_R5 * (R5 / K_R5)) *
                (alpha_R3_R6 * alpha_R4_R6 * alpha_R5_R6 * (R6 / K_R6))
            )
        )

    five_metab_bound =
        alpha_R1_R2 *
        alpha_R1_R3 *
        alpha_R1_R4 *
        alpha_R1_R5 *
        alpha_R2_R3 *
        alpha_R2_R4 *
        alpha_R2_R5 *
        alpha_R3_R4 *
        alpha_R3_R5 *
        alpha_R4_R5 *
        (R1 / K_R1) *
        (R2 / K_R2) *
        (R3 / K_R3) *
        (R4 / K_R4) *
        (R5 / K_R5) +
        alpha_R1_R2 *
        alpha_R1_R3 *
        alpha_R1_R4 *
        alpha_R1_R6 *
        alpha_R2_R3 *
        alpha_R2_R4 *
        alpha_R2_R6 *
        alpha_R3_R4 *
        alpha_R3_R6 *
        alpha_R4_R6 *
        (R1 / K_R1) *
        (R2 / K_R2) *
        (R3 / K_R3) *
        (R4 / K_R4) *
        (R6 / K_R6) +
        alpha_R1_R2 *
        alpha_R1_R3 *
        alpha_R1_R5 *
        alpha_R1_R6 *
        alpha_R2_R3 *
        alpha_R2_R5 *
        alpha_R2_R6 *
        alpha_R3_R5 *
        alpha_R3_R6 *
        alpha_R5_R6 *
        (R1 / K_R1) *
        (R2 / K_R2) *
        (R3 / K_R3) *
        (R5 / K_R5) *
        (R6 / K_R6) +
        alpha_R1_R2 *
        alpha_R1_R4 *
        alpha_R1_R5 *
        alpha_R1_R6 *
        alpha_R2_R4 *
        alpha_R2_R5 *
        alpha_R2_R6 *
        alpha_R4_R5 *
        alpha_R4_R6 *
        alpha_R5_R6 *
        (R1 / K_R1) *
        (R2 / K_R2) *
        (R4 / K_R4) *
        (R5 / K_R5) *
        (R6 / K_R6) +
        alpha_R1_R3 *
        alpha_R1_R4 *
        alpha_R1_R5 *
        alpha_R1_R6 *
        alpha_R3_R4 *
        alpha_R3_R5 *
        alpha_R3_R6 *
        alpha_R4_R5 *
        alpha_R4_R6 *
        alpha_R5_R6 *
        (R1 / K_R1) *
        (R3 / K_R3) *
        (R4 / K_R4) *
        (R5 / K_R5) *
        (R6 / K_R6) +
        alpha_R2_R3 *
        alpha_R2_R4 *
        alpha_R2_R5 *
        alpha_R2_R6 *
        alpha_R3_R4 *
        alpha_R3_R5 *
        alpha_R3_R6 *
        alpha_R4_R5 *
        alpha_R4_R6 *
        alpha_R5_R6 *
        (R2 / K_R2) *
        (R3 / K_R3) *
        (R4 / K_R4) *
        (R5 / K_R5) *
        (R6 / K_R6)

    six_metab_bound =
        alpha_R1_R2 *
        alpha_R1_R3 *
        alpha_R1_R4 *
        alpha_R1_R5 *
        alpha_R1_R6 *
        alpha_R2_R3 *
        alpha_R2_R4 *
        alpha_R2_R5 *
        alpha_R2_R6 *
        alpha_R3_R4 *
        alpha_R3_R5 *
        alpha_R3_R6 *
        alpha_R4_R5 *
        alpha_R4_R6 *
        alpha_R5_R6 *
        (R1 / K_R1) *
        (R2 / K_R2) *
        (R3 / K_R3) *
        (R4 / K_R4) *
        (R5 / K_R5) *
        (R6 / K_R6)

    return 1 +
           one_metab_bound +
           two_metab_bound +
           three_metab_bound +
           four_metab_bound +
           five_metab_bound +
           six_metab_bound
end

"Generate the names of the parameters for the rate equation using the same input as @derive_general_mwc_rate_eq"
function generate_param_names(processed_input)
    param_names = (:L, :Vmax_a, :Vmax_i)
    for metab in [processed_input[:substrates]..., processed_input[:products]...]
        param_names = (param_names..., Symbol("K_a_", metab), Symbol("K_i_", metab))
    end
    for metab in processed_input[:regulators]
        param_names =
            (param_names..., Symbol("K_a_", metab, "_reg"), Symbol("K_i_", metab, "_reg"))
    end
    for substrate in processed_input[:substrates]
        for metabolite in [processed_input[:products]..., processed_input[:regulators]...]
            param_names = (param_names..., Symbol("alpha_", substrate, "_", metabolite))
        end
    end
    for product in processed_input[:products]
        for regulator in processed_input[:regulators]
            param_names = (param_names..., Symbol("alpha_", product, "_", regulator))
        end
    end
    for (i, regulator1) in enumerate(processed_input[:regulators][1:(end-1)])
        for regulator2 in processed_input[:regulators][(i+1):end]
            param_names = (param_names..., Symbol("alpha_", regulator1, "_", regulator2))
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
