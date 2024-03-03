#=
CODE FOR RATE EQUATION DERIVATION
=#

@doc raw"""
    derive_general_mwc_rate_eq(metabs_and_regulators_kwargs...)

Derive a function that calculates the rate of a reaction using the general MWC rate equation given the list of substrates, products, and regulators that bind to specific cat or reg sites.

The general MWC rate equation is given by:

```math
Rate = \frac{{V_{max}^a \prod_{i=1}^{n} \left(\frac{S_i}{K_{a, i}}\right) - V_{max}^a_{rev} \prod_{i=1}^{n} \left(\frac{P_i}{K_{a, i}}\right) \cdot Z_{a, cat}^{n-1} \cdot Z_{a, reg}^n + L \left(V_{max}^i \prod_{i=1}^{n} \left(\frac{S_i}{K_{i, i}}\right) - V_{max}^i_{rev} \prod_{i=1}^{n} \left(\frac{P_i}{K_{i, i}}\right)\right) \cdot Z_{i, cat}^{n-1} \cdot Z_{i, reg}^n}}{Z_{a, cat}^n \cdot Z_{a, reg}^n + L \cdot Z_{i, cat}^n \cdot Z_{i, reg}^n}
```

where:
- ``V_{max}^a`` is the maximum rate of the forward reaction
- ``V_{max}^a_{rev}`` is the maximum rate of the reverse reaction
- ``V_{max}^i`` is the maximum rate of the forward reaction
- ``V_{max}^i_{rev}`` is the maximum rate of the reverse reaction
- ``S_i`` is the concentration of the ``i^{th}`` substrate
- ``P_i`` is the concentration of the ``i^{th}`` product
- ``K_{a, i}`` is the Michaelis constant for the ``i^{th}`` substrate
- ``K_{i, i}`` is the Michaelis constant for the ``i^{th}`` product
- ``Z_{a, cat}`` is the allosteric factor for the catalytic site
- ``Z_{i, cat}`` is the allosteric factor for the catalytic site
- ``Z_{a, reg}`` is the allosteric factor for the regulatory site
- ``Z_{i, reg}`` is the allosteric factor for the regulatory site
- ``L`` is the concentration of the ligand
- ``n`` is the oligomeric state of the enzyme

# Arguments
- `metabs_and_regulators_kwargs...`: keyword arguments that specify the substrates, products, catalytic sites, regulatory sites, and other parameters of the reaction.

# Returns
- A function that calculates the rate of the reaction using the general MWC rate equation
- A tuple of the names of the metabolites and parameters used in the rate equation
"""
macro derive_general_mwc_rate_eq(metabs_and_regulators_kwargs)

    # processed_input = NamedTuple()
    # for expr in metabs_and_regulators_kwargs
    #     expr.args[1] ∈ expected_input_kwargs || error(
    #         "invalid keyword: ",
    #         expr.args[1],
    #         ". The only supported keywords are: ",
    #         expected_input_kwargs,
    #     )
    #     println(expr)
    #     println(expr.args[1])
    #     println(typeof(expr.args[1]))
    #     println(expr.args[2])
    #     println(typeof(expr.args[2]))
    #     processed_input = merge(processed_input, (; expr.args[1] => eval(expr)))
    # end
    processed_input = getfield(__module__, Symbol(metabs_and_regulators_kwargs))
    expected_input_kwargs = [
        :substrates,
        :products,
        :cat1,
        :cat2,
        :cat3,
        :cat4,
        :reg1,
        :reg2,
        :reg3,
        :Keq,
        :oligomeric_state,
        :rate_equation_name
    ]
    for kwarg in keys(processed_input)
        kwarg ∈ expected_input_kwargs || error(
            "invalid keyword: ",
            kwarg,
            ". The only supported keywords are: ",
            expected_input_kwargs,
        )
    end
    @assert 0 < length(processed_input.substrates) <= 4 "At least 1 and no more that 4 substrates are supported"
    @assert 0 < length(processed_input.products) <= 4 "At least 1 and no more that 4 products are supported"
    cat_site_metabs = [processed_input.cat1...]
    for field in [:cat2, :cat3, :cat4]
        if hasproperty(processed_input, field)
            push!(cat_site_metabs, processed_input[field]...)
        end
    end
    # println(cat_site_metabs)
    for metab_name in [processed_input.substrates..., processed_input.products...]
        @assert metab_name ∈ cat_site_metabs "Each substrate and product has to be assigned to one of the catalytic sites"
    end
    for field in [:cat1, :cat2, :cat3, :cat4]
        hasproperty(processed_input, field) &&
            @assert length(processed_input[field]) <= 2 "No more that 2 metabolites can bind to one catalytic site. One substrate and one product"
    end
    for field in [:reg1, :reg2, :reg3]
        hasproperty(processed_input, field) &&
            @assert length(processed_input[field]) <= 3 "No more that 3 catalytic sites are supported"
    end
    metab_names = generate_metab_names(processed_input)
    param_names = generate_param_names(processed_input)

    enz = NamedTuple()
    for field in propertynames(processed_input)
        if field == :cat1
            for metab_name in processed_input[field]
                if metab_name ∈ processed_input.substrates
                    enz = merge(enz, (; Symbol(:S, 1, "_", field) => metab_name))
                elseif metab_name ∈ processed_input.products
                    enz = merge(enz, (; Symbol(:P, 1, "_", field) => metab_name))
                end
            end
        elseif hasproperty(processed_input, :cat2) && field == :cat2
            for metab_name in processed_input[field]
                if metab_name ∈ processed_input.substrates
                    enz = merge(enz, (; Symbol(:S, 2, "_", field) => metab_name))
                elseif metab_name ∈ processed_input.products
                    enz = merge(enz, (; Symbol(:P, 2, "_", field) => metab_name))
                end
            end
        elseif hasproperty(processed_input, :cat3) && field == :cat3
            for metab_name in processed_input[field]
                if metab_name ∈ processed_input.substrates
                    enz = merge(enz, (; Symbol(:S, 3, "_", field) => metab_name))
                elseif metab_name ∈ processed_input.products
                    enz = merge(enz, (; Symbol(:P, 3, "_", field) => metab_name))
                end
            end
        elseif hasproperty(processed_input, :cat4) && field == :cat4
            for metab_name in processed_input[field]
                if metab_name ∈ processed_input.substrates
                    enz = merge(enz, (; Symbol(:S, 4, "_", field) => metab_name))
                elseif metab_name ∈ processed_input.products
                    enz = merge(enz, (; Symbol(:P, 4, "_", field) => metab_name))
                end
            end
        elseif hasproperty(processed_input, :reg1) && field == :reg1
            for (i, metab_name) in enumerate(processed_input[field])
                enz = merge(enz, (; Symbol(:R, i, "_", field) => metab_name))
            end
        elseif hasproperty(processed_input, :reg2) && field == :reg2
            for (i, metab_name) in enumerate(processed_input[field])
                enz = merge(enz, (; Symbol(:R, i, "_", field) => metab_name))
            end
        elseif hasproperty(processed_input, :reg3) && field == :reg3
            for (i, metab_name) in enumerate(processed_input[field])
                enz = merge(enz, (; Symbol(:R, i, "_", field) => metab_name))
            end
        end
    end

    #TODO: use Base.method_argnames(methods(general_mwc_rate_equation)[1])[2:end] to get args
    mwc_rate_eq_args = [
        :S1_cat1,
        :S2_cat2,
        :S3_cat3,
        :S4_cat4,
        :P1_cat1,
        :P2_cat2,
        :P3_cat3,
        :P4_cat4,
        :R1_reg1,
        :R2_reg1,
        :R3_reg1,
        :R1_reg2,
        :R2_reg2,
        :R3_reg2,
        :R1_reg3,
        :R2_reg3,
        :R3_reg3,
    ]
    missing_keys = filter(x -> !haskey(enz, x), mwc_rate_eq_args)
    for key in missing_keys
        enz = merge(enz, (; key => nothing))
    end
    # qualified_name = esc(GlobalRef(Main, :rate_equation))
    function_name = (hasproperty(processed_input, :rate_equation_name) ? esc(processed_input.rate_equation_name) : esc(:rate_equation))
    oligomeric_state = esc(processed_input.oligomeric_state)
    return quote
        @inline function $(function_name)(metabs, params, Keq)
            general_mwc_rate_equation(
                $(enz.S1_cat1 isa Symbol) ? metabs.$(enz.S1_cat1) : 1.0,
                $(enz.S2_cat2 isa Symbol) ? metabs.$(enz.S2_cat2) : 1.0,
                $(enz.S3_cat3 isa Symbol) ? metabs.$(enz.S3_cat3) : 1.0,
                $(enz.S4_cat4 isa Symbol) ? metabs.$(enz.S4_cat4) : 1.0,
                $(enz.P1_cat1 isa Symbol) ? metabs.$(enz.P1_cat1) : 1.0,
                $(enz.P2_cat2 isa Symbol) ? metabs.$(enz.P2_cat2) : 1.0,
                $(enz.P3_cat3 isa Symbol) ? metabs.$(enz.P3_cat3) : 1.0,
                $(enz.P4_cat4 isa Symbol) ? metabs.$(enz.P4_cat4) : 1.0,
                params.L,
                params.Vmax_a,
                params.Vmax_i,
                $(enz.S1_cat1 isa Symbol) ? params.$(Symbol("K_a_", enz.S1_cat1, "_cat1")) :
                1.0,
                $(enz.S1_cat1 isa Symbol) ? params.$(Symbol("K_i_", enz.S1_cat1, "_cat1")) :
                1.0,
                $(enz.S2_cat2 isa Symbol) ? params.$(Symbol("K_a_", enz.S2_cat2, "_cat2")) :
                1.0,
                $(enz.S2_cat2 isa Symbol) ? params.$(Symbol("K_i_", enz.S2_cat2, "_cat2")) :
                1.0,
                $(enz.S3_cat3 isa Symbol) ? params.$(Symbol("K_a_", enz.S3_cat3, "_cat3")) :
                1.0,
                $(enz.S3_cat3 isa Symbol) ? params.$(Symbol("K_i_", enz.S3_cat3, "_cat3")) :
                1.0,
                $(enz.S4_cat4 isa Symbol) ? params.$(Symbol("K_a_", enz.S4_cat4, "_cat4")) :
                1.0,
                $(enz.S4_cat4 isa Symbol) ? params.$(Symbol("K_i_", enz.S4_cat4, "_cat4")) :
                1.0,
                $(enz.P1_cat1 isa Symbol) ? params.$(Symbol("K_a_", enz.P1_cat1, "_cat1")) :
                1.0,
                $(enz.P1_cat1 isa Symbol) ? params.$(Symbol("K_i_", enz.P1_cat1, "_cat1")) :
                1.0,
                $(enz.P2_cat2 isa Symbol) ? params.$(Symbol("K_a_", enz.P2_cat2, "_cat2")) :
                1.0,
                $(enz.P2_cat2 isa Symbol) ? params.$(Symbol("K_i_", enz.P2_cat2, "_cat2")) :
                1.0,
                $(enz.P3_cat3 isa Symbol) ? params.$(Symbol("K_a_", enz.P3_cat3, "_cat3")) :
                1.0,
                $(enz.P3_cat3 isa Symbol) ? params.$(Symbol("K_i_", enz.P3_cat3, "_cat3")) :
                1.0,
                $(enz.P4_cat4 isa Symbol) ? params.$(Symbol("K_a_", enz.P4_cat4, "_cat4")) :
                1.0,
                $(enz.P4_cat4 isa Symbol) ? params.$(Symbol("K_i_", enz.P4_cat4, "_cat4")) :
                1.0,
                #Z_a_cat
                calculate_z_cat(
                    $(enz.S1_cat1 isa Symbol) ? metabs.$(enz.S1_cat1) : 0.0,
                    $(enz.S2_cat2 isa Symbol) ? metabs.$(enz.S2_cat2) : 0.0,
                    $(enz.S3_cat3 isa Symbol) ? metabs.$(enz.S3_cat3) : 0.0,
                    $(enz.S4_cat4 isa Symbol) ? metabs.$(enz.S4_cat4) : 0.0,
                    $(enz.P1_cat1 isa Symbol) ? metabs.$(enz.P1_cat1) : 0.0,
                    $(enz.P2_cat2 isa Symbol) ? metabs.$(enz.P2_cat2) : 0.0,
                    $(enz.P3_cat3 isa Symbol) ? metabs.$(enz.P3_cat3) : 0.0,
                    $(enz.P4_cat4 isa Symbol) ? metabs.$(enz.P4_cat4) : 0.0,
                    $(enz.S1_cat1 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.S1_cat1, "_cat1")) : Inf,
                    $(enz.S2_cat2 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.S2_cat2, "_cat2")) : Inf,
                    $(enz.S3_cat3 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.S3_cat3, "_cat3")) : Inf,
                    $(enz.S4_cat4 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.S4_cat4, "_cat4")) : Inf,
                    $(enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.P1_cat1, "_cat1")) : Inf,
                    $(enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.P2_cat2, "_cat2")) : Inf,
                    $(enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.P3_cat3, "_cat3")) : Inf,
                    $(enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.P4_cat4, "_cat4")) : Inf,
                    $(enz.S1_cat1 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S1_cat1 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P3_cat3)) : 0.0,
                    $(enz.S1_cat1 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P3_cat3)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P3_cat3)) : 0.0,
                ),
                #Z_i_cat
                calculate_z_cat(
                    $(enz.S1_cat1 isa Symbol) ? metabs.$(enz.S1_cat1) : 0.0,
                    $(enz.S2_cat2 isa Symbol) ? metabs.$(enz.S2_cat2) : 0.0,
                    $(enz.S3_cat3 isa Symbol) ? metabs.$(enz.S3_cat3) : 0.0,
                    $(enz.S4_cat4 isa Symbol) ? metabs.$(enz.S4_cat4) : 0.0,
                    $(enz.P1_cat1 isa Symbol) ? metabs.$(enz.P1_cat1) : 0.0,
                    $(enz.P2_cat2 isa Symbol) ? metabs.$(enz.P2_cat2) : 0.0,
                    $(enz.P3_cat3 isa Symbol) ? metabs.$(enz.P3_cat3) : 0.0,
                    $(enz.P4_cat4 isa Symbol) ? metabs.$(enz.P4_cat4) : 0.0,
                    $(enz.S1_cat1 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.S1_cat1, "_cat1")) : Inf,
                    $(enz.S2_cat2 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.S2_cat2, "_cat2")) : Inf,
                    $(enz.S3_cat3 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.S3_cat3, "_cat3")) : Inf,
                    $(enz.S4_cat4 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.S4_cat4, "_cat4")) : Inf,
                    $(enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.P1_cat1, "_cat1")) : Inf,
                    $(enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.P2_cat2, "_cat2")) : Inf,
                    $(enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.P3_cat3, "_cat3")) : Inf,
                    $(enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.P4_cat4, "_cat4")) : Inf,
                    $(enz.S1_cat1 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S1_cat1 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P3_cat3)) : 0.0,
                    $(enz.S1_cat1 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S1_cat1, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P3_cat3)) : 0.0,
                    $(enz.S2_cat2 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S2_cat2, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S3_cat3 isa Symbol && enz.P4_cat4 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S3_cat3, "_", enz.P4_cat4)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P1_cat1 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P1_cat1)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P2_cat2 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P2_cat2)) : 0.0,
                    $(enz.S4_cat4 isa Symbol && enz.P3_cat3 isa Symbol) ?
                    params.$(Symbol("alpha_", enz.S4_cat4, "_", enz.P3_cat3)) : 0.0,
                ),
                #Z_a_reg
                calculate_z_reg(
                    $(enz.R1_reg1 isa Symbol) ? metabs.$(enz.R1_reg1) : 0.0,
                    $(enz.R2_reg1 isa Symbol) ? metabs.$(enz.R2_reg1) : 0.0,
                    $(enz.R3_reg1 isa Symbol) ? metabs.$(enz.R3_reg1) : 0.0,
                    $(enz.R1_reg2 isa Symbol) ? metabs.$(enz.R1_reg2) : 0.0,
                    $(enz.R2_reg2 isa Symbol) ? metabs.$(enz.R2_reg2) : 0.0,
                    $(enz.R3_reg2 isa Symbol) ? metabs.$(enz.R3_reg2) : 0.0,
                    $(enz.R1_reg3 isa Symbol) ? metabs.$(enz.R1_reg3) : 0.0,
                    $(enz.R2_reg3 isa Symbol) ? metabs.$(enz.R2_reg3) : 0.0,
                    $(enz.R3_reg3 isa Symbol) ? metabs.$(enz.R3_reg3) : 0.0,
                    $(enz.R1_reg1 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R1_reg1, "_reg1")) : Inf,
                    $(enz.R2_reg1 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R2_reg1, "_reg1")) : Inf,
                    $(enz.R3_reg1 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R3_reg1, "_reg1")) : Inf,
                    $(enz.R1_reg2 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R1_reg2, "_reg2")) : Inf,
                    $(enz.R2_reg2 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R2_reg2, "_reg2")) : Inf,
                    $(enz.R3_reg2 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R3_reg2, "_reg2")) : Inf,
                    $(enz.R1_reg3 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R1_reg3, "_reg3")) : Inf,
                    $(enz.R2_reg3 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R2_reg3, "_reg3")) : Inf,
                    $(enz.R3_reg3 isa Symbol) ?
                    params.$(Symbol("K_a_", enz.R3_reg3, "_reg3")) : Inf,
                ),
                #Z_i_reg
                calculate_z_reg(
                    $(enz.R1_reg1 isa Symbol) ? metabs.$(enz.R1_reg1) : 0.0,
                    $(enz.R2_reg1 isa Symbol) ? metabs.$(enz.R2_reg1) : 0.0,
                    $(enz.R3_reg1 isa Symbol) ? metabs.$(enz.R3_reg1) : 0.0,
                    $(enz.R1_reg2 isa Symbol) ? metabs.$(enz.R1_reg2) : 0.0,
                    $(enz.R2_reg2 isa Symbol) ? metabs.$(enz.R2_reg2) : 0.0,
                    $(enz.R3_reg2 isa Symbol) ? metabs.$(enz.R3_reg2) : 0.0,
                    $(enz.R1_reg3 isa Symbol) ? metabs.$(enz.R1_reg3) : 0.0,
                    $(enz.R2_reg3 isa Symbol) ? metabs.$(enz.R2_reg3) : 0.0,
                    $(enz.R3_reg3 isa Symbol) ? metabs.$(enz.R3_reg3) : 0.0,
                    $(enz.R1_reg1 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R1_reg1, "_reg1")) : Inf,
                    $(enz.R2_reg1 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R2_reg1, "_reg1")) : Inf,
                    $(enz.R3_reg1 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R3_reg1, "_reg1")) : Inf,
                    $(enz.R1_reg2 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R1_reg2, "_reg2")) : Inf,
                    $(enz.R2_reg2 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R2_reg2, "_reg2")) : Inf,
                    $(enz.R3_reg2 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R3_reg2, "_reg2")) : Inf,
                    $(enz.R1_reg3 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R1_reg3, "_reg3")) : Inf,
                    $(enz.R2_reg3 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R2_reg3, "_reg3")) : Inf,
                    $(enz.R3_reg3 isa Symbol) ?
                    params.$(Symbol("K_i_", enz.R3_reg3, "_reg3")) : Inf,
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
    S1_cat1::T,
    S2_cat2::T,
    S3_cat3::T,
    S4_cat4::T,
    P1_cat1::T,
    P2_cat2::T,
    P3_cat3::T,
    P4_cat4::T,
    L::T,
    Vmax_a::T,
    Vmax_i::T,
    K_a_S1_cat1::T,
    K_i_S1_cat1::T,
    K_a_S2_cat2::T,
    K_i_S2_cat2::T,
    K_a_S3_cat3::T,
    K_i_S3_cat3::T,
    K_a_S4_cat4::T,
    K_i_S4_cat4::T,
    K_a_P1_cat1::T,
    K_i_P1_cat1::T,
    K_a_P2_cat2::T,
    K_i_P2_cat2::T,
    K_a_P3_cat3::T,
    K_i_P3_cat3::T,
    K_a_P4_cat4::T,
    K_i_P4_cat4::T,
    Z_a_cat::T,
    Z_i_cat::T,
    Z_a_reg::T,
    Z_i_reg::T,
    Keq::T,
    n::Int,
) where {T<:Float64}
    Vmax_a = 1.0
    Vmax_a_rev = ifelse(
        (K_a_P1_cat1 * K_a_P2_cat2 * K_a_P3_cat3 * K_a_P4_cat4) != Inf,
        Vmax_a * K_a_P1_cat1 * K_a_P2_cat2 * K_a_P3_cat3 * K_a_P4_cat4 /
        (Keq * K_a_S1_cat1 * K_a_S2_cat2 * K_a_S3_cat3 * K_a_S4_cat4),
        0.0,
    )
    Vmax_i_rev = ifelse(
        (K_i_P1_cat1 * K_i_P2_cat2 * K_i_P3_cat3 * K_i_S4_cat4) != Inf,
        Vmax_i * K_i_P1_cat1 * K_i_P2_cat2 * K_i_P3_cat3 * K_i_P4_cat4 /
        (Keq * K_i_S1_cat1 * K_i_S2_cat2 * K_i_S3_cat3 * K_i_S4_cat4),
        0.0,
    )
    Rate =
        (
            (
                Vmax_a *
                (S1_cat1 / K_a_S1_cat1) *
                (S2_cat2 / K_a_S2_cat2) *
                (S3_cat3 / K_a_S3_cat3) *
                (S4_cat4 / K_a_S4_cat4) -
                Vmax_a_rev *
                (P1_cat1 / K_a_P1_cat1) *
                (P2_cat2 / K_a_P2_cat2) *
                (P3_cat3 / K_a_P3_cat3) *
                (P4_cat4 / K_a_P4_cat4)
            ) *
            (Z_a_cat^(n - 1)) *
            (Z_a_reg^n) +
            L *
            (
                Vmax_i *
                (S1_cat1 / K_i_S1_cat1) *
                (S2_cat2 / K_i_S2_cat2) *
                (S3_cat3 / K_i_S3_cat3) *
                (S4_cat4 / K_i_S4_cat4) -
                Vmax_i_rev *
                (P1_cat1 / K_i_P1_cat1) *
                (P2_cat2 / K_i_P2_cat2) *
                (P3_cat3 / K_i_P3_cat3) *
                (P4_cat4 / K_i_P4_cat4)
            ) *
            (Z_i_cat^(n - 1)) *
            (Z_i_reg^n)
        ) / ((Z_a_cat^n) * (Z_a_reg^n) + L * (Z_i_cat^n) * (Z_i_reg^n))

    return Rate
end

@inline function calculate_z_reg(
    R1_reg1,
    R2_reg1,
    R3_reg1,
    R1_reg2,
    R2_reg2,
    R3_reg2,
    R1_reg3,
    R2_reg3,
    R3_reg3,
    K_R1_reg1,
    K_R2_reg1,
    K_R3_reg1,
    K_R1_reg2,
    K_R2_reg2,
    K_R3_reg2,
    K_R1_reg3,
    K_R2_reg3,
    K_R3_reg3,
)
    Z_reg = (
        (1 + R1_reg1 / K_R1_reg1 + R2_reg1 / K_R2_reg1 + R3_reg1 / K_R3_reg1) *
        (1 + R1_reg2 / K_R1_reg2 + R2_reg2 / K_R2_reg2 + R3_reg2 / K_R3_reg2) *
        (1 + R1_reg3 / K_R1_reg3 + R2_reg3 / K_R2_reg3 + R3_reg3 / K_R3_reg3)
    )
    return Z_reg
end

@inline function calculate_z_cat(
    S1_cat1,
    S2_cat2,
    S3_cat3,
    S4_cat4,
    P1_cat1,
    P2_cat2,
    P3_cat3,
    P4_cat4,
    K_S1_cat1,
    K_S2_cat2,
    K_S3_cat3,
    K_S4_cat4,
    K_P1_cat1,
    K_P2_cat2,
    K_P3_cat3,
    K_P4_cat4,
    alpha_S1_P2,
    alpha_S1_P3,
    alpha_S1_P4,
    alpha_S2_P1,
    alpha_S2_P3,
    alpha_S2_P4,
    alpha_S3_P1,
    alpha_S3_P2,
    alpha_S3_P4,
    alpha_S4_P1,
    alpha_S4_P2,
    alpha_S4_P3,
)
    Z_cat = (
        (
            (1 + S1_cat1 / K_S1_cat1) *
            (1 + S2_cat2 / K_S2_cat2) *
            (1 + S3_cat3 / K_S3_cat3) *
            (1 + S4_cat4 / K_S4_cat4)
        ) - 1 +
        (
            (1 + P1_cat1 / K_P1_cat1) *
            (1 + P2_cat2 / K_P2_cat2) *
            (1 + P3_cat3 / K_P3_cat3) *
            (1 + P4_cat4 / K_P4_cat4)
        ) +
        (S1_cat1 / K_S1_cat1) *
        (1 + alpha_S1_P2 * (P2_cat2 / K_P2_cat2)) *
        (1 + alpha_S1_P3 * (P3_cat3 / K_P3_cat3)) *
        (1 + alpha_S1_P4 * (P4_cat4 / K_P4_cat4)) - (S1_cat1 / K_S1_cat1) +
        (S2_cat2 / K_S2_cat2) *
        (1 + alpha_S2_P1 * (P1_cat1 / K_P1_cat1)) *
        (1 + alpha_S2_P3 * (P3_cat3 / K_P3_cat3)) *
        (1 + alpha_S2_P4 * (P4_cat4 / K_P4_cat4)) - (S2_cat2 / K_S2_cat2) +
        (S3_cat3 / K_S3_cat3) *
        (1 + alpha_S3_P1 * (P1_cat1 / K_P1_cat1)) *
        (1 + alpha_S3_P2 * (P2_cat2 / K_P2_cat2)) *
        (1 + alpha_S3_P4 * (P4_cat4 / K_P4_cat4)) - (S3_cat3 / K_S3_cat3) +
        (S4_cat4 / K_S4_cat4) *
        (1 + alpha_S4_P1 * (P1_cat1 / K_P1_cat1)) *
        (1 + alpha_S4_P2 * (P2_cat2 / K_P2_cat2)) *
        (1 + alpha_S4_P3 * (P3_cat3 / K_P3_cat3)) - (S4_cat4 / K_S4_cat4) +
        (P1_cat1 / K_P1_cat1) *
        (1 + alpha_S2_P1 * (S2_cat2 / K_S2_cat2)) *
        (1 + alpha_S3_P1 * (S3_cat3 / K_S3_cat3)) *
        (1 + alpha_S4_P1 * (S4_cat4 / K_S4_cat4)) - (P1_cat1 / K_P1_cat1) +
        (P2_cat2 / K_P2_cat2) *
        (1 + alpha_S1_P2 * (S1_cat1 / K_S1_cat1)) *
        (1 + alpha_S3_P2 * (S3_cat3 / K_S3_cat3)) *
        (1 + alpha_S4_P2 * (S4_cat4 / K_S4_cat4)) - (P2_cat2 / K_P2_cat2) +
        (P3_cat3 / K_P3_cat3) *
        (1 + alpha_S1_P3 * (S1_cat1 / K_S1_cat1)) *
        (1 + alpha_S2_P3 * (S2_cat2 / K_S2_cat2)) *
        (1 + alpha_S4_P3 * (S4_cat4 / K_S4_cat4)) - (P3_cat3 / K_P3_cat3) +
        (P4_cat4 / K_P4_cat4) *
        (1 + alpha_S1_P4 * (S1_cat1 / K_S1_cat1)) *
        (1 + alpha_S2_P4 * (S2_cat2 / K_S2_cat2)) *
        (1 + alpha_S3_P4 * (S3_cat3 / K_S3_cat3)) - (P4_cat4 / K_P4_cat4) +
        (
            alpha_S1_P3 *
            alpha_S1_P4 *
            alpha_S2_P3 *
            alpha_S2_P4 *
            (S1_cat1 / K_S1_cat1) *
            (S2_cat2 / K_S2_cat2) *
            (P3_cat3 / K_P3_cat3) *
            (P4_cat4 / K_P4_cat4)
        ) +
        (
            alpha_S1_P2 *
            alpha_S1_P4 *
            alpha_S3_P2 *
            alpha_S3_P4 *
            (S1_cat1 / K_S1_cat1) *
            (P2_cat2 / K_P2_cat2) *
            (S3_cat3 / K_S3_cat3) *
            (P4_cat4 / K_P4_cat4)
        ) +
        (
            alpha_S1_P2 *
            alpha_S1_P3 *
            alpha_S4_P2 *
            alpha_S4_P3 *
            (S1_cat1 / K_S1_cat1) *
            (P2_cat2 / K_P2_cat2) *
            (P3_cat3 / K_P3_cat3) *
            (S4_cat4 / K_S4_cat4)
        ) +
        (
            alpha_S2_P1 *
            alpha_S2_P4 *
            alpha_S3_P1 *
            alpha_S3_P4 *
            (P1_cat1 / K_P1_cat1) *
            (S2_cat2 / K_S2_cat2) *
            (S3_cat3 / K_S3_cat3) *
            (P4_cat4 / K_P4_cat4)
        ) +
        (
            alpha_S2_P1 *
            alpha_S2_P3 *
            alpha_S4_P1 *
            alpha_S4_P3 *
            (P1_cat1 / K_P1_cat1) *
            (S2_cat2 / K_S2_cat2) *
            (P3_cat3 / K_P3_cat3) *
            (S4_cat4 / K_S4_cat4)
        ) +
        (
            alpha_S3_P1 *
            alpha_S3_P2 *
            alpha_S4_P1 *
            alpha_S4_P2 *
            (P1_cat1 / K_P1_cat1) *
            (P2_cat2 / K_P2_cat2) *
            (S3_cat3 / K_S3_cat3) *
            (S4_cat4 / K_S4_cat4)
        )
    )
    return Z_cat
end

"Generate the names of the parameters for the rate equation using the same input as @derive_general_mwc_rate_eq"
function generate_param_names(processed_input)
    metab_binding_sites = [site for site in keys(processed_input) if site ∈ [:cat1, :cat2, :cat3, :cat4, :cat5, :reg1, :reg2, :reg3]]
    param_names = (:L, :Vmax_a, :Vmax_i)
    for site in metab_binding_sites
        for metab in processed_input[site]
            param_names = (param_names..., Symbol("K_a_", metab, "_", site), Symbol("K_i_", metab, "_", site))
        end
    end
    catalytic_sites = [site for site in keys(processed_input) if site ∈ [:cat1, :cat2, :cat3, :cat4, :cat5]]
    for site in catalytic_sites
        for substrate in processed_input[site]
            if substrate ∈ processed_input[:substrates]
                for product in processed_input[:products]
                    if product ∉ processed_input[site]
                        param_names = (param_names..., Symbol("alpha_", substrate, "_", product))
                    end
                end
            end
        end
    end
    return param_names
end

"Generate the names of the metabolites for the rate equation using the same input as @derive_general_mwc_rate_eq"
function generate_metab_names(processed_input)
    metab_names = ()
    for site in keys(processed_input)
        if site ∈ [:cat1, :cat2, :cat3, :cat4, :cat5, :reg1, :reg2, :reg3]
            metab_names = (metab_names..., processed_input[site]...)
        end
    end
    return Tuple(unique(metab_names))
end
