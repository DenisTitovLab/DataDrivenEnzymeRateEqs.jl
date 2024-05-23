using Symbolics

"""
    display_rate_equation(
    rate_equation::Function,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    nt_param_removal_code = nothing
)

Return the symbolic rate equation for the given `rate_equation` function.

# Arguments
- `rate_equation::Function`: The rate equation function.
- `metab_names::Tuple{Symbol,Vararg{Symbol}}`: The names of the metabolites.
- `param_names::Tuple{Symbol,Vararg{Symbol}}`: The names of the parameters.
- `nt_param_removal_code::NamedTuple`: The named tuple of the parameters to remove from the rate equation.
"""
function display_rate_equation(
    rate_equation::Function,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    nt_param_removal_code = nothing
)

    metab_name_syms = []
    for x in metab_names
        push!(metab_name_syms, (@variables $(x))...)
    end
    nt_metab_names = NamedTuple{metab_names}(metab_name_syms)

    param_name_syms = []
    for x in param_names
        push!(param_name_syms, (@variables $(x))...)
    end
    if !isnothing(nt_param_removal_code)
        param_name_syms =
            param_subset_select(param_name_syms, param_names, nt_param_removal_code)
    end
    nt_param_names = NamedTuple{param_names}(param_name_syms)

    @variables Keq
    sym_rate_eqn = rate_equation(nt_metab_names, nt_param_names, Keq)
    return sym_rate_eqn
end
