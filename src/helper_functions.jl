# WORK IN PROGRESS

# using Symbolics

function display_rate_equation(
    rate_equation::Function,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
)
    vect = []
    for x in metab_names
        push!(vect, (@variables $(x))...)
    end
    nt_metab_names = NamedTuple{metab_names}(vect)
    vect = []
    for x in param_names
        push!(vect, (@variables $(x))...)
    end
    nt_param_names = NamedTuple{param_names}(vect)
    @variables Keq
    sym_rate_eqn = rate_equation(nt_metab_names, nt_param_names, Keq)

    # Display the rate equation.
    display(rate_eqn_expr)
end
