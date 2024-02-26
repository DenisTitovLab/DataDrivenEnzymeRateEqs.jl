module EnzymeFitting
include("general_rate_equation_derivation.jl")
include("rate_equation_fitting.jl")
include("optimal_rate_equation_selection.jl")

export @derive_general_mwc_rate_eq
export fit_rate_equation
export optimal_rate_equation_selection

end
