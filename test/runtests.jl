using DataDrivenEnzymeRateEqs
using Test
using SafeTestsets

@time begin
    @safetestset "MWC Rate Eq Derivation" begin
        include("tests_for_mwc_general_rate_eq_derivation.jl")
    end
    @safetestset "MWC Rate Eq Fitting" begin
        include("tests_for_rate_eq_fitting.jl")
    end
    @safetestset "MWC Rate Eq Subset Selection" begin
        include("tests_for_optimal_rate_eq_selection.jl")
    end
    @safetestset "Rapid Equilibrium Rate Eq Derivation" begin
        include("tests_for_rapid_equilibrium_general_rate_eq_derivation.jl")
    end
end
