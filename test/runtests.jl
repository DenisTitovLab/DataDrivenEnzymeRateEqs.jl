using EnzymeFitting
using Test
using SafeTestsets

@time begin
    @safetestset "MWC Rate Eq Derivation" begin
        include("tests_for_general_rate_eq_derivation.jl")
    end
    @safetestset "MWC Rate Eq Fitting" begin
        include("tests_for_rate_eq_fitting.jl")
    end
    @safetestset "MWC Rate Eq Subset Selection" begin
        include("tests_for_optimal_rate_eq_selection.jl")
    end
end
