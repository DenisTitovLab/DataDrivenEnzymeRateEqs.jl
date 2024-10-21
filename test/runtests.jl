using DataDrivenEnzymeRateEqs
using Test
using SafeTestsets

@time begin
    @safetestset "MWC Rate Eq Derivation" begin
        include("tests_for_mwc_general_rate_eq_derivation_macro.jl")
    end
    @safetestset "QSSA Rate Eq Derivation" begin
        include("tests_for_qssa_general_rate_eq_derivation_macro.jl")
    end
    @safetestset "MWC and QSSA Rate Eq Fitting" begin
        include("tests_mwc_qssa_for_rate_eq_fitting.jl")
    end
    @safetestset "Plotting" begin
        include("tests_for_plotting.jl")
    end
    @safetestset "MWC Rate Eq Subset Selection" begin
        include("tests_for_optimal_rate_eq_selection.jl")
    end
end
