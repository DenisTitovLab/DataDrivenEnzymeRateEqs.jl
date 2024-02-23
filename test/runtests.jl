using EnzymeFitting
using Test
using SafeTestsets

@time begin
    @time @safetestset "MWC Rate Eq Derivation" begin
        include("tests_for_rate_eq_derivation.jl")
    end
end
