# using TestEnv
# TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test, BenchmarkTools

##
#=
TODO: Add tests for the following:

Use Symbolics.jl to test that the rate equation has all possible combinations of
substrates, products, and regulators with corresponding alphas
=#

##
#tests for rate function speed and allocations
substrates = [Symbol(:S, i) for i = 1:rand(1:3)]
products = [Symbol(:P, i) for i = 1:rand(1:3)]
inhibitors = [Symbol(:I, i) for i = 1:rand(1:2)]
regulators = [Symbol(:R, i) for i = 1:rand(1:5)]
Keq = 20_000.0
enzyme_parameters = (;
    substrates = substrates,
    products = products,
    inhibitors = inhibitors,
    regulators = regulators,
    Keq = Keq,
    oligomeric_state = rand(1:4),
    rate_equation_name = :rand_enz_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}(rand(length(metab_names)))
benchmark_result = @benchmark rand_enz_rate_equation($(metabs_nt), $(params_nt), $(Keq))
@test mean(benchmark_result.times) <= 800 #ns
@test benchmark_result.allocs == 0
benchmark_result = @benchmark rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test mean(benchmark_result.times) <= 800 #ns
@test benchmark_result.allocs <= 1

#test Rate < 0 when [Substrates] = 0
metabs_nt = NamedTuple{metab_names}((
    zeros(length(substrates))...,
    rand(length(products))...,
    rand(length(inhibitors))...,
    rand(length(regulators))...,
))
@test rand_enz_rate_equation(metabs_nt, params_nt, enzyme_parameters.Keq) < 0.0

#test Rate > 0 when [Products] = 0
metabs_nt = NamedTuple{metab_names}((
    rand(length(substrates))...,
    zeros(length(products))...,
    rand(length(inhibitors))...,
    rand(length(regulators))...,
))
@test rand_enz_rate_equation(metabs_nt, params_nt, Keq) > 0.0

#test Rate = Vmax_a when [Substrates] and [Activators] -> Inf and Vmax_a_rev when [Products] and [Activators] -> Inf
Vmax_a = rand()
Vmax_i = rand()
params_vec = []
for param_name in propertynames(params_nt)
    if startswith(string(param_name), "K_a")
        push!(params_vec, 1e-3)
    elseif startswith(string(param_name), "alpha_")
        push!(params_vec, rand([0.0, 1.0]))
    elseif param_name == :Vmax_a
        push!(params_vec, Vmax_a)
    elseif param_name == :Vmax_i
        push!(params_vec, Vmax_i)
    else
        push!(params_vec, 1.0)
    end
end
params_nt = NamedTuple{param_names}(params_vec)
metabs_nt = NamedTuple{metab_names}((
    (1e3 .* ones(length(substrates)))...,
    zeros(length(products))...,
    zeros(length(inhibitors))...,
    1e3 .* ones(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax_a, rtol = 1e-2)
global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_params = [Symbol("K_a_", Symbol(:S, i)) for i = 1:length(substrates)]
product_params = [Symbol("K_a_", Symbol(:P, i)) for i = 1:length(products)]
for param in propertynames(params_nt)
    if param ∈ substrate_params
        global substrate_Ks *= params_nt[param]
    elseif param ∈ product_params
        global product_Ks *= params_nt[param]
    end
end
Vmax_a_rev = Vmax_a * product_Ks / (Keq * substrate_Ks)
metabs_nt = NamedTuple{metab_names}((
    (zeros(length(substrates)))...,
    1e3 .* ones(length(products))...,
    zeros(length(inhibitors))...,
    1e3 .* ones(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_a_rev, rtol = 1e-2)

#test Rate = Vmax_i when [Substrates] and [Allosteric Inhibitors] -> Inf and Vmax_i_rev when [Products] and [Allosteric Inhibitors] -> Inf
Vmax_a = rand()
Vmax_i = rand()
params_vec = []
for param_name in propertynames(params_nt)
    if startswith(string(param_name), "K_i")
        push!(params_vec, 1e-3)
    elseif param_name == :Vmax_a
        push!(params_vec, Vmax_a)
    elseif param_name == :Vmax_i
        push!(params_vec, Vmax_i)
    elseif startswith(string(param_name), "alpha_")
        push!(params_vec, rand([0.0, 1.0]))
    else
        push!(params_vec, 1.0)
    end
end
params_nt = NamedTuple{param_names}(params_vec)
metabs_nt = NamedTuple{metab_names}((
    (1e3 .* ones(length(substrates)))...,
    zeros(length(products))...,
    zeros(length(inhibitors))...,
    1e3 .* ones(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), Vmax_i, rtol = 1e-2)
global substrate_Ks = 1.0
global product_Ks = 1.0
substrate_params = [Symbol("K_i_", Symbol(:S, i)) for i = 1:length(substrates)]
product_params = [Symbol("K_i_", Symbol(:P, i)) for i = 1:length(products)]
for param in propertynames(params_nt)
    if param ∈ substrate_params
        global substrate_Ks *= params_nt[param]
    elseif param ∈ product_params
        global product_Ks *= params_nt[param]
    end
end
Vmax_i_rev = Vmax_i * product_Ks / (Keq * substrate_Ks)
metabs_nt = NamedTuple{metab_names}((
    (zeros(length(substrates)))...,
    1e3 .* ones(length(products))...,
    zeros(length(inhibitors))...,
    1e3 .* ones(length(regulators))...,
))
@test isapprox(rand_enz_rate_equation(metabs_nt, params_nt, Keq), -Vmax_i_rev, rtol = 1e-2)

#test Rate = 0 when [Inhibitors] -> Inf
params_nt = NamedTuple{param_names}(rand(length(param_names)))
metabs_nt = NamedTuple{metab_names}((
    rand(length(substrates))...,
    rand(length(products))...,
    1e3 .* ones(length(inhibitors))...,
    rand(length(regulators))...,
))
rand_enz_rate_equation(metabs_nt, params_nt, Keq)
@test isapprox(1.0 - rand_enz_rate_equation(metabs_nt, params_nt, Keq), 1.0, atol = 1e-2)

#test Rate is unchanged regardless of regulators levels when Vmax_a = Vmax_i and Ka = Ki for all reactants
Vmax_val = rand()
params_vec = []
# First pass: set all parameters with random values or Vmax_val
for param_name in propertynames(params_nt)
    if param_name == :Vmax_a || param_name == :Vmax_i
        push!(params_vec, Vmax_val)
    else
        push!(params_vec, rand())
    end
end
# Second pass: set Ki = Ka for substrates, products, and inhibitors only
for (i, param_name) in enumerate(param_names)
    if startswith(string(param_name), "K_i_")
        metab_suffix = string(param_name)[5:end]
        if metab_suffix[1] in ['S', 'P', 'I']  # not regulators
            corresponding_Ka = Symbol("K_a_", metab_suffix)
            Ka_index = findfirst(x -> x == corresponding_Ka, param_names)
            if Ka_index !== nothing
                params_vec[i] = params_vec[Ka_index]
            end
        end
    end
end
params_nt = NamedTuple{param_names}(params_vec)

# Test with one set of regulator concentrations
substr_conc = rand(length(substrates))
prod_conc = rand(length(products))
inhib_conc = rand(length(inhibitors))
metabs_nt_first_reg = NamedTuple{metab_names}((
    substr_conc...,
    prod_conc...,
    inhib_conc...,
    rand(length(regulators))...,
))
rate_first_reg = rand_enz_rate_equation(metabs_nt_first_reg, params_nt, Keq)

# Test with second set of regulator concentrations
metabs_nt_second_reg = NamedTuple{metab_names}((
    substr_conc...,
    prod_conc...,
    inhib_conc...,
    rand(length(regulators))...,
))
rate_second_reg = rand_enz_rate_equation(metabs_nt_second_reg, params_nt, Keq)

@test isapprox(rate_low_reg, rate_high_reg, rtol = 1e-6)
