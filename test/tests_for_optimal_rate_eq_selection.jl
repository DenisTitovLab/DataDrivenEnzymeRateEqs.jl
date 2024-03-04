# using TestEnv
# TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools

#test forward_selection_next_param_removal_codes
num_metabolites = rand(4:8)
n_alphas = rand(1:4)
num_previous_step_params = rand((2+num_metabolites):(3+2*num_metabolites))
num_params = num_previous_step_params - 1
param_names = (
    :L,
    :Vmax_a,
    :Vmax_i,
    [Symbol(:K_a, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:K_i, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:alpha, "_$(i)") for i = 1:n_alphas]...,
)
param_removal_code_names = (
    [
        Symbol(replace(string(param_name), "_a" => "")) for
        param_name in param_names if !contains(string(param_name), "_i")
    ]...,
)
all_param_removal_codes = DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes(param_names)
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]
nt_funct_output_param_subset_codes = DataDrivenEnzymeRateEqs.forward_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names
)
funct_output_param_subset_codes = [values(nt) for nt in nt_funct_output_param_subset_codes]
#ensure that funct_output_param_subset_codes have one less parameter than previous_param_removal_codes
@test all(
    length(param_names) - n_alphas -
    sum(funct_output_param_subset_code[1:(end-n_alphas)] .> 0) ==
    (num_previous_step_params - 1) for
    funct_output_param_subset_code in funct_output_param_subset_codes
)
#ensure that non-zero elements from previous_param_removal_codes are present in > 1 of the funct_output_param_subset_code but less than the max_matches
count_matches = []
non_zero_code_combos_per_param = ()
for param_name in param_names
    if param_name == :L
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 1)
    elseif occursin("Vmax_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 2)
    elseif occursin("K_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 3)
    elseif occursin("alpha", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 1)
    end
end
for funct_output_param_subset_code in funct_output_param_subset_codes
    count = 0
    max_matches_vect = Int[]
    for previous_param_removal_code in previous_param_removal_codes
        count +=
            funct_output_param_subset_code[1:end-n_alphas] .* previous_param_removal_code[1:end-n_alphas] ==
            previous_param_removal_code[1:end-n_alphas] .^ 2
        push!(max_matches_vect, sum((previous_param_removal_code[1:end-n_alphas] .== 0) .* non_zero_code_combos_per_param[1:end-n_alphas]))
    end
    max_matches = maximum(max_matches_vect)
    push!(count_matches, max_matches >= count > 0)
end
@test all(count_matches)

#test reverse_selection_next_param_removal_codes
num_metabolites = rand(4:8)
n_alphas = rand(1:4)
num_previous_step_params = rand((2+num_metabolites):(3+2*num_metabolites))
num_params = num_previous_step_params + 1
param_names = (
    :L,
    :Vmax_a,
    :Vmax_i,
    [Symbol(:K_a, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:K_i, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:alpha, "_$(i)") for i = 1:n_alphas]...,
)
param_removal_code_names = (
    [
        Symbol(replace(string(param_name), "_a" => "")) for
        param_name in param_names if !contains(string(param_name), "_i")
    ]...,
)
all_param_removal_codes = DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes(param_names)
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]

nt_funct_output_param_subset_codes = DataDrivenEnzymeRateEqs.reverse_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
)
funct_output_param_subset_codes = [values(nt) for nt in nt_funct_output_param_subset_codes]
#ensure that funct_output_param_subset_codes have one more parameter than previous_param_removal_codes
@test all(
    length(param_names) - n_alphas -
    sum(funct_output_param_subset_code[1:(end-n_alphas)] .> 0) ==
    (num_previous_step_params + 1) for
    funct_output_param_subset_code in funct_output_param_subset_codes
)
#ensure that non-zero elements from each funct_output_param_subset_codes are present in > 1 of the previous_param_removal_codes but less than the max_matches
count_matches = []
non_zero_code_combos_per_param = ()
for param_name in param_names
    if param_name == :L
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 1)
    elseif occursin("Vmax_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 2)
    elseif occursin("K_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 3)
    elseif occursin("alpha", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 1)
    end
end
for funct_output_param_subset_code in funct_output_param_subset_codes
    count = 0
    for previous_param_removal_code in previous_param_removal_codes
        count +=
            funct_output_param_subset_code[1:end-n_alphas] == previous_param_removal_code[1:end-n_alphas] .* (funct_output_param_subset_code[1:end-n_alphas] .!= 0)
    end
    max_matches = sum((funct_output_param_subset_code[1:end-n_alphas] .== 0) .* non_zero_code_combos_per_param[1:end-n_alphas])
    push!(count_matches, max_matches >= count > 0)
end
@test all(count_matches)


# #
# #test the ability of `data_driven_rate_equation_selection` to recover the rate_equation and params used to generated data for an arbitrary enzyme
data_gen_rate_equation(metabs, params) = params.Vmax * (metabs.S / params.K_S - metabs.P / params.K_P) / (1 + metabs.S / params.K_S + metabs.P / params.K_P)
param_names = (:Vmax, :K_S, :K_P)
metab_names = (:S, :P)
params = (Vmax=10.0, K_S=1.0, K_P=5.0)
data = DataFrame(S=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
data.P = [0.0 for row in eachrow(data)]
noise_sd = 0.2
data.Rate = [data_gen_rate_equation(row, params) * (1 + noise_sd * randn()) for row in eachrow(data)]
data.source = ["Figure1" for i in 1:nrow(data)]
fit_result = fit_rate_equation(data_gen_rate_equation, data, metab_names, param_names; n_iter=20)
@test isapprox(fit_result.params.K_S, params.K_S, rtol=3 * noise_sd)

enzyme_parameters = (; substrates=[:S,], products=[:P], cat1=[:S, :P], reg1=[], reg2=[], Keq=1.0, oligomeric_state=1, rate_equation_name=:derived_rate_equation)
metab_names, param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
nt_params = NamedTuple{param_names}(rand(length(param_names)))
nt_metabs = NamedTuple{metab_names}(rand(length(metab_names)))
derived_rate_equation(nt_metabs, nt_params) = derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
fit_result = fit_rate_equation(derived_rate_equation, data, metab_names, param_names; n_iter=20)
selection_result = data_driven_rate_equation_selection(derived_rate_equation, data, metab_names, param_names, (3, 7), true)
