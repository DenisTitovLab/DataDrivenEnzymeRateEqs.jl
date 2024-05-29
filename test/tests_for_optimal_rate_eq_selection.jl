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


##
#test the ability of `data_driven_rate_equation_selection` to recover the MWC rate_equation and params used to generated data for an arbitrary enzyme

data_gen_rate_equation_Keq = 1.0
data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) = (1 / params.K_a_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P)
data_gen_param_names = (:Vmax_a, :K_a_S, :K_a_P)
metab_names = (:S, :P)
params = (Vmax=10.0, K_a_S=1e-3, K_a_P=5e-3)
#create DataFrame of simulated data
num_datapoints = 10
num_figures = 4
S_concs = Float64[]
P_concs = Float64[]
sources = String[]

for i in 1:num_figures
    if i < num_figures ÷ 2
        for S in range(0, rand(1:10) * params.K_a_S, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, S)
            push!(P_concs, 0.0)
            push!(sources, "Figure$i")
        end
    else
        for P in range(0, rand(1:10) * params.K_a_P, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, 0.0)
            push!(P_concs, P)
            push!(sources, "Figure$i")
        end
    end
end
data = DataFrame(S=S_concs, P=P_concs, source=sources)
noise_sd = 0.2
data.Rate = [data_gen_rate_equation(row, params, data_gen_rate_equation_Keq) * (1 + noise_sd * randn()) for row in eachrow(data)]
data

enzyme_parameters = (; substrates=[:S,], products=[:P], regulators=[], Keq=1.0, oligomeric_state=1, rate_equation_name=:mwc_derived_rate_equation)

metab_names, derived_param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
mwc_derived_rate_equation_no_Keq(nt_metabs, nt_params) = mwc_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(mwc_derived_rate_equation_no_Keq, data, metab_names, derived_param_names, (3, 7), true)

#Display best equation with 3 parameters. Compare with data_gen_rate_equation with Vmax=1
#TODO: remove the filtering for 3 parameters after we add the automatic determination of the best number of parameters
nt_param_removal_code = filter(x -> x.num_params .== 3, selection_result.test_results).nt_param_removal_codes[1]

using Symbolics
selected_sym_rate_equation = display_rate_equation(mwc_derived_rate_equation, metab_names, derived_param_names; nt_param_removal_code=nt_param_removal_code)
original_sym_rate_equation = display_rate_equation(data_gen_rate_equation, metab_names, data_gen_param_names)
println("Selected MWC rate equation:")
println(simplify(selected_sym_rate_equation))
println("Original MWC rate equation:")
println(simplify(original_sym_rate_equation))
@test simplify(original_sym_rate_equation - selected_sym_rate_equation) == 0

##
#test the ability of `data_driven_rate_equation_selection` to recover the QSSA rate_equation and params used to generated data for an arbitrary enzyme

data_gen_rate_equation_Keq = 1.0
data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) = (1 / params.K_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (1 + metabs.S / params.K_S + metabs.P / params.K_P)
data_gen_param_names = (:Vmax, :K_S, :K_P)
metab_names = (:S, :P)
params = (Vmax=10.0, K_S=1e-3, K_P=5e-3)
#create DataFrame of simulated data
num_datapoints = 10
num_figures = 4
S_concs = Float64[]
P_concs = Float64[]
sources = String[]

for i in 1:num_figures
    if i < num_figures ÷ 2
        for S in range(0, rand(1:10) * params.K_S, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, S)
            push!(P_concs, 0.0)
            push!(sources, "Figure$i")
        end
    else
        for P in range(0, rand(1:10) * params.K_P, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, 0.0)
            push!(P_concs, P)
            push!(sources, "Figure$i")
        end
    end
end
data = DataFrame(S=S_concs, P=P_concs, source=sources)
noise_sd = 0.2
data.Rate = [data_gen_rate_equation(row, params, data_gen_rate_equation_Keq) * (1 + noise_sd * randn()) for row in eachrow(data)]
data

enzyme_parameters = (; substrates=[:S,], products=[:P], regulators=[], Keq=1.0, rate_equation_name=:qssa_derived_rate_equation)

metab_names, derived_param_names = @derive_general_qssa_rate_eq(enzyme_parameters)
qssa_derived_rate_equation_no_Keq(nt_metabs, nt_params) = qssa_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(qssa_derived_rate_equation_no_Keq, data, metab_names, derived_param_names, (1, 4), true)

#Display best equation with 3 parameters. Compare with data_gen_rate_equation with Vmax=1
#TODO: remove the filtering for 3 parameters after we add the automatic determination of the best number of parameters
nt_param_removal_code = filter(x -> x.num_params .== 3, selection_result.test_results).nt_param_removal_codes[1]

using Symbolics
selected_sym_rate_equation = display_rate_equation(qssa_derived_rate_equation, metab_names, derived_param_names; nt_param_removal_code=nt_param_removal_code)
original_sym_rate_equation = display_rate_equation(data_gen_rate_equation, metab_names, data_gen_param_names)
println("Selected QSSA rate equation:")
println(simplify(selected_sym_rate_equation))
println("Original QSSA rate equation:")
println(simplify(original_sym_rate_equation))
@test simplify(original_sym_rate_equation - selected_sym_rate_equation) == 0
