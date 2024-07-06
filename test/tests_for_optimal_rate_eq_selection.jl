using TestEnv
TestEnv.activate()

##
using DataDrivenEnzymeRateEqs, Test
using CMAEvolutionStrategy, DataFrames, CSV, Statistics
using BenchmarkTools

#test forward_selection_next_param_removal_codes
num_metabolites = rand(4:8)
metab_names = Tuple(Symbol("Metabolite$(i)") for i = 1:num_metabolites)
n_alphas = rand(1:4)
max_zero_alpha = 1 + ceil(Int, length(metab_names) / 2)
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
        Symbol(replace(string(param_name), "_a_" => "_allo_", "Vmax_a" => "Vmax_allo"))
        for param_name in param_names if
        !contains(string(param_name), "_i") && param_name != :Vmax && param_name != :L
    ]...,
)
practically_unidentifiable_params = ()
all_param_removal_codes =
    DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes(param_names, ())
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]
nt_previous_param_removal_codes =
    [NamedTuple{param_removal_code_names}(x) for x in previous_param_removal_codes]
param_removal_code_names

nt_funct_output_param_subset_codes =
    DataDrivenEnzymeRateEqs.forward_selection_next_param_removal_codes(
        nt_previous_param_removal_codes,
        metab_names,
        practically_unidentifiable_params,
        n_alphas,
        max_zero_alpha,
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
    if occursin("Vmax_a", string(param_name))
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
            funct_output_param_subset_code[1:(end-n_alphas)] .*
            previous_param_removal_code[1:(end-n_alphas)] ==
            previous_param_removal_code[1:(end-n_alphas)] .^ 2
        push!(
            max_matches_vect,
            sum(
                (previous_param_removal_code[1:(end-n_alphas)] .== 0) .*
                non_zero_code_combos_per_param[1:(end-n_alphas)],
            ),
        )
    end
    max_matches = maximum(max_matches_vect)
    push!(count_matches, max_matches >= count > 0)
end
@test all(count_matches)

#test reverse_selection_next_param_removal_codes
num_metabolites = rand(4:8)
metab_names = Tuple(Symbol("Metabolite$(i)") for i = 1:num_metabolites)
n_alphas = rand(1:4)
max_zero_alpha = 1 + ceil(Int, length(metab_names) / 2)
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
        Symbol(replace(string(param_name), "_a_" => "_allo_")) for
        param_name in param_names if
        !contains(string(param_name), "_i") && param_name != :Vmax && param_name != :L
    ]...,
)
practically_unidentifiable_params = ()
all_param_removal_codes =
    DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes(param_names, ())
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]
nt_previous_param_removal_codes =
    [NamedTuple{param_removal_code_names}(x) for x in previous_param_removal_codes]
nt_funct_output_param_subset_codes =
    DataDrivenEnzymeRateEqs.reverse_selection_next_param_removal_codes(
        nt_previous_param_removal_codes,
        metab_names,
        practically_unidentifiable_params,
        n_alphas,
        max_zero_alpha,
    )
funct_output_param_subset_codes = [values(nt) for nt in nt_funct_output_param_subset_codes]
#ensure that funct_output_param_subset_codes have one more parameter than nt_previous_param_removal_codes
@test all(
    length(param_names) - n_alphas -
    sum(funct_output_param_subset_code[1:(end-n_alphas)] .> 0) ==
    (num_previous_step_params + 1) for
    funct_output_param_subset_code in funct_output_param_subset_codes
)
#ensure that non-zero elements from each funct_output_param_subset_codes are present in > 1 of the nt_previous_param_removal_codes but less than the max_matches
count_matches = []
non_zero_code_combos_per_param = ()
for param_name in param_names
    if occursin("Vmax_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 2)
    elseif occursin("K_a", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 3)
    elseif occursin("alpha", string(param_name))
        global non_zero_code_combos_per_param = (non_zero_code_combos_per_param..., 1)
    end
end
for funct_output_param_subset_code in funct_output_param_subset_codes
    count = 0
    for previous_param_removal_code in values.(nt_previous_param_removal_codes)
        count +=
            funct_output_param_subset_code[1:end-n_alphas] ==
            previous_param_removal_code[1:end-n_alphas] .*
            (funct_output_param_subset_code[1:end-n_alphas] .!= 0)
    end
    max_matches = sum(
        (funct_output_param_subset_code[1:end-n_alphas] .== 0) .*
        non_zero_code_combos_per_param[1:end-n_alphas],
    )
    push!(count_matches, max_matches >= count > 0)
end
@test all(count_matches)

#test calculate_all_parameter_removal_codes_w_num_params
num_metabolites = rand(4:8)
metab_names = Tuple(Symbol("Metabolite$(i)") for i = 1:num_metabolites)
n_alphas = rand(1:4)
max_zero_alpha = 1 + ceil(Int, length(metab_names) / 2)
param_names = (
    :L,
    :Vmax_a,
    :Vmax_i,
    [Symbol(:K_a, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:K_i, "_Metabolite$(i)") for i = 1:num_metabolites]...,
    [Symbol(:alpha, "_$(i)") for i = 1:n_alphas]...,
)
practically_unidentifiable_params = ()
param_removal_code_names = (
    [
        Symbol(replace(string(param_name), "_a_" => "_allo_")) for
        param_name in param_names if
        !contains(string(param_name), "_i") && param_name != :Vmax && param_name != :L
    ]...,
)
all_param_removal_codes =
    DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes(param_names, ())

num_params = rand(
    (length(param_names)-length(param_removal_code_names)):length(param_names)-n_alphas,
)
nt_param_subset_codes_w_num_params =
    DataDrivenEnzymeRateEqs.calculate_all_parameter_removal_codes_w_num_params(
        num_params,
        all_param_removal_codes,
        param_names,
        param_removal_code_names,
        metab_names,
        practically_unidentifiable_params,
        n_alphas,
        max_zero_alpha,
    )
#ensure that funct_output_param_subset_codes have the correct number of parameters
@test all(
    length(param_names) - n_alphas - sum(values(subset_code)[1:(end-n_alphas)] .> 0) ==
    num_params for subset_code in nt_param_subset_codes_w_num_params
)

# test find_practically_unidentifiable_params on PKM2 and LDH data
#Load and process data
PKM2_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/PKM2_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
PKM2_data_for_fit.source = PKM2_data_for_fit.Article .* "_" .* PKM2_data_for_fit.Fig
data = PKM2_data_for_fit
PKM2_enzyme = (;
    substrates = [:PEP, :ADP],
    products = [:Pyruvate, :ATP],
    regulators = [:F16BP, :Phenylalanine],
    Keq = 20_000.0,
    oligomeric_state = 4,
    rate_equation_name = :pkm2_rate_equation,
)
metab_names, param_names = @derive_general_mwc_rate_eq(PKM2_enzyme)
unidentifiable_params =
    DataDrivenEnzymeRateEqs.find_practically_unidentifiable_params(data, param_names)
@test all(
    unidentifiable_param ∈ (
        :alpha_PEP_Pyruvate,
        :alpha_ADP_Pyruvate,
        :alpha_Pyruvate_F16BP,
        :alpha_Pyruvate_Phenylalanine,
        :alpha_ATP_Phenylalanine,
    ) for unidentifiable_param in unidentifiable_params
)

#test filter_param_removal_codes_to_prevent_wrong_param_combos
nt_param_removal_codes = [
    NamedTuple{(:Vmax, :K_S1, :K_S2, :K_S3, :K_S1_S2, :K_S1_S3, :K_S2_S3, :K_S1_S2_S3)}(
        combo,
    ) for combo in Iterators.product(
        [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]...,
    )
]
metab_names = (:S1, :S2, :S3)
filtered_nt =
    DataDrivenEnzymeRateEqs.filter_param_removal_codes_to_prevent_wrong_param_combos(
        nt_param_removal_codes,
        metab_names,
    )
for filtered_code in filtered_nt
    for param in [:K_S1, :K_S2, :K_S3]
        if filtered_code[param] == 1
            for param_combo in [:K_S1_S2, :K_S1_S3, :K_S2_S3, :K_S1_S2_S3]
                if occursin(string(param), string(param_combo))
                    @test filtered_code[param_combo] != 2
                end
            end
        end
    end
end

#test filter_param_removal_codes_for_max_zero_alpha
#TODO: add a test with practically_unidentifiable_alphas
num_alpha_params = rand(5:10)
num_non_alpha_params = rand(5:10)
max_zero_alpha = rand(1:3)
param_names = (
    [Symbol("param$(i)") for i = 1:num_non_alpha_params]...,
    [Symbol("alpha$(i)") for i = 1:num_alpha_params]...,
)
param_names[end-num_alpha_params+1]
nt_param_removal_codes = [
    NamedTuple{(param_names)}(combo) for
    combo in Iterators.product([[0, 1] for _ in param_names]...)
]
filtered_nt = DataDrivenEnzymeRateEqs.filter_param_removal_codes_for_max_zero_alpha(
    nt_param_removal_codes,
    practically_unidentifiable_params,
    max_zero_alpha,
)
sum_alpha = [sum(values(nt)[end-num_alpha_params+1:end]) for nt in filtered_nt]
@test all(num_alpha_params .- sum_alpha .<= max_zero_alpha)
@test any(sum_alpha .== num_alpha_params - max_zero_alpha)
@test all(sum_alpha .!= num_alpha_params - max_zero_alpha - 1)
@test any(sum_alpha .== num_alpha_params)
@test all(sum_alpha .!= num_alpha_params + 1)

#Load and process data
LDH_data_for_fit = CSV.read(joinpath(@__DIR__, "Data_for_tests/LDH_data.csv"), DataFrame)
#Add source column that uniquely identifies a figure from publication
LDH_data_for_fit.source = LDH_data_for_fit.Article .* "_" .* LDH_data_for_fit.Fig
data = LDH_data_for_fit
LDH_enzyme = (;
    substrates = [:Pyruvate, :NADH],
    products = [:Lactate, :NAD],
    regulators = [],
    Keq = 20_000.0,
    rate_equation_name = :ldh_rate_equation,
)
metab_names, param_names = @derive_general_qssa_rate_eq(LDH_enzyme)
unidentifiable_params =
    DataDrivenEnzymeRateEqs.find_practically_unidentifiable_params(data, param_names)
@test all(
    unidentifiable_param ∈ (:K_Pyruvate_NADH_Lactate_NAD,) for
    unidentifiable_param in unidentifiable_params
)


##
#test the ability of `data_driven_rate_equation_selection` to recover the MWC rate_equation and params used to generated data for an arbitrary enzyme
data_gen_rate_equation_Keq = 1.0
mwc_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) =
    (1 / params.K_a_S) * (
        metabs.S * (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P) -
        metabs.P * (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P) /
        data_gen_rate_equation_Keq
    ) / (
        (1 + metabs.S / params.K_a_S + metabs.P / params.K_a_P)^2 +
        params.L * (1 + metabs.P / params.K_a_P)^2
    )
mwc_alternative_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) =
    (1 / params.K_a_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (
        1 +
        metabs.S / params.K_a_S +
        metabs.P / params.K_a_P +
        metabs.S * metabs.P / (params.K_a_S * params.K_a_P + params.L)
    )

data_gen_param_names = (:Vmax_a, :L, :K_a_S, :K_a_P)
metab_names = (:S, :P)
params = (Vmax = 10.0, L = 10000, K_a_S = 1e-3, K_a_P = 5e-3)
#create DataFrame of simulated data
num_datapoints = 10
num_figures = 4
S_concs = Float64[]
P_concs = Float64[]
sources = String[]

for i = 1:num_figures
    if i < num_figures ÷ 2
        for S in
            range(0, rand(1:10) * params.K_a_S, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, S)
            push!(P_concs, 0.0)
            push!(sources, "Figure$i")
        end
    else
        for P in
            range(0, rand(1:10) * params.K_a_P, rand(num_datapoints÷2:num_datapoints*2))
            push!(S_concs, 0.0)
            push!(P_concs, P)
            push!(sources, "Figure$i")
        end
    end
end
data = DataFrame(S = S_concs, P = P_concs, source = sources)
noise_sd = 0.2
data.Rate = [
    mwc_data_gen_rate_equation(row, params, data_gen_rate_equation_Keq) *
    (1 + noise_sd * randn()) for row in eachrow(data)
]
data

enzyme_parameters = (;
    substrates = [:S],
    products = [:P],
    regulators = [],
    Keq = 1.0,
    oligomeric_state = 2,
    rate_equation_name = :mwc_derived_rate_equation,
)

metab_names, derived_param_names = @derive_general_mwc_rate_eq(enzyme_parameters)
mwc_derived_rate_equation_no_Keq(nt_metabs, nt_params) =
    mwc_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(
    mwc_derived_rate_equation_no_Keq,
    data,
    metab_names,
    derived_param_names,
)
reverse_selection_result = @time data_driven_rate_equation_selection(
    mwc_derived_rate_equation_no_Keq,
    data,
    metab_names,
    derived_param_names;
    forward_model_selection = false,
)

#Display best equation with 4 parameters. Compare with data_gen_rate_equation with Vmax=1
#TODO: remove the filtering for 4 parameters after we add the automatic determination of the best number of parameters
nt_param_removal_code =
    filter(x -> x.num_params .== 4, selection_result.test_results).nt_param_removal_codes[1]
nt_reverse_param_removal_code = filter(
    x -> x.num_params .== 4,
    reverse_selection_result.test_results,
).nt_param_removal_codes[1]

using Symbolics
selected_sym_rate_equation = display_rate_equation(
    mwc_derived_rate_equation,
    metab_names,
    derived_param_names;
    nt_param_removal_code = nt_param_removal_code,
)
reverse_selected_sym_rate_equation = display_rate_equation(
    mwc_derived_rate_equation,
    metab_names,
    derived_param_names;
    nt_param_removal_code = nt_reverse_param_removal_code,
)
original_sym_rate_equation =
    display_rate_equation(mwc_data_gen_rate_equation, metab_names, data_gen_param_names)
alrenative_original_sym_rate_equation = display_rate_equation(
    mwc_alternative_data_gen_rate_equation,
    metab_names,
    data_gen_param_names,
)

println("Selected MWC rate equation:")
println(selected_sym_rate_equation)
println("Reverse Selected MWC rate equation:")
println(reverse_selected_sym_rate_equation)
println("Original MWC rate equation:")
println(original_sym_rate_equation)
forward_is_reverse =
    simplify(selected_sym_rate_equation - reverse_selected_sym_rate_equation) == 0
# @test forward_is_reverse
#equation with S*P term and without it is equally likely to be selected as there's no data with S and P present. Hence the OR condition below
selected_is_original =
    simplify(original_sym_rate_equation - selected_sym_rate_equation) == 0
selected_is_original = selected_is_original isa Bool ? selected_is_original : false
selected_is_alternative =
    simplify(alrenative_original_sym_rate_equation - selected_sym_rate_equation) == 0
selected_is_alternative = selected_is_alternative isa Bool ? selected_is_alternative : false
# @test selected_is_original || selected_is_alternative

##
#test the ability of `data_driven_rate_equation_selection` to recover the QSSA rate_equation and params used to generated data for an arbitrary enzyme

data_gen_rate_equation_Keq = 1.0
qssa_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) =
    (1 / params.K_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) /
    (1 + metabs.S / params.K_S + metabs.P / params.K_P)
qssa_alternative_data_gen_rate_equation(metabs, params, data_gen_rate_equation_Keq) =
    (1 / params.K_S) * (metabs.S - metabs.P / data_gen_rate_equation_Keq) / (
        1 +
        metabs.S / params.K_S +
        metabs.P / params.K_P +
        metabs.S * metabs.P / (params.K_S * params.K_P)
    )

data_gen_param_names = (:Vmax, :K_S, :K_P)
metab_names = (:S, :P)
params = (Vmax = 10.0, K_S = 1e-3, K_P = 5e-3)
#create DataFrame of simulated data
num_datapoints = 10
num_figures = 4
S_concs = Float64[]
P_concs = Float64[]
sources = String[]

for i = 1:num_figures
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
data = DataFrame(S = S_concs, P = P_concs, source = sources)
noise_sd = 0.2
data.Rate = [
    qssa_data_gen_rate_equation(row, params, data_gen_rate_equation_Keq) *
    (1 + noise_sd * randn()) for row in eachrow(data)
]
data

enzyme_parameters = (;
    substrates = [:S],
    products = [:P],
    regulators = [],
    Keq = 1.0,
    rate_equation_name = :qssa_derived_rate_equation,
)

metab_names, derived_param_names = @derive_general_qssa_rate_eq(enzyme_parameters)
qssa_derived_rate_equation_no_Keq(nt_metabs, nt_params) =
    qssa_derived_rate_equation(nt_metabs, nt_params, enzyme_parameters.Keq)
selection_result = @time data_driven_rate_equation_selection(
    qssa_derived_rate_equation_no_Keq,
    data,
    metab_names,
    derived_param_names,
)
reverse_selection_result = @time data_driven_rate_equation_selection(
    qssa_derived_rate_equation_no_Keq,
    data,
    metab_names,
    derived_param_names;
    forward_model_selection = false,
)

#Display best equation with 3 parameters. Compare with data_gen_rate_equation with Vmax=1
#TODO: remove the filtering for 3 parameters after we add the automatic determination of the best number of parameters
nt_param_removal_code =
    filter(x -> x.num_params .== 3, selection_result.test_results).nt_param_removal_codes[1]
nt_reverse_param_removal_code = filter(
    x -> x.num_params .== 3,
    reverse_selection_result.test_results,
).nt_param_removal_codes[1]

using Symbolics
selected_sym_rate_equation = display_rate_equation(
    qssa_derived_rate_equation,
    metab_names,
    derived_param_names;
    nt_param_removal_code = nt_param_removal_code,
)
reverse_selected_sym_rate_equation = display_rate_equation(
    qssa_derived_rate_equation,
    metab_names,
    derived_param_names;
    nt_param_removal_code = nt_param_removal_code,
)
original_sym_rate_equation =
    display_rate_equation(qssa_data_gen_rate_equation, metab_names, data_gen_param_names)

println("Selected QSSA rate equation:")
println(selected_sym_rate_equation)
println("Reverse Selected QSSA rate equation:")
println(reverse_selected_sym_rate_equation)
println("Original QSSA rate equation:")
println(original_sym_rate_equation)
#equation with S*P term and without it is equally likely to be selected as there's no data with S and P present. Hence the OR condition below
forward_is_reverse =
    simplify(selected_sym_rate_equation - reverse_selected_sym_rate_equation) == 0
@test forward_is_reverse
selected_is_original =
    simplify(original_sym_rate_equation - selected_sym_rate_equation) == 0
selected_is_original = selected_is_original isa Bool ? selected_is_original : false
@test selected_is_original
