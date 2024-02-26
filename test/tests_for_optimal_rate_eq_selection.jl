# using TestEnv
# TestEnv.activate()

##
using EnzymeFitting, Test, BenchmarkTools

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
all_param_removal_codes = EnzymeFitting.calculate_all_parameter_removal_codes(param_names)
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]
funct_output_param_subset_codes = EnzymeFitting.forward_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
)
#test1
@test all(
    length(param_names) - n_alphas -
    sum(funct_output_param_subset_code[1:(end-n_alphas)] .> 0) ==
    (num_previous_step_params - 1) for
    funct_output_param_subset_code in funct_output_param_subset_codes
)
#test2
count_matches = Bool[]
for funct_output_param_subset_code in funct_output_param_subset_codes
    count = 0
    for previous_param_removal_code in previous_param_removal_codes
        count +=
            funct_output_param_subset_code .* [previous_param_removal_code...] ==
            [previous_param_removal_code...] .^ 2
    end
    push!(count_matches, count > 0)
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
all_param_removal_codes = EnzymeFitting.calculate_all_parameter_removal_codes(param_names)
param_subset_codes_with_num_params = [
    x for x in all_param_removal_codes if
    length(param_names) - sum(values(x[1:(end-n_alphas)]) .> 0) - n_alphas ==
    num_previous_step_params
]
previous_param_removal_codes =
    [rand(param_subset_codes_with_num_params) for i = 1:rand(1:20)]

funct_output_param_subset_codes = EnzymeFitting.reverse_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
)

#test1
@assert all(
    length(param_names) - n_alphas -
    sum(funct_output_param_subset_code[1:(end-n_alphas)] .> 0) ==
    (num_previous_step_params + 1) for
    funct_output_param_subset_code in funct_output_param_subset_codes
)

# #test2
# count_matches = []
# for previous_param_removal_code in previous_param_removal_codes
#     count = 0
#     for funct_output_param_subset_code in funct_output_param_subset_codes
#         count +=
#             funct_output_param_subset_code[1:end-n_alphas] .* [previous_param_removal_code[1:end-n_alphas]...] ==
#             [funct_output_param_subset_code[1:end-n_alphas]...] .^ 2
#     end
#     push!(count_matches, count)
# end
# funct_output_param_subset_codes
# previous_param_removal_codes
# funct_output_param_subset_codes[2] .* [previous_param_removal_codes[1]...] ==
# [funct_output_param_subset_codes[2]...] .^ 2
# count_matches
# @assert all(count_matches)
