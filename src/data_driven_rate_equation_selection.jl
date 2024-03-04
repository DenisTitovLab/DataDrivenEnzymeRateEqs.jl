using Dates, CSV, DataFrames

"""
    data_driven_rate_equation_selection(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}},
        range_number_params::Tuple{Int,Int},
        forward_model_selection::Bool,
    )

This function is used to perform data-driven rate equation selection using a general rate equation and data. The function will select the best rate equation by iteratively removing parameters from the general rate equation and finding an equation that yield best test scores on data not used for fitting.

# Arguments
- `general_rate_equation::Function`: Function that takes a NamedTuple of metabolite concentrations (with `metab_names` keys) and parameters (with `param_names` keys) and returns an enzyme rate.
- `data::DataFrame`: DataFrame containing the data with column `Rate` and columns for each `metab_names` where each row is one measurement. It also needs to have a column `source` that contains a string that identifies the source of the data. This is used to calculate the weights for each figure in the publication.
- `metab_names::Tuple`: Tuple of metabolite names that correspond to the metabolites of `rate_equation` and column names in `data`.
- `param_names::Tuple`: Tuple of parameter names that correspond to the parameters of `rate_equation`.
- `range_number_params::Tuple{Int,Int}`: A tuple of integers representing the range of the number of parameters of general_rate_equation to search over.
- `forward_model_selection::Bool`: A boolean indicating whether to use forward model selection (true) or reverse model selection (false).

# Returns nothing, but saves a csv file for each `num_params` with the results of the training for each combination of parameters tested and a csv file with test results for top 10% of the best results with each number of parameters tested.

"""
function data_driven_rate_equation_selection(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool,
)
    #check that range_number_params within bounds of minimal and maximal number of parameters
    @assert range_number_params[1] >=
            (1 + sum([occursin("K_a_", string(param_name)) for param_name in param_names]))
    @assert range_number_params[2] <= length(param_names)
    println("Past assertions")
    #generate param_removal_code_names by converting each mirror parameter for a and i into one name
    #(e.g. K_a_Metabolite1 and K_i_Metabolite1 into K_Metabolite1)
    param_removal_code_names = (
        [
            Symbol(replace(string(param_name), "_a" => "")) for
            param_name in param_names if !contains(string(param_name), "_i")
        ]...,
    )

    #generate all possible combination of parameter removal codes
    all_param_removal_codes = calculate_all_parameter_removal_codes(param_names)

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    if forward_model_selection
        num_param_range = (range_number_params[2]-1):-1:range_number_params[1]
        starting_param_removal_codes = [
            x for
            x in all_param_removal_codes if length(param_names) - num_alpha_params -
            sum(values(x[1:(end-num_alpha_params)]) .> 0) == range_number_params[2]
        ]
    elseif !forward_model_selection
        num_param_range = (range_number_params[1]+1):1:range_number_params[2]
        starting_param_removal_codes = [
            x for
            x in all_param_removal_codes if length(param_names) - num_alpha_params -
            sum(values(x[1:(end-num_alpha_params)]) .> 0) == range_number_params[1]
        ]
    end

    previous_param_removal_codes = starting_param_removal_codes
    println("About to start loop with num_params: $num_param_range")
    for num_params in num_param_range
        println("Running loop with num_params: $num_params")

        #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
        if forward_model_selection
            nt_param_removal_codes = forward_selection_next_param_removal_codes(
                all_param_removal_codes,
                previous_param_removal_codes,
                num_params,
                param_names,
                param_removal_code_names,
            )
        elseif !forward_model_selection
            nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                all_param_removal_codes,
                previous_param_removal_codes,
                num_params,
                param_names,
                param_removal_code_names,
            )
        end
        #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
        #TODO: change to pmap
        results_array = map(
            nt_param_removal_code -> train_rate_equation(
                general_rate_equation,
                data,
                metab_names,
                param_names;
                n_iter = 20,
                nt_param_removal_code = nt_param_removal_code,
            ),
            nt_param_removal_codes,
        )

        #convert results_array to DataFrame and save in csv file
        df_results = DataFrame(results_array)
        df_results.nt_param_removal_codes = nt_param_removal_codes
        df_results
        # CSV.write(
        #     "$(Dates.format(now(),"mmddyy"))_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
        #     df_results,
        # )
        #store top 10% for next loop as `previous_param_removal_codes`
        filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)
        previous_param_removal_codes = values.(df_results.nt_param_removal_codes)

        #calculate test loss for top 10% subsets for each `num_params`
        #TODO: loop over all figures and calculate test loss for each figure
        #TODO: consider looping over all figures and calculating test loss separately from train loss calculations
        #TODO: change to pmap
        test_loss = map(
            nt_fitted_params -> test_rate_equation(
                general_rate_equation,
                data,
                nt_fitted_params,
                metab_names,
                param_names,
            ),
            df_results.params,
        )
        #store rescaled results

    end
    println("Finished loop")

    #return train loss and params for all tested subsets, test loss for all tested subsets
end

"""Function to calculate loss for a given `rate_equation` and `nt_fitted_params` on `data` that was not used for training"""
function test_rate_equation(
    rate_equation::Function,
    data,
    nt_fitted_params::NamedTuple,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
)
    # Add a new column to data to assign an integer to each source/figure from publication
    data.fig_num = vcat(
        [
            i * ones(Int64, count(==(unique(data.source)[i]), data.source)) for
            i = 1:length(unique(data.source))
        ]...,
    )

    # Convert DF to NamedTuple for better type stability / speed
    rate_data_nt =
        Tables.columntable(data[.!isnan.(data.Rate), [:Rate, metab_names..., :fig_num]])

    # Make a vector containing indexes of points corresponding to each figure
    fig_point_indexes = [findall(data.fig_num .== i) for i in unique(data.fig_num)]
    fitted_params = values(nt_fitted_params)
    test_loss = loss_rate_equation(
        fitted_params,
        rate_equation::Function,
        rate_data_nt::NamedTuple,
        param_names::Tuple{Symbol,Vararg{Symbol}},
        fig_point_indexes::Vector{Vector{Int64}};
        rescale_params_from_0_10_scale = false,
        nt_param_removal_code = nothing,
    )
    return test_loss
end

"""Generate all possibles codes for ways that mirror params for a and i states of MWC enzyme can be removed from the rate equation"""
function calculate_all_parameter_removal_codes(param_names::Tuple{Symbol,Vararg{Symbol}})
    feasible_param_subset_codes = ()
    for param_name in param_names
        if param_name == :L
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif occursin("Vmax_a", string(param_name))
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2])
        elseif occursin("K_a", string(param_name))
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2, 3])
        elseif occursin("alpha", string(param_name))
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        end
    end
    return collect(Iterators.product(feasible_param_subset_codes...))
end

"""
Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_removal_code
"""
function param_subset_select(params, param_names, nt_param_removal_code)
    @assert length(params) == length(param_names)
    params_dict =
        Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

    for param_choice in keys(nt_param_removal_code)
        if startswith(string(param_choice), "L") && nt_param_removal_code[param_choice] == 1
            params_dict[:L] = 0.0
        elseif startswith(string(param_choice), "Vmax") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[:Vmax_i] = params_dict[:Vmax_a]
        elseif startswith(string(param_choice), "Vmax") &&
               nt_param_removal_code[param_choice] == 2
            global params_dict[:Vmax_i] = 0.0
        elseif startswith(string(param_choice), "K") &&
               nt_param_removal_code[param_choice] == 1
            K_i = Symbol("K_i_" * string(param_choice)[3:end])
            K_a = Symbol("K_a_" * string(param_choice)[3:end])
            params_dict[K_i] = params_dict[K_a]
        elseif startswith(string(param_choice), "K") &&
               nt_param_removal_code[param_choice] == 2
            K_a = Symbol("K_a_" * string(param_choice)[3:end])
            params_dict[K_a] = Inf
        elseif startswith(string(param_choice), "K") &&
               nt_param_removal_code[param_choice] == 3
            K_i = Symbol("K_i_" * string(param_choice)[3:end])
            params_dict[K_i] = Inf
        elseif startswith(string(param_choice), "alpha") &&
               nt_param_removal_code[param_choice] == 0
            params_dict[param_choice] = 0.0
        elseif startswith(string(param_choice), "alpha") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[param_choice] = 1.0
        end
    end

    new_params_sorted = [params_dict[param_name] for param_name in param_names]
    return new_params_sorted
end

"""
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
)

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        length(param_names) - num_alpha_params -
        sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params + 1 for
        param_removal_code in previous_param_removal_codes
    ])
    previous_param_subset_masks = unique([
        (
            mask = (
                (previous_param_removal_code[1:(end-num_alpha_params)] .== 0)...,
                zeros(Int64, num_alpha_params)...,
            ),
            non_zero_params = previous_param_removal_code .*
                              (previous_param_removal_code .!= 0),
        ) for previous_param_removal_code in previous_param_removal_codes
    ])

    #select all param_removal_codes that yield equations with `num_params` number of parameters
    all_param_codes_w_num_params = [
        param_removal_codes for param_removal_codes in all_param_removal_codes if (
            length(param_names) - num_alpha_params -
            sum(param_removal_codes[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                param_code_w_num_params .* previous_param_subset_mask.mask .+
                previous_param_subset_mask.non_zero_params for
                param_code_w_num_params in all_param_codes_w_num_params if (
                    length(param_names) - num_alpha_params - sum(
                        (param_code_w_num_params.*previous_param_subset_mask.mask.+previous_param_subset_mask.non_zero_params)[1:(end-num_alpha_params)] .>
                        0,
                    )
                ) == num_params
            ])...,
        )
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(param_removal_codes)]
    return nt_param_removal_codes
end

"""
Calculate `param_removal_codes` with `num_params` including zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params+1`
"""
function reverse_selection_next_param_removal_codes(
    all_param_removal_codes,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
)

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        length(param_names) - num_alpha_params -
        sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params - 1 for
        param_removal_code in previous_param_removal_codes
    ])
    previous_param_subset_masks = unique([
        (
            mask = [
                (previous_param_removal_code[1:(end-num_alpha_params)] .== 0)...,
                zeros(Int64, num_alpha_params)...,
            ],
            non_zero_params = previous_param_removal_code .*
                              (previous_param_removal_code .!= 0),
        ) for previous_param_removal_code in previous_param_removal_codes
    ])

    #select all codes that yield equations with `num_params` number of parameters
    all_param_codes_w_num_params = [
        param_removal_codes for param_removal_codes in all_param_removal_codes if (
            length(param_names) - num_alpha_params -
            sum(param_removal_codes[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                previous_param_subset_mask.non_zero_params .*
                (param_code_w_num_params .!= 0) for
                param_code_w_num_params in all_param_codes_w_num_params if (
                    length(param_names) - num_alpha_params - sum(
                        (previous_param_subset_mask.non_zero_params.*(param_code_w_num_params.!=0))[1:(end-num_alpha_params)] .>
                        0,
                    )
                ) == num_params
            ])...,
        )
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(param_removal_codes)]
    return nt_param_removal_codes
end
