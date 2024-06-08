using Dates, CSV, DataFrames, Distributed

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


    #generate param_removal_code_names by converting each mirror parameter for a and i into one name
    #(e.g. K_a_Metabolite1 and K_i_Metabolite1 into K_allo_Metabolite1)
    param_removal_code_names = (
        [
            Symbol(replace(string(param_name), "_a_" => "_allo_", "Vmax_a" => "Vmax_allo")) for
            param_name in param_names if
            !contains(string(param_name), "_i") && param_name != :Vmax
        ]...,
    )

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    #check that range_number_params within bounds of minimal and maximal number of parameters
    @assert range_number_params[1] >= length(param_names) - length(param_removal_code_names) "starting range_number_params cannot be below $(length(param_names) - length(param_removal_code_names))"
    @assert range_number_params[2] <= length(param_names) "ending range_number_params cannot be above $(length(param_names))"

    if forward_model_selection
        num_param_range = (range_number_params[2]):-1:range_number_params[1]
    elseif !forward_model_selection
        num_param_range = (range_number_params[1]):1:range_number_params[2]
    end

    #calculate starting_param_removal_codes num_param_range[1] parameters
    all_param_removal_codes = calculate_all_parameter_removal_codes(param_names)
    starting_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
        num_param_range[1],
        all_param_removal_codes,
        param_names,
        num_alpha_params,
    )

    nt_param_removal_codes = starting_param_removal_codes
    nt_previous_param_removal_codes = similar(nt_param_removal_codes)
    println("About to start loop with num_params: $num_param_range")
    df_train_results = DataFrame()
    df_test_results = DataFrame()
    for num_params in num_param_range
        println("Running loop with num_params: $num_params")

        #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
        if num_params != num_param_range[1]
            if forward_model_selection
                nt_param_removal_codes = forward_selection_next_param_removal_codes(
                    nt_previous_param_removal_codes,
                    num_alpha_params,
                )
            elseif !forward_model_selection
                nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                    nt_previous_param_removal_codes,
                    num_alpha_params,
                )
            end
        end
        #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
        results_array = pmap(
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

        #convert results_array to DataFrame
        df_results = DataFrame(results_array)
        df_results.num_params = fill(num_params, nrow(df_results))
        df_results.nt_param_removal_codes = nt_param_removal_codes
        df_train_results = vcat(df_train_results, df_results)

        # Optinally consider saving results to csv file for long running calculation of cluster
        # CSV.write(
        #     "$(Dates.format(now(),"mmddyy"))_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
        #     df_results,
        # )

        #if all train_loss are Inf, then skip to next loop
        if all(df_results.train_loss .== Inf)
            nt_previous_param_removal_codes = [
                NamedTuple{param_removal_code_names}(x) for
                x in values.(df_results.nt_param_removal_codes)
            ]
            continue
        end

        #store top 10% for next loop as `previous_param_removal_codes`
        filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)
        nt_previous_param_removal_codes = [
            NamedTuple{param_removal_code_names}(x) for
            x in values.(df_results.nt_param_removal_codes)
        ]
        #calculate loocv test loss for top subset for each `num_params`
        best_nt_param_removal_code =
            df_results.nt_param_removal_codes[argmin(df_results.train_loss)]
        test_results = pmap(
            removed_fig -> loocv_rate_equation(
                removed_fig,
                general_rate_equation,
                data,
                metab_names,
                param_names;
                n_iter = 20,
                nt_param_removal_code = best_nt_param_removal_code,
            ),
            unique(data.source),
        )
        df_results = DataFrame(test_results)
        df_results.num_params = fill(num_params, nrow(df_results))
        df_results.nt_param_removal_codes =
            fill(best_nt_param_removal_code, nrow(df_results))
        df_test_results = vcat(df_test_results, df_results)
    end
    return (train_results = df_train_results, test_results = df_test_results)
end

"function to calculate train loss without a figure and test loss on removed figure"
function loocv_rate_equation(
    fig,
    rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    n_iter = 20,
    nt_param_removal_code = nothing,
)
    # Drop selected figure from data
    train_data = data[data.source.!=fig, :]
    test_data = data[data.source.==fig, :]
    # Calculate fit
    train_res = train_rate_equation(
        rate_equation,
        train_data,
        metab_names,
        param_names;
        n_iter = n_iter,
        nt_param_removal_code = nt_param_removal_code,
    )
    test_loss = test_rate_equation(
        rate_equation,
        test_data,
        train_res.params,
        metab_names,
        param_names,
    )
    return (
        dropped_fig = fig,
        train_loss_wo_fig = train_res.train_loss,
        test_loss_leftout_fig = test_loss,
        params = train_res.params,
    )
end

"""Function to calculate loss for a given `rate_equation` and `nt_fitted_params` on `data` that was not used for training"""
function test_rate_equation(
    rate_equation::Function,
    data,
    nt_fitted_params::NamedTuple,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
)
    filtered_data = data[.!isnan.(data.Rate), [:Rate, metab_names..., :source]]
    #Only include Rate > 0 because otherwise log_ratio_predict_vs_data() will have to divide by 0
    filter!(row -> row.Rate != 0, filtered_data)
    # Add a new column to data to assign an integer to each source/figure from publication
    filtered_data.fig_num = vcat(
        [
            i * ones(
                Int64,
                count(==(unique(filtered_data.source)[i]), filtered_data.source),
            ) for i = 1:length(unique(filtered_data.source))
        ]...,
    )
    # Add a column containing indexes of points corresponding to each figure
    fig_point_indexes =
        [findall(filtered_data.fig_num .== i) for i in unique(filtered_data.fig_num)]
    # Convert DF to NamedTuple for better type stability / speed
    rate_data_nt = Tables.columntable(filtered_data)

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

"""Generate all possibles codes for ways that params can be removed from the rate equation"""
function calculate_all_parameter_removal_codes(param_names::Tuple{Symbol,Vararg{Symbol}})
    feasible_param_subset_codes = ()
    for param_name in param_names
        if param_name == :L
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(string(param_name), "Vmax_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2])
        elseif startswith(string(param_name), "K_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2, 3])
        elseif startswith(string(param_name), "K_") &&
               !startswith(string(param_name), "K_i") &&
               !startswith(string(param_name), "K_a") &&
               length(split(string(param_name), "_")) == 2
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(string(param_name), "K_") &&
               !startswith(string(param_name), "K_i") &&
               !startswith(string(param_name), "K_a") &&
               length(split(string(param_name), "_")) > 2
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1, 2])
        elseif startswith(string(param_name), "alpha")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        end
    end
    # return collect(Iterators.product(feasible_param_subset_codes...))
    return Iterators.product(feasible_param_subset_codes...)
end

"""Generate codes for ways that params can be removed from the rate equation but still leave `num_params`"""
function calculate_all_parameter_removal_codes_w_num_params(
    num_params,
    all_param_removal_codes,
    param_names,
    num_alpha_params,
)
    codes_with_num_params = Tuple[]
    num_non_zero_in_each_code = Int[]
    for code in all_param_removal_codes
        sum_non_zero = 0
        for i = 1:(length(code)-num_alpha_params)
            if code[i] > 0
                sum_non_zero += 1
            end
        end
        push!(num_non_zero_in_each_code, sum_non_zero)
    end
    num_params_in_each_code =
        length(param_names) .- num_alpha_params .- num_non_zero_in_each_code
    for (i, code) in enumerate(all_param_removal_codes)
        if num_params_in_each_code[i] == num_params
            push!(codes_with_num_params, code)
        end
    end
    return codes_with_num_params
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
        elseif startswith(string(param_choice), "Vmax_allo") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[:Vmax_i] = params_dict[:Vmax_a]
        elseif startswith(string(param_choice), "Vmax_allo") &&
               nt_param_removal_code[param_choice] == 2
            global params_dict[:Vmax_i] = 0.0
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 1
            K_i = Symbol("K_i_" * string(param_choice)[8:end])
            K_a = Symbol("K_a_" * string(param_choice)[8:end])
            params_dict[K_i] = params_dict[K_a]
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 2
            K_a = Symbol("K_a_" * string(param_choice)[8:end])
            params_dict[K_a] = Inf
        elseif startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 3
            K_i = Symbol("K_i_" * string(param_choice)[8:end])
            params_dict[K_i] = Inf
        elseif startswith(string(param_choice), "K_") &&
               !startswith(string(param_choice), "K_allo") &&
               nt_param_removal_code[param_choice] == 1
            params_dict[param_choice] = Inf
        elseif startswith(string(param_choice), "K_") &&
               !startswith(string(param_choice), "K_allo") &&
               length(split(string(param_choice), "_")) > 2 &&
               nt_param_removal_code[param_choice] == 2
            params_dict[param_choice] =
                prod([
                    params_dict[Symbol("K_" * string(metab))] for
                    metab in split(string(param_choice), "_")[2:end]
                ])^(1 / (length(split(string(param_choice), "_")[2:end])))
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
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `nt_previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    nt_previous_param_removal_codes::Vector{T} where {T<:NamedTuple},
    num_alpha_params,
)
    param_removal_code_names = keys(nt_previous_param_removal_codes[1])
    next_param_removal_codes = Vector{Vector{Int}}()
    for previous_param_removal_code in nt_previous_param_removal_codes
        i_cut_off = length(previous_param_removal_code) - num_alpha_params
        for (i, code_element) in enumerate(previous_param_removal_code)
            if i <= i_cut_off && code_element == 0
                if param_removal_code_names[i] == :L
                    feasible_param_subset_codes = [1]
                elseif startswith(string(param_removal_code_names[i]), "Vmax_allo")
                    feasible_param_subset_codes = [1, 2]
                elseif startswith(string(param_removal_code_names[i]), "K_allo")
                    feasible_param_subset_codes = [1, 2, 3]
                elseif startswith(string(param_removal_code_names[i]), "K_") &&
                       !startswith(string(param_removal_code_names[i]), "K_allo") &&
                       length(split(string(param_removal_code_names[i]), "_")) == 2
                    feasible_param_subset_codes = [1]
                elseif startswith(string(param_removal_code_names[i]), "K_") &&
                       !startswith(string(param_removal_code_names[i]), "K_allo") &&
                       length(split(string(param_removal_code_names[i]), "_")) > 2
                    feasible_param_subset_codes = [1, 2]
                end
                for code_element in feasible_param_subset_codes
                    next_param_removal_code = collect(Int, previous_param_removal_code)
                    next_param_removal_code[i] = code_element
                    push!(next_param_removal_codes, next_param_removal_code)
                end
            end
        end
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(next_param_removal_codes)]
    return nt_param_removal_codes
end

"""
Use `nt_previous_param_removal_codes` to calculate `nt_next_param_removal_codes` that have one additional zero elements except for for elements <= `num_alpha_params` from the end
"""
function reverse_selection_next_param_removal_codes(
    nt_previous_param_removal_codes::Vector{T} where {T<:NamedTuple},
    num_alpha_params::Int,
)
    param_removal_code_names = keys(nt_previous_param_removal_codes[1])
    next_param_removal_codes = Vector{Vector{Int}}()
    for previous_param_removal_code in nt_previous_param_removal_codes
        i_cut_off = length(previous_param_removal_code) - num_alpha_params
        for (i, code_element) in enumerate(previous_param_removal_code)
            if i <= i_cut_off && code_element != 0
                next_param_removal_code = collect(Int, previous_param_removal_code)
                next_param_removal_code[i] = 0
                push!(next_param_removal_codes, next_param_removal_code)
            end
        end
    end
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(next_param_removal_codes)]
    return nt_param_removal_codes
end
