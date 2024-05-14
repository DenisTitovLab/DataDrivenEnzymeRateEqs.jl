using Dates, CSV, DataFrames, Distributed
include("rate_equation_fitting.jl")


function prepare_data(data::DataFrame)

    # Check if the column source exists and add it if it doesn't
    if !hasproperty(data, :source)
           #Add source column that uniquely identifies a figure from publication
        data.source .= data.Article .* "_" .* data.Fig
    end

    # Remove Na's
    data = data[.!isnan.(data.Rate), [:Rate, metab_names..., :source]]

    #Only include Rate > 0 because otherwise log_ratio_predict_vs_data() will have to divide by 0
    filter!(row -> row.Rate != 0, data)

    return data
end



function data_driven_rate_equation_selection(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool,
    n_repetiotions_opt::Int,
    maxiter_opt::Int;
    model_selection_method = "denis",
)
    
    data = prepare_data(data)
    
    #check that range_number_params within bounds of minimal and maximal number of parameters
    @assert range_number_params[1] >=
            (1 + sum([occursin("K_a_", string(param_name)) for param_name in param_names]))
    @assert range_number_params[2] <= length(param_names)

    
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

    # keep for each number of params: all the subsets with this number
    param_subsets_tuple = [(length(param_names) - num_alpha_params -  sum(values(x[1:(end-num_alpha_params)]) .> 0) , values(x)) 
       for x in all_param_removal_codes]
    param_subsets_per_n_params = Dict{Int, Vector}()
    for (key, value) in param_subsets_tuple
        if haskey(param_subsets_per_n_params, key)
            push!(param_subsets_per_n_params[key], value)
        else
            param_subsets_per_n_params[key] = [value]
        end
    end

    if model_selection_method == "denis"
        results = fit_rate_equation_selection_per_fig(
            general_rate_equation,
            data,
            metab_names,
            param_names,
            param_removal_code_names, 
            range_number_params,
            forward_model_selection,
            n_repetiotions_opt,
            maxiter_opt,
            param_subsets_per_n_params,
            all_param_removal_codes;
            dropped_fig = nothing
            )
    end
    
end


function fit_rate_equation_selection_per_fig(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}},
        param_removal_code_names, 
        range_number_params::Tuple{Int,Int},
        forward_model_selection::Bool,
        n_repetiotions_opt::Int,
        maxiter_opt::Int,
        param_subsets_per_n_params,
        all_param_removal_codes;
        dropped_fig = nothing
        )


        if forward_model_selection
            num_param_range = (range_number_params[2]):-1:range_number_params[1]
            starting_param_removal_codes = param_subsets_per_n_params[range_number_params[2]]
        elseif !forward_model_selection
            num_param_range = (range_number_params[1]):1:range_number_params[2]
            starting_param_removal_codes = param_subsets_per_n_params[range_number_params[1]]
        end
    
        previous_param_removal_codes = starting_param_removal_codes
        println("About to start loop with num_params: $num_param_range")
        
        df_train_results = DataFrame()
        df_test_results = DataFrame()
        for num_params in num_param_range
            println("Running loop with num_params: $num_params")
    
            #calculate param_removal_codes for `num_params` given `all_param_removal_codes` and fixed params from previous `num_params`
            if forward_model_selection
                nt_param_removal_codes = forward_selection_next_param_removal_codes(
                    param_subsets_per_n_params,
                    previous_param_removal_codes,
                    num_params,
                    param_names,
                    param_removal_code_names,
                )
            elseif !forward_model_selection
                nt_param_removal_codes = reverse_selection_next_param_removal_codes(
                    param_subsets_per_n_params,
                    previous_param_removal_codes,
                    num_params,
                    param_names,
                    param_removal_code_names,
                )
            end
    
            #pmap over nt_param_removal_codes for a given `num_params` return rescaled and nt_param_subset added
            results_array = map(
                nt_param_removal_code -> train_rate_equation(
                    general_rate_equation,
                    data,
                    metab_names,
                    param_names;
                    n_iter = n_repetiotions_opt,
                    maxiter_opt = maxiter_opt,
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
            #store top 10% for next loop as `previous_param_removal_codes`
            filter!(row -> row.train_loss < 1.1 * minimum(df_results.train_loss), df_results)
            previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
    
            #calculate loocv test loss for top subset for each `num_params`
            #TODO: change to pmap
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

            if dropped_fig !== nothing
                df_test_results.dropped_fig = fill(dropped_fig, nrow(df_test_results))
                df_train_results.dropped_fig = fill(dropped_fig, nrow(df_train_results))
            end
    
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

function param_subset_select_may(params, param_names, nt_param_removal_code)
    @assert length(params) == length(param_names)
    params_dict =
        Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

    # for param_choice in keys(nt_param_removal_code)
    for (name, choice) in pairs(nt_param_removal_code)
        name_str = string(name)
        choice_str = string(choice)

        # handle K params
        if startswith(uppercase(name_str), "K")
            K_a = replace(name_str, "K_" => "K_a_")
            K_i = replace(name_str, "K_" => "K_i_")
      
            if choice > 0
                if choice == 1
                    params_dict[Symbol(K_i)] = params_dict[Symbol(K_a)] 

                elseif choice == 2
                    params_dict[Symbol(K_a)] = Inf

                elseif choice == 3
                    params_dict[Symbol(K_i)] = Inf
                end
            end
        
        elseif startswith(name_str, "alpha")
            if choice == 0
                params_dict[Symbol(name_str)] = 0.0
            elseif choice == 1
                params_dict[Symbol(name_str)] = 1.0
            end
            
        elseif name_str == "Vmax"
            if choice == 1
                params_dict[Symbol(name_str , "_i")] = 1.0
            elseif choice == 2
                #TODO: check why it's appear with global in denis's code
                params_dict[Symbol(name_str)] = 0.0
            end

        elseif name_str == "L"
            if choice == 1
                params_dict[Symbol(name_str)] = 0.0
            end

        end
    end
    new_params_sorted = [params_dict[param_name] for param_name in param_names]
    return new_params_sorted
end
"""
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    param_subsets_per_n_params,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
    )

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params + 1
        ) || (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params
        ) for param_removal_code in previous_param_removal_codes
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
    all_param_codes_w_num_params = param_subsets_per_n_params[num_params]

    # #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                param_code_w_num_params .* previous_param_subset_mask.mask .+
                previous_param_subset_mask.non_zero_params for
                param_code_w_num_params in all_param_codes_w_num_params
            ])...,
        )
    end
    nt_param_removal_codes = [
        NamedTuple{param_removal_code_names}(x) for
        x in unique(param_removal_codes) if (
            length(param_names) - num_alpha_params - sum(x[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    return nt_param_removal_codes
end

"""
Calculate `param_removal_codes` with `num_params` including zero term combinations for codes (excluding alpha terms) in each `previous_param_removal_codes` that has `num_params+1`
"""
function reverse_selection_next_param_removal_codes(
    param_subsets_per_n_params,
    previous_param_removal_codes,
    num_params,
    param_names,
    param_removal_code_names,
)

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    @assert all([
        (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params - 1
        ) || (
            length(param_names) - num_alpha_params -
            sum(param_removal_code[1:(end-num_alpha_params)] .> 0) == num_params
        ) for param_removal_code in previous_param_removal_codes
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
    all_param_codes_w_num_params = param_subsets_per_n_params[num_params]

    #choose param_removal_codes with n_removed_params number of parameters removed that also contain non-zero elements from previous_param_removal_codes
    param_removal_codes = []
    for previous_param_subset_mask in previous_param_subset_masks
        push!(
            param_removal_codes,
            unique([
                previous_param_subset_mask.non_zero_params .*
                (param_code_w_num_params .!= 0) for
                param_code_w_num_params in all_param_codes_w_num_params
            ])...,
        )
    end
    nt_param_removal_codes = [
        NamedTuple{param_removal_code_names}(x) for
        x in unique(param_removal_codes) if (
            length(param_names) - num_alpha_params - sum(x[1:(end-num_alpha_params)] .> 0)
        ) == num_params
    ]
    return nt_param_removal_codes
end