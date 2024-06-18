using Dates, CSV, DataFrames, Distributed, HypothesisTests
include("rate_equation_fitting.jl")


function prepare_data(data::DataFrame, metab_names)

    # Check if the column source exists and add it if it doesn't
    if !hasproperty(data, :source)
           #Add source column that uniquely identifies a figure from publication
        data.source .= data.Article .* "_" .* data.Fig
    end

    # Remove Na's
    data = data[.!isnan.(data.Rate), [:Rate, metab_names..., :source]]

    #Only include Rate > 0 because otherwise log_ratio_predict_vs_data() will have to divide by 0
    filter!(row -> row.Rate != 0, data)

    # Check if all values in metab_names are columns in the data
    missing_columns = setdiff(metab_names, Symbol.(names(data)))
    @assert isempty(missing_columns) "The following metab columns are missing from the data: $(join(missing_columns, ", "))"
    
    return data
end



function data_driven_rate_equation_selection(
    general_rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    range_number_params::Tuple{Int,Int},
    forward_model_selection::Bool;
    n_reps_opt::Int = 20,
    maxiter_opt::Int = 50_000,
    model_selection_method = "denis",
    p_val_threshold = .4,
    save_train_results::Bool = false,
    enzyme_name::String = "Enzyme",
)
    
    data = prepare_data(data, metab_names)
    
    #generate param_removal_code_names by converting each mirror parameter for a and i into one name
    #(e.g. K_a_Metabolite1 and K_i_Metabolite1 into K_Metabolite1)
    param_removal_code_names = (
        [
            Symbol(replace(string(param_name), "_a_" => "_allo_")) for
            param_name in param_names if
            !contains(string(param_name), "_i") && param_name != :Vmax
        ]...,
    )
    
    #check that range_number_params within bounds of minimal and maximal number of parameters
    @assert range_number_params[1] >= length(param_names) - length(param_removal_code_names) "starting range_number_params cannot be below $(length(param_names) - length(param_removal_code_names))"
    @assert range_number_params[2] <= length(param_names) "ending range_number_params cannot be above $(length(param_names))"


    #generate all possible combination of parameter removal codes
    # param_subsets_per_n_params = calculate_all_parameter_removal_codes(param_names, range_number_params)
    all_param_removal_codes = calculate_all_parameter_removal_codes(param_names)

    if model_selection_method == "denis"
        results = fit_rate_equation_selection_denis(
            general_rate_equation,
            data,
            metab_names,
            param_names,
            param_removal_code_names, 
            range_number_params,
            forward_model_selection,
            n_reps_opt,
            maxiter_opt,
            all_param_removal_codes, 
            save_train_results, 
            enzyme_name
            )

        best_n_params = find_optimal_n_params(results.test_results, p_val_threshold)
        best_subset = get_nt_subset(results.test_results, best_n_params)
        println("Best subset")
        println(best_subset)

        # find best_subset row in train_results
        best_subset_row = filter(row -> row.nt_param_removal_codes == best_subset, results.train_results)
        println("best subset row")
        println(best_subset_row)


    elseif model_selection_method == "cv_denis"
        figs = unique(data.source) 
        results_figs_df = pmap(
            dropped_fig -> fit_rate_equation_selection_per_fig(
                general_rate_equation,
                data,
                metab_names,
                param_names,
                param_removal_code_names,
                range_number_params,
                forward_model_selection,
                n_reps_opt,
                maxiter_opt,
                all_param_removal_codes,
                dropped_fig
                ), 
            figs
        )
        train_results = [res.train_results for res in results_figs_df]
        test_results = [res.test_results for res in results_figs_df]
        combined_train_results = vcat(train_results...)
        combined_test_results = vcat(test_results...)
        results = (train_results =combined_train_results, test_results =combined_test_results )

        best_n_params = find_optimal_n_params(results.test_results, p_val_threshold)

        best_subset_row = train_and_choose_best_subset(
            general_rate_equation, 
            data,
            all_param_removal_codes, 
            best_n_params,
            metab_names, 
            param_names, 
            param_removal_code_names, 
            n_reps_opt, 
            maxiter_opt, 
            save_train_results, 
            enzyme_name
        )
        println("best subset row")
        println(best_subset_row)

    elseif model_selection_method == "cv_all_subsets"
        results =  fit_rate_equation_selection_all_subsets(
            general_rate_equation,
            data,
            all_param_removal_codes, 
            meta_names,
            param_names,
            param_removal_code_names,
            n_reps_opt, 
            maxiter_opt
            )

        best_n_params = find_optimal_n_params(results, p_val_threshold)

        best_subset_row = train_and_choose_best_subset(
            general_rate_equation, 
            data,
            all_param_removal_codes, 
            best_n_params,
            metab_names, 
            param_names, 
            param_removal_code_names, 
            n_reps_opt, 
            maxiter_opt, 
            save_train_results, 
            enzyme_name
        )
        println("best subset row")
        println(best_subset_row)

    end

    return (results = results, best_n_params = best_n_params, best_subset_row = best_subset_row)
end

function get_nt_subset(df, num)
    # Filter the DataFrame where n_params equals num
    filtered_df = filter(row -> row.num_params == num, df)

    return filtered_df.nt_param_removal_codes[1]

end

function find_optimal_n_params(df_results::DataFrame, p_value_threshold::Float64) :: Int
    # Group by number of parameters and calculate average test loss
    grouped = groupby(df_results, :num_params)
    avg_losses = combine(grouped, :test_loss => mean => :avg_test_loss)
    # Sort by number of parameters
    sort!(avg_losses, :num_params)
    println("Avg CV error for each n params:")
    println(avg_losses)
    # Find the row with the minimum average test loss
    idx_min_loss = argmin(avg_losses.avg_test_loss)
    n_param_minimal_loss = avg_losses[idx_min_loss, :num_params]
    losses_minimal_loss = filter(row -> row.num_params == n_param_minimal_loss, df_results).test_loss

    current_n_params = n_param_minimal_loss
    # Start checking from the model just below the minimal average loss model downwards
    for i in idx_min_loss-1:-1:1        
        current_n_params = avg_losses[i, :num_params]
        # Perform Wilcoxon signed-rank test on test losses
        losses_current = filter(row -> row.num_params == current_n_params, df_results).test_loss
        # compare with best n params: 
        test_result = SignedRankTest(losses_current, losses_minimal_loss)
        pval = pvalue(test_result)

        # If the difference is not significant, continue; else, stop and return last non-significant model's params
        if pval <= p_value_threshold
            current_n_params = avg_losses[i+1, :num_params]
            break  # Stop if a significant difference is found
        end
    end
    
    best_n_params = current_n_params

    return best_n_params
end



function fit_rate_equation_selection_denis(
        general_rate_equation::Function,
        data::DataFrame,
        metab_names::Tuple{Symbol,Vararg{Symbol}},
        param_names::Tuple{Symbol,Vararg{Symbol}},
        param_removal_code_names, 
        range_number_params::Tuple{Int,Int},
        forward_model_selection::Bool,
        n_repetiotions_opt::Int,
        maxiter_opt::Int,
        all_param_removal_codes, 
        save_train_results::Bool, 
        enzyme_name::String
        )

        num_alpha_params = count(occursin.("alpha", string.([param_names...])))

        if forward_model_selection
            num_param_range = (range_number_params[2]):-1:range_number_params[1]
        elseif !forward_model_selection
            num_param_range = (range_number_params[1]):1:range_number_params[2]
        end

        starting_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
            num_param_range[1],
            all_param_removal_codes,
            param_names,
            param_removal_code_names,
            num_alpha_params,
        )

        # starting_param_removal_codes = param_subsets_per_n_params[num_param_range[1]]
    
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
            if save_train_results
                CSV.write(
                    "$(Dates.format(now(),"mmddyy"))_$(enzyme_name)_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
                    df_results,
                )
            end   

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
            # previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
            nt_previous_param_removal_codes = [
                NamedTuple{param_removal_code_names}(x) for
                x in values.(df_results.nt_param_removal_codes)
            ]

            # save best subset for each `num_params` (afterwards loocv test loss will be calculated)
            best_nt_param_removal_code =
                df_results.nt_param_removal_codes[argmin(df_results.train_loss)]
            
            df_results = DataFrame(:num_params => [num_params], :nt_param_removal_codes => [best_nt_param_removal_code])
            df_test_results = vcat(df_test_results, df_results)
        end
       
        # calculate loocv test loss for top subsets:
        # Prepare the data for pmap
        subsets_to_fit = [(row.nt_param_removal_codes, removed_fig, row.num_params) for row in eachrow(df_test_results) for removed_fig in unique(data.source)]

        results = pmap(
            subset -> loocv_rate_equation(
                subset[2], #removed_fig
                general_rate_equation,
                data,
                metab_names,
                param_names;
                n_iter = n_repetiotions_opt,
                maxiter_opt = maxiter_opt,
                nt_param_removal_code = subset[1],
            ), 
        subsets_to_fit
        )
        # arrange test result ds
        result_dfs = DataFrame[]
        for (res, subset) in zip(results, subsets_to_fit)
            res_df = DataFrame([res])
            res_df[!, :nt_param_removal_codes] = [subset[1]]
            res_df[!, :num_params] = [subset[3]]
            push!(result_dfs, res_df)
        end

        df_test_results = vcat(result_dfs...)

        return (train_results = df_train_results, test_results = df_test_results)
    

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
    all_param_removal_codes,
    test_fig
    )

    train_data = data[data.source.!=test_fig, :]
    test_data = data[data.source.==test_fig, :]

    num_alpha_params = count(occursin.("alpha", string.([param_names...])))

    if forward_model_selection
        num_param_range = (range_number_params[2]):-1:range_number_params[1]
    elseif !forward_model_selection
        num_param_range = (range_number_params[1]):1:range_number_params[2]
    end

    starting_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
        num_param_range[1],
        all_param_removal_codes,
        param_names,
        param_removal_code_names,
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
                train_data,
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
        df_results.dropped_fig = fill(test_fig, nrow(df_results))
        df_results.nt_param_removal_codes = nt_param_removal_codes
        df_train_results = vcat(df_train_results, df_results)
        
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
        # previous_param_removal_codes = values.(df_results.nt_param_removal_codes)
        nt_previous_param_removal_codes = [
            NamedTuple{param_removal_code_names}(x) for
            x in values.(df_results.nt_param_removal_codes)
        ]

        # Save the best subset for each num_params. afterwards, test loss will be calculated using test_fig
        idx_min_row = argmin(df_results.train_loss)
        best_nt_param_removal_code = df_results[idx_min_row, :nt_param_removal_codes]
        best_subset_rescaled_params =  df_results[idx_min_row, :params]

        df_results = DataFrame(:num_params => num_params,
        :nt_param_removal_codes => best_nt_param_removal_code,
        :params => best_subset_rescaled_params)
        
        df_test_results = vcat(df_test_results, df_results)
    end
    
    # calculate test loss for top subsets:
    # Prepare the data for pmap
    subsets_to_test = [(row.params, row.nt_param_removal_codes,row.num_params) for row in eachrow(df_test_results)]

    test_results = pmap(
        best_subset_params -> test_rate_equation(
            general_rate_equation,
            test_data,
            best_subset_params[1], #rescaled params 
            metab_names, 
            param_names
        ), 
        subsets_to_test
    )

    result_dfs = DataFrame[]
    for (res, subset) in zip(test_results, subsets_to_test)
        res_df = DataFrame(
            test_loss = res,          
            num_params = subset[3],          
            nt_param_removal_code =subset[2], 
            test_fig =test_fig,
            params = subset[1]           
        )
        push!(result_dfs, res_df)
    end

    df_test_results = vcat(result_dfs...)
    # df_results =  DataFrame(
    #     test_loss = test_loss,          
    #     num_params = num_params,          
    #     nt_param_removal_code =best_nt_param_removal_code, 
    #     test_fig =test_fig,
    #     params = best_subset_rescaled_params           
    # )
        
    # df_test_results = vcat(df_test_results, df_results)
    

    return (train_results = df_train_results, test_results = df_test_results)

end

function fit_rate_equation_selection_all_subsets(
    general_rate_equation::Function,
    data::DataFrame,
    all_param_removal_codes, 
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    n_reps_opt::Int,
    maxiter_opt::Int,
    )

    # create param_subsets_per_n_params
    len_param_subset = length(first(all_param_removal_codes))
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))
    n = length(param_names)

    # keep for each number of params: all the subsets with this number
    param_subsets_per_n_params = Dict{Int, Vector{NTuple{len_param_subset, Int}}}()
    # for x in Iterators.take(all_param_removal_codes, 30000)
    for x in all_param_removal_codes
        n_param = n - num_alpha_params - sum(x[1:end-num_alpha_params] .> 0)
        if !haskey(param_subsets_per_n_params, n_param)
            param_subsets_per_n_params[n_param] = Vector{NTuple{len_param_subset, Int}}()
        end
        push!(param_subsets_per_n_params[n_param], x)
    end

    figs = unique(data.source)

    # Initialize an empty list for the combined results
    all_subsets_figs_to_fit = []
    lengths = []

    for (n_params, subsets) in param_subsets_per_n_params
        nt_param_subsets = [
            NamedTuple{param_removal_code_names}(x) for
            x in unique(param_removal_codes)
        ]

        # Create the product for this particular number of parameters
        temp_product = collect(Iterators.product(nt_param_subsets, figs))
        # Append the product to the main list
        append!(all_subsets_figs_to_fit, temp_product)
        # Record the length of the product
        push!(lengths, length(temp_product))
    end
    
    # Create the parameter mapping using the recorded lengths
    n_params_mapping = Int[]
    for (n_params, length) in zip(keys(param_subsets_per_n_params), lengths)
        append!(n_params_mapping, fill(n_params, length))
    end

    results_array = pmap(
        subset_fig_to_fit -> loocv_rate_equation(
            subset_fig_to_fit[2],
            general_rate_equation,
            data,
            metab_names,
            param_names;
            n_iter = n_reps_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = subset_fig_to_fit[1],
        ),
        all_subsets_figs_to_fit, 
    )

    df_results = DataFrame(results_array)
    df_results.num_params = n_params_mapping
    all_subsets  = [item[1] for item in all_subsets_figs_to_fit]
    df_results.nt_param_removal_codes = all_subsets
       
    return (train_test_results = df_results)

end


"function to calculate train loss without a figure and test loss on removed figure"
function loocv_rate_equation(
    fig,
    rate_equation::Function,
    data::DataFrame,
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}};
    n_iter = 20,
    maxiter_opt = 50_000,
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
        maxiter_opt = maxiter_opt,
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
        train_loss = train_res.train_loss,
        test_loss = test_loss,
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
        param_name_str = string(param_name)
        if param_name == :L
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(param_name_str, "Vmax_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2])
        elseif startswith(param_name_str, "K_a")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2,3])
        elseif startswith(param_name_str, "K_") &&
               !startswith(param_name_str, "K_i") &&
               !startswith(param_name_str, "K_a") &&
               length(split(param_name_str, "_")) == 2
               feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        elseif startswith(param_name_str, "K_") &&
               !startswith(param_name_str, "K_i") &&
               !startswith(param_name_str, "K_a") &&
               length(split(param_name_str, "_")) > 2
               feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1,2])
        elseif startswith(string(param_name), "alpha")
            feasible_param_subset_codes = (feasible_param_subset_codes..., [0, 1])
        end
    end
    return Iterators.product(feasible_param_subset_codes...)
end

"""Generate NamedTuple of codes for ways that params can be removed from the rate equation but still leave `num_params`"""
function calculate_all_parameter_removal_codes_w_num_params(
    num_params::Int,
    all_param_removal_codes,
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names::Tuple{Symbol,Vararg{Symbol}},
    num_alpha_params::Int,
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
    nt_param_removal_codes =
        [NamedTuple{param_removal_code_names}(x) for x in unique(codes_with_num_params)]
    return nt_param_removal_codes
end

"""
Function to convert parameter vector to vector where some params are equal to 0, Inf or each other based on nt_param_removal_code
"""
# function param_subset_select_denis(params, param_names, nt_param_removal_code)
#     @assert length(params) == length(param_names)
#     params_dict =
#         Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

#     for param_choice in keys(nt_param_removal_code)
#         if startswith(string(param_choice), "L") && nt_param_removal_code[param_choice] == 1
#             params_dict[:L] = 0.0
#         elseif startswith(string(param_choice), "Vmax") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[:Vmax_i] = params_dict[:Vmax_a]
#         elseif startswith(string(param_choice), "Vmax") &&
#                nt_param_removal_code[param_choice] == 2
#             global params_dict[:Vmax_i] = 0.0
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 1
#             K_i = Symbol("K_i_" * string(param_choice)[8:end])
#             K_a = Symbol("K_a_" * string(param_choice)[8:end])
#             params_dict[K_i] = params_dict[K_a]
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 2
#             K_a = Symbol("K_a_" * string(param_choice)[8:end])
#             params_dict[K_a] = Inf
#         elseif startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 3
#             K_i = Symbol("K_i_" * string(param_choice)[8:end])
#             params_dict[K_i] = Inf
#         elseif startswith(string(param_choice), "K_") &&
#                !startswith(string(param_choice), "K_allo") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[param_choice] = Inf
#         elseif startswith(string(param_choice), "K_") &&
#                !startswith(string(param_choice), "K_allo") &&
#                length(split(string(param_choice), "_")) > 2 &&
#                nt_param_removal_code[param_choice] == 2
#             params_dict[param_choice] =
#                 prod([
#                     params_dict[Symbol("K_" * string(metab))] for
#                     metab in split(string(param_choice), "_")[2:end]
#                 ])^(1 / (length(split(string(param_choice), "_")[2:end])))
#         elseif startswith(string(param_choice), "alpha") &&
#                nt_param_removal_code[param_choice] == 0
#             params_dict[param_choice] = 0.0
#         elseif startswith(string(param_choice), "alpha") &&
#                nt_param_removal_code[param_choice] == 1
#             params_dict[param_choice] = 1.0
#         end
#     end

#     new_params_sorted = [params_dict[param_name] for param_name in param_names]
#     return new_params_sorted
# end

function param_subset_select(params, param_names, nt_param_removal_code)
    @assert length(params) == length(param_names)
    params_dict =
        Dict(param_name => params[i] for (i, param_name) in enumerate(param_names))

    # for param_choice in keys(nt_param_removal_code)
    for (param, choice) in pairs(nt_param_removal_code)
        param_str = string(param)

        # handle K params
        if startswith(param_str, "K_allo")
            param_name = split(param_str, "K_allo_")[2]
            K_i = Symbol("K_i_" * param_name)
            K_a = Symbol("K_a_" * param_name)
      
            if choice > 0
                if choice == 1
                    params_dict[K_i] = params_dict[K_a] 

                elseif choice == 2
                    params_dict[K_a] = Inf

                elseif choice == 3
                    params_dict[K_i] = Inf
                end
            end

        elseif startswith(param_str, "K_") && !startswith(param_str, "K_allo")
            if choice == 1
                params_dict[Symbol(param_str)] = Inf
            elseif length(split(param_str, "_")) > 2 && choice == 2
                metabs = split(param_str, "_")[2:end]
                params_dict[Symbol(param_str)] = prod(params_dict[Symbol("K_" * metab)] for metab in metabs) ^ (1 / length(metabs))
            end
        
        elseif startswith(param_str, "alpha")
            if choice == 0
                params_dict[Symbol(param_str)] = 0.0
            elseif choice == 1
                params_dict[Symbol(param_str)] = 1.0
            end
            
        elseif startswith(param_str, "Vmax")
            if choice == 1
                params_dict[:Vmax_i] = params_dict[:Vmax_a]
            elseif choice == 2
                #TODO: check why it's appear with global in denis's code
                params_dict[:Vmax_i] = 0.0
            end

        elseif startswith(param_str, "L")
            if choice == 1
                params_dict[:L] = 0.0
            end

        end
    end
    new_params_sorted = [params_dict[param_name] for param_name in param_names]
    return new_params_sorted
end

"""
Calculate `nt_param_removal_codes` with `num_params` including non-zero term combinations for codes (excluding alpha terms) in each `nt_previous_param_removal_codes` that has `num_params-1`
"""
function forward_selection_next_param_removal_codes(
    nt_previous_param_removal_codes::Vector{T} where T<:NamedTuple,
    num_alpha_params::Int,
)
    feasible_param_subset_codes = []
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
    nt_previous_param_removal_codes::Vector{T} where T<:NamedTuple,
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


function train_and_choose_best_subset(
    general_rate_equation::Function,
    data::DataFrame,
    all_param_removal_codes, 
    best_n_params::Int, 
    metab_names::Tuple{Symbol,Vararg{Symbol}},
    param_names::Tuple{Symbol,Vararg{Symbol}},
    param_removal_code_names, 
    n_reps_opt::Int, 
    maxiter_opt::Int, 
    save_train_results::Bool, 
    enzyme_name::String
)
    num_alpha_params = count(occursin.("alpha", string.([param_names...])))

    nt_param_removal_codes = calculate_all_parameter_removal_codes_w_num_params(
        best_n_params,
        all_param_removal_codes,
        param_names,
        param_removal_code_names,
        num_alpha_params,
    )

    results_array = pmap(
        nt_param_removal_code -> train_rate_equation(
            general_rate_equation,
            data,
            metab_names,
            param_names;
            n_iter = n_reps_opt,
            maxiter_opt = maxiter_opt,
            nt_param_removal_code = nt_param_removal_code,
        ),
        nt_param_removal_codes,
    )

    #convert results_array to DataFrame
    df_results = DataFrame(results_array)
    df_results.num_params = fill(best_n_params, nrow(df_results))
    df_results.nt_param_removal_codes = nt_param_removal_codes

    # Optinally consider saving results to csv file for long running calculation of cluster
    if save_train_results
        CSV.write(
            "$(Dates.format(now(),"mmddyy"))_$(enzyme_name)_$(forward_model_selection ? "forward" : "reverse")_model_select_results_$(num_params)_num_params.csv",
            df_results,
        )
    end  

    best_param_subset = DataFrame(df_results[argmin(df_results.train_loss),:])
    println("Best subset: $(best_param_subset.nt_param_removal_codes)")

    return best_param_subset
end



